#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

// #include <algorithm>

#include <boost/algorithm/string/join.hpp>
#include <iterator>
#include <functional>
#include <Rcpp.h>

// mz, intensities, isotope strings, formulas, ft id
struct datatable {
  std::vector< float > mzs; //mass to charge ratios from scan
  std::vector< float > its; //intensity values from scan
  std::vector< std::string > isos; //isotopes any isotopes less than mztol
  std::vector< std::string > formulas; //formula (if any) for closest mass
  int tid; //target feature table row ID
};

// convert float to string with target precision (number of decimals)
std::string tos(float number, int precision){
  std::stringstream stream;
  stream << std::fixed << std::setprecision(precision) << number;
  std::string s = stream.str();
  return s;
}

// construct the isotope string ex: 33S1(0.79%;0.00099Da)
std::string isostring(std::string sym, float comp, float ddmz){
  std::stringstream stream;
  stream << sym << "(";
  stream << tos(comp,2) << "%;";
  stream << tos(ddmz,5) << ")";
  std::string s = stream.str();
  return s;
}

// Example output: 33S1(0.79%;0.00099Da)_15N1(0.37%;0.00136Da)_13C1(1.08%;0.00496Da)
std::string get_iso_string(
    std::vector<float> rm,
    std::vector<std::string> sym,
    std::vector<float> comp,
    float dmz,
    float mztol
  ){
  // The goal is to return a string of possible isotopes sorted by distance.
  std::vector<float> drms; //stores (rm - (tmz - mz))
  std::vector<int> ids; //ids from iso table, used to sort

  // store all of the distances to the relative masses for each isotope
  for(int i = 0; i < rm.size(); i++){
    drms.push_back( std::abs( rm[i] - dmz ) );
    ids.push_back( i );
  }

  //sort with ascending relative mass distance
  std::sort( begin(ids), end(ids), [&](int a, int b){ return drms[a] < drms[b]; } );
  std::vector<std::string> isos;
  // get each of the relative mass distances less than the tolerance
  for(int i = 0; i < ids.size(); i++){
    float drm = drms[ ids[i] ];
    if( mztol < drm ) break;
    isos.push_back( isostring(sym[ids[i]], comp[ids[i]], drm) );
  }
  std::string s = boost::algorithm::join(isos, "_"); // may need to "(apt|dnf) install boost-devel"
  return s;
}

datatable process_scan_ft(
    std::vector<float> rm,
    std::vector<std::string> sym,
    std::vector<float> comp,
    Rcpp::NumericMatrix scan,
    float tmz,
    std::string tformula,
    float mzWindowLow,
    float mzWindowHigh,
    float mztol
  ){
  // create structure to store data
  struct datatable t;

  float min_dmz = 10000.0; // store the closest match to the target mz
  int min_dmz_id = 0;      // store the index of closest match
  std::vector<float> dmzs;
  for(int i = 0; i < scan.nrow(); i++)
  {
    float mz = scan(i, 1);
    float it = scan(i, 2);
    // get everything that's within the mz window and has a non zero intensity
    if( (mz + mzWindowLow < tmz) && (tmz < mz + mzWindowHigh) && (0 < it) )
    {
      t.mzs.push_back( mz );
      t.its.push_back( it );
      dmzs.push_back( std::abs(tmz - mz) );
      // get the minimum id to set "M"
      if( dmzs.back() < min_dmz ){
        min_dmz    = dmzs.back();
        min_dmz_id = dmzs.size()-1;
      }
    }
  }

  // get the possible isotope strings for all after "M"
  for(int i = 0; i < t.mzs.size(); i++){
    std::string s = "";
    std::string f = "";
    if(i == min_dmz_id) {
      s = "M";
      f = tformula;
    }
    if(min_dmz_id < i) s = get_iso_string(rm, sym, comp, dmzs[i], mztol);
    t.isos.push_back(s);
    t.formulas.push_back(f);
  }

  // it's possible that "M" isn't close enough, so skip reporting this feature table 
  struct datatable emptytable;
  return min_dmz < mztol ? t : emptytable;
}

void print_tables(
    std::vector< datatable > tables,
    std::string filenameOutput,
    std::string filenameSource
  ){
  std::fstream out;
  out.open(filenameOutput, std::ios::app);

  for(int i = 0; i < tables.size(); i++)
  {
    for(int j = 0; j < tables[i].mzs.size(); j++){
      // Feature ID, m/z, Intensity, Isotope, Formula, PredAbundance, File
      out << std::to_string(tables[i].tid) << ",";
      out << tos(tables[i].mzs[j],5) << ",";
      out << tos(tables[i].its[j],0) << ",";
      out << tables[i].isos[j]       << ",";
      out << tables[i].formulas[j]   << ",";
      out <<  ""  << ",";
      out << filenameSource << std::endl;
    }
  }
  out.close();
}

// [[Rcpp::export]]
void save_MS1s_to_file(
    Rcpp::List scans, 
    Rcpp::List FT,
    Rcpp::List Iso,
    float rtw,
    float dtw,
    std::string filenameOutput,
    std::string filenameSource,
    float mzWindowLow,
    float mzWindowHigh,
    float mztol,
    bool hasdt
  ){
  /*
    Inputs:
    scans - scan data in data frames sorted by RT and then DT
    FT    - feature table (data matrix, formulas)
    Iso   - isotopes, relative masses, symbols, compositions
    rtw   - retention time window
    dtw   - drift time window
    filenameOutput - the destination where data will be saved
    filenameSource - the sourse where the scan data came from
    mzWindowLow  - how much lower than target mz can be considered
    mzWindowHigh - how much higher than target mz can be considered
    mztol - mass to charge tolerance (if there's not an mz within this, skip the feature table)
    hasdt - allow this to proces data with(out) drift time
  */

  // unpack variables
  Rcpp::NumericMatrix ft = FT[0]; // mz, rt, id, dt(?)
  std::vector<std::string> ft_formula = FT[1]; // formulas from feature table
  std::vector<float> rm = Iso[0]; // relative masses
  std::vector<std::string> sym = Iso[1]; // symbols (ex: 13C1)
  std::vector<float> comp = Iso[2]; // composition/abundance, % found from generic sample ex:1.08%

  int ftidx = 0;
  std::vector<datatable> tables;
  for(int scanid = 0; scanid < scans.size(); scanid++){
    // traverse through every scan in the dataset
    Rcpp::NumericMatrix scan = scans[scanid];
    float scan_rt = scan(0,0);

    // go through a matching window where feature table retention time is close
    for (int i = ftidx; i < ft.nrow(); i++) {
      float trt = ft(i, 1);

      // possibly skip to (next scan | FT row)
      if( scan_rt + rtw < trt ){ break; }
      if( trt < scan_rt - rtw ){ ftidx++; continue; }

      float tmz = ft(i, 0);
      float tid = ft(i, 2);
      float tdt = ft.ncol() == 4 ? ft(i, 3) : 0.0;
      bool dtbool = hasdt ? ((scan(0,3) - dtw < tdt) && (tdt < scan(0,3) + dtw)) : true;
      std::string tformula = ft_formula[i];
      // rt (and dt) is in the window, create an MS1 table
      if( (scan_rt - rtw < trt) && (trt < scan_rt + rtw) && dtbool){
        struct datatable table_i = process_scan_ft(rm, sym, comp, scan, tmz
                                    , tformula, mzWindowLow, mzWindowHigh, mztol);
        table_i.tid = tid;
        tables.push_back( table_i );
        ftidx++;
      }
    }
  }

  print_tables(tables, filenameOutput, filenameSource);
}
