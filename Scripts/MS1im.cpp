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
struct scantable {
  std::vector< float > mzs; //mass to charge ratios from scan
  std::vector< float > rts; //mass to charge ratios from scan
  std::vector< float > dts; //mass to charge ratios from scan
  std::vector< float > its; //intensity values from scan
  std::vector< std::string > isos; //isotopes any isotopes less than mztol
  std::vector< std::string > formulas; //formula (if any) for closest mass
  int tid; //target feature table row ID
  int scanid;
  int scanrows;
  int originalid;
  float tmz;
  float trt;
  float tdt;
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

scantable process_scan_ft(
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
  struct scantable t;

  float min_dmz = 10000.0; // store the closest match to the target mz
  int min_dmz_id = 0;      // store the index of closest match
  std::vector<float> dmzs;
  for(int i = 0; i < scan.nrow(); i++)
  {
    float mz = scan(i, 1);
    float it = scan(i, 2);
    float rt = scan(i, 0);
    float dt = scan.ncol() == 4 ? scan(i, 3) : 0.0;
    // get everything that's within the mz window and has a non zero intensity
    if( (tmz + mzWindowLow < mz) && (mz < tmz + mzWindowHigh) && (0 < it) )
    {
      t.mzs.push_back( mz );
      t.its.push_back( it );
      t.rts.push_back( rt );
      t.dts.push_back( dt );
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
  // struct scantable emptytable;
  // return min_dmz < mztol ? t : emptytable;
  return t;
}

void print_tables(
    std::vector< scantable > tables,
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

      // out << std::to_string(tables[i].tid) << ",";
      // out << std::to_string(tables[i].scanid) << ",";
      // out << std::to_string(tables[i].scanrows) << ",";
      // out << std::to_string(tables[i].originalid) << ",";
      // out << tos(tables[i].mzs[j],5) << ",";
      // out << tos(tables[i].tmz,4) << ",";
      // out << tos(tables[i].rts[j],5) << ",";
      // out << tos(tables[i].trt,4) << ",";
      // out << tos(tables[i].dts[j],5) << ",";
      // out << tos(tables[i].tdt,4) << ",";
      // out << tos(tables[i].its[j],3) << ",";
      // out << tables[i].isos[j]       << ",";
      // out << tables[i].formulas[j]   << ",";
      // out <<  ""  << ",";
      // out << filenameSource << std::endl;
    }
  }
  out.close();
}


std::vector<int> get_closest_scan_ids(
    Rcpp::List scans, 
    Rcpp::NumericMatrix ft,
    float rttol,
    float dttol,
    bool hasdt
  ){

  std::vector<float> closest_scan_rts;
  std::vector<float> closest_scan_dts;
  std::vector<int> closest_scan_ids;
  for (int i = 0; i < ft.nrow(); i++) {
    closest_scan_rts.push_back(0.0); // store the closest rt
    closest_scan_dts.push_back(0.0); // store the closest dt
    closest_scan_ids.push_back(0); // store the '' index
  }

  int ftidx = 0;
  float ep = 0.0001; //account for fp error
  // traverse through every scan in the dataset
  for(int scanid = 0; scanid < scans.size(); scanid++){
    Rcpp::NumericMatrix scan = scans[scanid];
    float scan_rt = scan(0,0);
    float scan_dt = hasdt ? scan(0,3) : 0.0;

    // go through feature table within retention time tolerance
    for (int i = ftidx; i < ft.nrow(); i++) {
      float trt = ft(i, 1);
      float tdt = hasdt ? ft(i, 3) : 0.0;

      // possibly skip to (next scan | FT row)
      // This is so that we can skip over rows in either dataset
      if( scan_rt + rttol < trt ){ break; } // scan's rt isn't high enough
      if( trt < scan_rt - rttol ){ ftidx++; continue; } //scan's rt is too high, try next feature

      // get the distance to current retention time and closest rt (so far)
      float drt_scan = std::abs(trt - scan_rt);
      float drt_best = std::abs(trt - closest_scan_rts[i]);
      // enter if at new best rt
      if( drt_scan + ep < drt_best ){
        closest_scan_rts[i] = scan_rt;
        // set closest if it's not a drift time id
        if(!hasdt){
          closest_scan_ids[i] = scanid;
        }
        // std::cout << ft(i, 2) << ",  " << scanid << ",  " << tos(closest_scan_rts[i],4) << ",  " << tos(closest_scan_dts[i],4) << std::endl;
      }

      // enter if at current closest rt (within epsilon)
      if( drt_scan < drt_best + ep && hasdt ){
        float ddt_scan = std::abs(tdt - scan_dt);
        float ddt_best = std::abs(tdt - closest_scan_dts[i]);
        // enter if at current closest dt
        if( ddt_scan < ddt_best && ddt_scan < dttol ){
          closest_scan_dts[i] = scan_dt;
          closest_scan_ids[i] = scanid;
        }
      }
    }
  }

  // for (int i = 0; i < ft.nrow(); i++) {
  //   std::cout << ft(i, 2) << ",  " << closest_scan_ids[i] << ",  " << tos(closest_scan_rts[i],4) << ",  " << tos(closest_scan_dts[i],4) << std::endl;
  // }

  return closest_scan_ids;
}

std::vector<scantable> get_MS1_tables(
    Rcpp::List Iso,
    Rcpp::List scans,
    std::vector<int> closest_scan_ids,
    std::vector<int> mzrtdt_ids,
    Rcpp::NumericMatrix ft,
    std::vector<std::string> ft_formulas,
    float mzWindowLow,
    float mzWindowHigh,
    float mztol,
    float rttol,
    float dttol,
    bool hasdt
  ){

  std::vector<float> rm = Iso[0]; // relative masses
  std::vector<std::string> sym = Iso[1]; // symbols (ex: 13C1)
  std::vector<float> comp = Iso[2]; // composition/abundance, % found from generic sample ex:1.08%

  std::vector<scantable> MS1tables;
  for (int i = 0; i < ft.nrow(); i++) {
    int scanid = closest_scan_ids[i];
    Rcpp::NumericMatrix scan = scans[scanid];
    float scan_rt = scan(0,0);
    float scan_dt = hasdt ? scan(0,3) : 0.0;

    float tmz = ft(i, 0);
    float trt = ft(i, 1);
    int   tid = ft(i, 2);
    float tdt = ft.ncol() == 4 ? ft(i, 3) : 0.0;
    std::string tformula = ft_formulas[i];

    // rt (and dt) is within the tolerance, create an MS1 table
    // also support datasets which don't have drift time (non Ion Mobility)
    float rtbool = std::abs(trt - scan_rt) < rttol;
    float dtbool = hasdt ? std::abs(tdt - scan_dt) < dttol : true;

    // std::cout << "Processing feature " << i << " " << std::abs(trt - scan_rt) + std::abs(tdt - scan_dt) << std::endl;

    if( rtbool && dtbool){
      struct scantable table_i = process_scan_ft(rm, sym, comp, scan, tmz
                                  , tformula, mzWindowLow, mzWindowHigh, mztol);
      table_i.tmz = tmz;
      table_i.trt = trt;
      table_i.tdt = tdt;
      table_i.tid = tid;
      table_i.scanid = scanid;
      table_i.scanrows = scan.nrow();
      table_i.originalid = mzrtdt_ids[i];
      MS1tables.push_back( table_i );
    }
  }

  std::sort( begin(MS1tables), end(MS1tables), [&](scantable a, scantable b){ 
      return a.originalid < b.originalid;
    });

  return MS1tables;
}

// [[Rcpp::export]]
void save_MS1s_to_file(
    Rcpp::List Iso,
    Rcpp::List scans, 
    Rcpp::List FT,
    std::vector<int> mzrtdt_ids,
    float mzWindowLow,
    float mzWindowHigh,
    float mztol,
    float rttol,
    float dttol,
    bool hasdt,
    std::string filenameOutput,
    std::string filenameSource
  ){
  /*
    Inputs:
    scans - scan data in data frames sorted by RT and then DT
    FT    - feature table (data matrix, formulas)
    Iso   - isotopes, relative masses, symbols, compositions
    mzrtdt_ids - the indices which correspond to the feature table
    mztol - mass to charge tolerance (if there's not an mz within this, skip the feature table)
    rttol - retention time tolerance
    dttol - drift time tolerance
    filenameOutput - the destination where data will be saved
    filenameSource - the sourse where the scan data came from
    mzWindowLow  - how much lower than target mz can be considered
    mzWindowHigh - how much higher than target mz can be considered
    hasdt - allow this to proces data with(out) drift time
  */

  // unpack variables
  Rcpp::NumericMatrix ft = FT[0]; // mz, rt, id, dt(?)
  std::vector<std::string> ft_formulas = FT[1]; // formulas from feature table

  // Get the indices for which scan is closest in drt + ddt to the feature rows
  // This represents where the sum of the retention and drift time differences is minimized
  std::vector<int> closest_scan_ids = get_closest_scan_ids(scans, ft, rttol, dttol, hasdt);

  // once we have the closest scans,
  //   we can just directly access and process these from the scan data
  std::vector<scantable> MS1tables = get_MS1_tables(Iso, scans, closest_scan_ids, mzrtdt_ids,
    ft, ft_formulas, mzWindowLow, mzWindowHigh, mztol, rttol, dttol, hasdt);

  print_tables(MS1tables, filenameOutput, filenameSource);
}
