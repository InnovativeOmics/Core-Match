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

struct datatable {
  std::vector< float > mzs;
  std::vector< float > its;
  std::vector< std::string > isos;
  std::vector< std::string > formulas;
  int tid;
};

// float to string with precision
std::string tos(float number, int precision){
  std::stringstream stream;
  stream << std::fixed << std::setprecision(precision) << number;
  std::string s = stream.str();
  return s;
}

// float to string with precision
std::string isostring(std::string sym, float comp, float ddmz){
  std::stringstream stream;
  stream << sym << "(";
  stream << tos(comp,2) << "%;";
  stream << tos(ddmz,5) << ")";
  std::string s = stream.str();
  return s;
}

std::string get_iso_string(
    std::vector<float> rm,
    std::vector<std::string> sym,
    std::vector<float> comp,
    float dmz,
    float mztol
  ){
  // The goal is to return a string of possible isotopes sorted by distance.
  // Example: 33S1(0.79%;0.00099Da)_15N1(0.37%;0.00136Da)_13C1(1.08%;0.00496Da)
  std::vector<float> ddmzs;
  std::vector<int> ids; //ids from iso table

  for(int i = 0; i < rm.size(); i++){
    ddmzs.push_back( std::abs( rm[i] - dmz ) );
    ids.push_back( i );
  }

  std::sort( begin(ids), end(ids), [&](int a, int b){ return ddmzs[a] < ddmzs[b]; } );
  std::vector<std::string> isos;
  for(int i = 0; i < ids.size(); i++){
    float ddmz = ddmzs[ ids[i] ];
    if( mztol < ddmz ) break;
    isos.push_back( isostring(sym[ids[i]], comp[i], ddmz) );
  }
  std::string s = boost::algorithm::join(isos, "_");
  return s;
}

datatable process_scan_ft(
    std::vector<float> rm,
    std::vector<std::string> sym,
    std::vector<float> comp,
    Rcpp::NumericMatrix scan,
    float tmz,
    std::string tformula
  ){

  struct datatable t;
  float mztol = 0.01;
  float mzlow = 1.0;
  float mzhigh = 5.0;
  float min_dmz = 10.0;
  int min_dmz_id = 0;
  std::vector<float> dmzs;

  for(int i = 0; i < scan.nrow(); i++)
  {
    float mz = scan(i, 1);
    float it = scan(i, 2);
    if( (mz - mzlow < tmz) && (tmz < mz + mzhigh) && (0 < it) )
    {
      t.mzs.push_back( mz );
      t.its.push_back( it );
      dmzs.push_back( std::abs(tmz - mz) );
      if( dmzs.back() < min_dmz ){
        min_dmz    = dmzs.back();
        min_dmz_id = dmzs.size()-1;
      }
    }
  }

  if(mztol < min_dmz) return t;

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

  return t;
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
// Feature, m/z, Intensity, Isotope, Formula, PredAbundance, File
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
void test_ISO(
    Rcpp::List scans, 
    Rcpp::List Iso,
    Rcpp::List FT,
    float rtw,
    float dtw,
    std::string filenameOutput,
    std::string filenameSource
  ){
  /*
    Inputs:
    scans - scan data in data frames sorted by RT and then DT
    Iso   - isotopes, relative masses, symbols, compositions
    FT    - feature table
    rtw   - retention time window
  */

  int ftidx = 0;
  Rcpp::NumericMatrix ft = FT[0]; // mz, rt, id, dt
  Rcpp::StringVector ft_formula = FT[1];

  std::vector<float> rm = Iso[0];
  std::vector<std::string> sym = Iso[1];
  std::vector<float> comp = Iso[2];

  // print_ions(rm, sym, comp);
  std::vector<datatable> tables;

  for(int scanid = 0; scanid < scans.size(); scanid++){
    Rcpp::NumericMatrix scan = scans[scanid];

    for (int i = ftidx; i < ft.nrow(); i++) {
      float tmz = ft(i, 0);
      float trt = ft(i, 1);
      float tid = ft(i, 2);
      float tdt = ft(i, 3);
      std::string tformula = "NA";
      if( scan(0,0) + rtw < trt ){ break; }
      if( trt < scan(0,0) - rtw ){ ftidx++; continue; }

      if(  (scan(0,0) - rtw < trt) && (trt < scan(0,0) + rtw) 
        && (scan(0,3) - dtw < tdt) && (tdt < scan(0,3) + dtw)){
          struct datatable table_i = process_scan_ft(rm, sym, comp, scan, tmz, tformula);
          table_i.tid = tid;
          tables.push_back( table_i );
          ftidx++;

          std::cout << std::to_string(scanid) << " ";
          std::cout << std::to_string(i) << " ";
          std::cout << std::to_string(ftidx) << " ";
          std::cout << tos(tmz, 4) << " ";
          std::cout << std::endl;
      }
    }
  }

  print_tables(tables, filenameOutput, filenameSource);
}
