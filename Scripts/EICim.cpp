#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <map>
// #include <algorithm>

#include <iterator>
#include<functional>
#include <Rcpp.h>

// [[Rcpp::export]]
std::vector< std::vector<int> > get_eic_indices(
    Rcpp::NumericVector mz,
    Rcpp::NumericVector rt,
    Rcpp::NumericVector dt,
    Rcpp::NumericVector target_mz,
    Rcpp::NumericVector target_rt,
    Rcpp::NumericVector target_dt,
    float mztolerance,
    float rttolerance,
    float dttolerance
  ){
  
  // fill indices vector of EIC indices with initial empty vectors
  std::vector< std::vector<int> > indices;
  for (int i = 0; i < target_mz.size(); i++)
  {
    std::vector<int> tmp;
    indices.push_back(tmp);
  }
  
  // For each (sorted(mz), rt, dt) data point
  //   Check each feature table within the presorted mz window.
  //     If the point is after the mz window, skip to next data point
  //     If the point is before the mz window, skip to next feature
  //     If the point is within window
  //       If the point is within all windows for that feature, add that index.
  int jstart = 0;
  for(int i = 0; i<mz.size(); i++)
  {
    for(int j = jstart; j < target_mz.size(); j++)
    {
      if( mz[i] + mztolerance < target_mz[j] ) break;
      if( target_mz[j] < mz[i] - mztolerance ) {jstart++; continue;}
      if(  (mz[i] - mztolerance < target_mz[j]) && (target_mz[j] < mz[i] + mztolerance)
        && (rt[i] - rttolerance < target_rt[j]) && (target_rt[j] < rt[i] + rttolerance)
        && (dt[i] - dttolerance < target_dt[j]) && (target_dt[j] < dt[i] + dttolerance))
        {
          indices[j].push_back(i);
        }
    }
  }
  return indices;
}

// [[Rcpp::export]]
std::vector< std::vector<int> > get_max_indices(
    std::vector< std::vector<int> > indices,
    Rcpp::NumericVector times,
    Rcpp::NumericVector it
  ){

  std::vector< std::vector< int> > newIDs;

  for (int i = 0; i < indices.size(); i++)
  {
    std::sort( begin(indices[i]), end(indices[i]), [&](int a, int b){ return times[a] < times[b]; } );

    std::map<int, int> mapt; // time, index
    for(int j = 0; j < indices[i].size(); j++){
      int t = int( times[indices[i][j]]*100 ); //bin to 0.01
      if( 0 == mapt.count(t) ){
        mapt.insert({t, indices[i][j]});
      }else{
        if ( it[ mapt[t] ] < it[ indices[i][j] ] )
          mapt[t] = indices[i][j];
      }
    }

    std::vector<int> tmp;
    for(auto const& [key, val] : mapt){
      tmp.push_back(val);
    }
    newIDs.push_back(tmp);
  }

  return newIDs;
}

// float to string with precision
std::string tos(float number, int precision){
  std::stringstream stream;
  stream << std::fixed << std::setprecision(precision) << number;
  std::string s = stream.str();
  return s;
}

// [[Rcpp::export]]
void print_eics(
    std::vector< std::vector<int> > indices,
    Rcpp::NumericVector mz,
    Rcpp::NumericVector rt,
    Rcpp::NumericVector dt,
    Rcpp::NumericVector it,
    Rcpp::NumericVector target_id,
    std::string filenameOutput,
    std::string filenameSource
  ){
    std::fstream out;
    out.open(filenameOutput, std::ios::app);
    for(int i = 0; i < target_id.size(); i++)
    {
      for(int j = 0; j < indices[i].size(); j++){
        out << target_id[i] << ",";
        out << tos(mz[indices[i][j]],5) << ",";
        out << tos(rt[indices[i][j]],2) << ",";
        out << tos(dt[indices[i][j]],2) << ",";
        out << tos(it[indices[i][j]],0) << ",";
        out << filenameSource << std::endl;
      }
    }
    out.close();
}

// [[Rcpp::export]]
void write_eics(
    Rcpp::NumericVector mz,
    Rcpp::NumericVector rt,
    Rcpp::NumericVector dt,
    Rcpp::NumericVector it,
    Rcpp::NumericVector target_mz,
    Rcpp::NumericVector target_rt,
    Rcpp::NumericVector target_id,
    Rcpp::NumericVector target_dt,
    float mztolerance,
    float rttolerance,
    float dttolerance,
    std::string filenameOutput,
    std::string filenameSource,
    bool isRTeic
  ){
  // NumericVector avoids copying over the data from R.

  std::vector< std::vector<int> > indices = get_eic_indices(mz, rt, dt
    , target_mz, target_rt, target_dt
    , mztolerance, rttolerance, dttolerance);

  // keep only indices for max intensity
  indices = get_max_indices(indices, isRTeic ? rt : dt, it);
  print_eics( indices, mz, rt, dt, it, target_id, filenameOutput, filenameSource );
}
