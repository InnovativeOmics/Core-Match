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
  // gets the indices of all (mz,rt,dt) within the tolerances

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
  //     If the point is within all windows for that feature, add that index.
  int jstart = 0;
  for(int i = 0; i<mz.size(); i++)
  {
    for(int j = jstart; j < target_mz.size(); j++)
    {
      if( mz[i] + mztolerance < target_mz[j] ) break; // if scan's mz is too low, skip to next scan
      if( target_mz[j] < mz[i] - mztolerance ) {jstart++; continue;} // if scan's mz is too high, skip to next feature
      if( (std::abs(mz[i] - target_mz[j]) < mztolerance)
        &&(std::abs(rt[i] - target_rt[j]) < rttolerance)
        &&(std::abs(dt[i] - target_dt[j]) < dttolerance))
        {
          indices[j].push_back(i);
        }
    }
  }
  return indices;
}

// [[Rcpp::export]]
float get_min_change_in_time(
    std::vector<int> indices,
    Rcpp::NumericVector times2, // DTs or RTs
    float target_time // DTs or RTs
){
  float change_in_time = 1000000.0;
  for( int i = 0; i < indices.size(); i++ ){
    float change_in_time_i = std::abs( target_time - times2[ indices[i] ] );
    if( change_in_time_i < change_in_time ){
      change_in_time = change_in_time_i;
    }
  }
  return change_in_time;
}


// [[Rcpp::export]]
std::vector< std::vector<int> > get_max_indices(
    std::vector< std::vector<int> > indices,
    Rcpp::NumericVector times1, // RTs or DTs
    Rcpp::NumericVector times2, // DTs or RTs
    Rcpp::NumericVector target_times, // DTs or RTs
    Rcpp::NumericVector it
  ){
  // sorts indices by times
  // round times to int (but multiply by 100 first)
  //    this will be the key for a map
  // if change in time is ep from the min change in time
  //    if time key (t) isn't in map, add that scan index
  //      if t in map,
  //        if new intensity is bigger, save it
  // create new id vector (for chromatogram)
  //    push each (max) element from map
  //    add eic to vector of eic indices
  // return vector of eic indices

  std::vector< std::vector<int> > newIDs;
  float ep = 0.001; //fudge factor

  for (int i = 0; i < indices.size(); i++)
  {
    std::sort( begin(indices[i]), end(indices[i]), [&](int a, int b){ return times1[a] < times1[b]; } );
    float min_change_in_time = get_min_change_in_time(indices[i], times2, target_times[i]);

    std::map<int, int> mapt; // time, index
    for(int j = 0; j < indices[i].size(); j++){
      float change_in_time = std::abs( target_times[i] - times2[indices[i][j]] );
      if( change_in_time < min_change_in_time + ep){
        int t = int( times1[indices[i][j]]*100 ); //bin to 0.01
        if( 0 == mapt.count(t) ){
          mapt.insert({t, indices[i][j]});
        }else{
          if ( it[ mapt[t] ] < it[ indices[i][j] ] )
            mapt[t] = indices[i][j];
        }
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
  indices = get_max_indices(indices
    , isRTeic ? rt : dt
    , isRTeic ? dt : rt
    , isRTeic ? target_dt : target_rt
    , it);
  print_eics( indices, mz, rt, dt, it, target_id, filenameOutput, filenameSource );
}
