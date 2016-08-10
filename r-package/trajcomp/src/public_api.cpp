#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
using namespace Rcpp;
// abort is following http://gallery.rcpp.org/articles/using-rcppprogress/

#include<iostream>
#include<string>
#include<sstream>
#include <string>
#include <set>
#include <exception>
#include <iostream>
// for windows (c++0x instead of c++11)
#define TRAJCOMP_DISABLE_CPPCHECK


#define TRAJCOMP_DOUGLAS_PEUCKER_TRACK_CONTRIBUTIONS

#include<trajcomp/trajcomp.hpp>
#include<trajcomp/trajcomp_traclus.hpp>
#include"edit_distance.hpp"
#include"dtw_distance.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace pt = boost::property_tree;

#include "public_api.h"



std::vector<Settings> gSettings;
Settings &getSettings(int i)
{
	if (i<0 || i >= gSettings.size())
	  throw(std::runtime_error("Invalid Handle, you have to compile the pattern first..."));
	return gSettings[i];
}



/*
 * SECTION 1: Elementary String Features
 * 
 * 
 * */



// encode global direction
std::string encode_globaldirection(std::vector<std::vector<double>> t, char startChar='A',size_t numChar=8)
{
	std::string ret("");
	for (size_t i=1; i < t.size(); i++)
	{

		double a = (atan2(t[i][1]-t[i-1][1],t[i][0]-t[i-1][0]) + M_PI) / (2*M_PI); // angle 0..1
		if (a>1) a = 1;
		if (a<0) a=0;
		int step = (int) (numChar * a); // this is 0..numChar
		if(step == numChar) step -= numChar;
		ret += (char) (startChar+step);
	}
	return ret;
}


// encode local direction
std::string encode_localdirection(std::vector<std::vector<double>> t, char startChar='A',size_t numChar=8)
{
	std::string ret("");
	for (size_t i=2; i < t.size(); i++)
	{

		double a1 = (atan2(t[i][1]-t[i-1][1],t[i][0]-t[i-1][0]) + M_PI) / (2*M_PI); // angle 0..1
		double a2 = (atan2(t[i-1][1]-t[i-2][1],t[i-1][0]-t[i-2][0]) + M_PI) / (2*M_PI); // angle 0..1
		if (a2 < a1) a2 += 1;
		double a = a2 -a1;

		if (a>1) a = 1;
		if (a<0) a=0;
		int step = (int) (numChar * a); // this is 0..numChar
		if(step == numChar) step -= numChar;
		ret += (char) (startChar+step);
	}
	return ret;
}


// encode length
std::string encode_length(std::vector<std::vector<double>> traj, double maxmag=1,  char startChar='A',size_t numChar=16)
{

  std::vector<double> t;

  double max = 0;


  for (int i=0; i < traj.size()-1; i++){
    t.push_back(sqrt(((traj[i][0]-traj[i+1][0])*(traj[i][0]-traj[i+1][0]))
                  +((traj[i][1]-traj[i+1][1])*(traj[i][1]-traj[i+1][1]))));
    // std::cout << "t[i]: " << t[i] <<std::endl;
    if(t[i]>max){
      max = t[i];
    }
  }
  // std::cout << "max: " << max << std::endl;

	std::string ret("");
	for (size_t i=0; i < t.size(); i++)
	{
		double a = t[i] / max;
		if (a<0) a=0;
		int step = (int) (numChar * a); // this is 0..numChar
		if (step >= numChar) step = numChar-1;
		ret += (char) (startChar+step);
	}
	return ret;
}

// encode global direction with a fixed start point
std::string encode_globaldirectionFSP(std::vector<std::vector<double>> t, char startChar='A',size_t numChar=8)
{
  std::string ret("");
  for (size_t i=1; i < t.size(); i++)
  {
    
    double a = (atan2(t[0][1]-t[i-1][1],t[0][0]-t[i-1][0]) + M_PI) / (2*M_PI); // angle 0..1
    if (a>1) a = 1;
    if (a<0) a=0;
    int step = (int) (numChar * a); // this is 0..numChar
    if(step == numChar) step -= numChar;
    ret += (char) (startChar+step);
  }
  return ret;
}


// encode length with a fixed start point
std::string encode_lengthFSP(std::vector<std::vector<double>> traj, double maxmag=1,  char startChar='A',size_t numChar=16)
{
  
  std::vector<double> t;
  
  double max = 0;
  
  // distance between the start point and all other points
  for (int i=0; i < traj.size()-1; i++){
    t.push_back(sqrt(((traj[0][0]-traj[i+1][0])*(traj[0][0]-traj[i+1][0]))
                       +((traj[0][1]-traj[i+1][1])*(traj[0][1]-traj[i+1][1]))));
    // std::cout << "t[i]: " << t[i] <<std::endl;
    if(t[i]>max){
      max = t[i];
    }
  }
  // std::cout << "max: " << max << std::endl;
  
  std::string ret("");
  for (size_t i=0; i < t.size(); i++)
  {
    double a = t[i] / max;
    if (a<0) a=0;
    int step = (int) (numChar * a); // this is 0..numChar
    if (step >= numChar) step = numChar-1;
    ret += (char) (startChar+step);
  }
  return ret;
}



/*
 * SECTION 1.2: Derived / Complex Features
 * 
 * 
 * */


// encode the given trajectory according to requested features and zip them 
std::string zip_features(std::vector<std::vector<double>> &q, std::vector<std::string> featureNames){
  std::vector<std::string> features;
  for(size_t i=0; i < featureNames.size(); i++){
    features.push_back(feature_dispatchbyname(q,featureNames[i]));
  }
  std::string ret;
  for (size_t i=0; i < q.size()-1; i++) // @TODO: This assumes at most |q| characters in encodings
  {
    for (size_t j=0; j < features.size(); j++)
    {
      if (features[j].size() > i){
        ret += features[j][i];
      }else{
        ret += '0';	// is before 'A' in ASCII
      }
    }
  }
  return ret;
}





/*
 * SECTION 2: Dispatching
 * 
 * 
 * */

/*
 * SECTION 2.1: String Feature Dispatching
 * 
 * 
 * */

std::string feature_dispatchbyname(std::vector<std::vector<double>> &q, std::string name){
  if(name == "global_direction"){
    return encode_globaldirection(q);
  }
  if(name == "global_direction_FSP"){
    return encode_globaldirectionFSP(q);
  }
  if(name == "local_direction"){
    return encode_localdirection(q);
  }
  if(name == "length"){
    return encode_length(q);
  }
  if(name == "length_FSP"){
    return encode_lengthFSP(q);
  }
//   if (name == "zipall")
//   {
// 	  std::vector<std::string> features = {encode_globaldirection(q),encode_localdirection(q),encode_length(q)};
// 	  std::string ret;
// 	  for (size_t i=0; i < q.size(); i++) // @TODO: This assumes at most |q| characters in encodings
// 	  {
// 		  for (size_t j=0; j < features.size(); j++)
// 		  {
// 			  if (features[j].size() > i){
// 				   ret += features[j][i];
// 			   }else{
// 				   ret += '0';	// is before 'A' in ASCII
// 			   }
// 		  }
// 	  }
// 	  return ret;
//   }

  throw(std::runtime_error("String Encoding "+name+" has not been defined "));
  return "";
}



/*
 * SECTION 2.2: Distance Dispatching
 * 
 * 
 * */




double distance_dispatchbyname(std::vector<std::vector<double>> &t, std::vector<std::vector<double>> &q, std::string name)
{
	if(name == "discrete_frechet"){
		return trajcomp::discrete_frechet(t,q);
	}
	if(name == "closest_pair"){
		return trajcomp::closest_pair(t,q);
	}
	if (name == "sum_of_pairs")
	{
	  double ret =-1;
	  try{ret =  trajcomp::sum_of_pairs(t,q);}
	  catch(...){}
	  return ret;
	}
	if(name == "dtw"){
	  return dtw_distance::dtw(t,q);
	}
	if (name == "discrete_hausdorff")
	{
		return trajcomp::discrete_hausdorff(t,q);
	}

	return -1;

}





/*
 * APPENDIX A: Debugging, Miscellanous, ...
 * 
 * 
 * */

// auxilarity function to call the R function for AAA_distance from C++
double AAA_distance_function(std::vector<std::string> t, std::vector<std::string> q,Function f){

  double dist = as<double>(f(t,q));
  return dist;
}

// [[Rcpp::export]]
double string_distance_dispatchbyname(std::string t, std::string q, std::string name, int k)
{
  if(name == "edit_distance"){
    return edit_distance::ed(edit_distance::explode_const(t,k),edit_distance::explode_const(q,k));
  }
  if(name == "AAA_distance"){

    Function f("AAA_distance");
    auto ret = AAA_distance_function(edit_distance::explode_const(t,k), edit_distance::explode_const(q,k),f);

    return ret;
  }
  throw(std::runtime_error("String Distance "+name+" has not been defined "));

}




// read setting file
Settings readSettings(std::string filename){
  
  // priority tree
  pt::ptree tree;
  pt::read_xml(filename, tree);
  
  Settings s;
  
  for(auto &child: tree.get_child("group")){
    
    if(child.first == "distance"){
      s.distances.push_back(child.second.get<std::string>("<xmlattr>.type"));
      s.lengths.push_back(child.second.get<int>("<xmlattr>.length"));
      
      std::vector<std::string> features;
      
      if(child.second.get<int>("<xmlattr>.length") != 0){
        for(auto &child2: child.second){
          if(child2.first == "feature"){
            features.push_back(child2.second.get<std::string>("<xmlattr>.type"));
          }
        }
        s.features.push_back(features);
      }
      else{
        features.push_back("null");
        s.features.push_back(features);
      }  
    }
  }
  
//   for(int i=0; i < s.distances.size(); i++){
//     std::cout << "distances: " << s.distances[i] << std::endl;
//   }
//   for(int i=0; i < s.lengths.size(); i++){
//     std::cout << "lengths: " << s.lengths[i] << std::endl;
//   }
//   for(int i=0; i < s.features.size(); i++){
//     for(int j = 0; j < s.features[i].size(); j++){
//       std::cout << "features" << i << ": " << s.features[i][j] << std::endl;
//     }
//   }
//   
  return s;
}

//' compiles a valid XML distance description into a handle
//' 
//' @param setting an XML string describing the distance to be used
//' @param setting an XML string describing the distance to be used
//' @return An integer handle to pass the (compiled) XML settings to other functions
//' @details
//' Distance functions are described via a specific XML format. Consider reading 
//' the documentation at \url{trajectorycomputing.com/libtrajcomp-xml} 
//' @export
// [[Rcpp::export]]
size_t compileSettings(std::string setting)
{

	Settings s;
	std::stringstream ss;
	ss << setting;
		
	
    pt::read_xml(ss, s.tree);
	auto &tree = s.tree;  //@todo: later refactor tree => s.tree
	// Parsing helper fields
  for(auto &child: tree.get_child("group")){
    
    if(child.first == "distance"){
      s.distances.push_back(child.second.get<std::string>("<xmlattr>.type"));
      s.lengths.push_back(child.second.get<int>("<xmlattr>.length"));
      
      std::vector<std::string> features;
      
      if(child.second.get<int>("<xmlattr>.length") != 0){
        for(auto &child2: child.second){
          if(child2.first == "feature"){
            features.push_back(child2.second.get<std::string>("<xmlattr>.type"));
          }
        }
        s.features.push_back(features);
      }
      else{
        features.push_back("null");
        s.features.push_back(features);
      }  
    }
  } // for child
  gSettings.push_back(s);
  return (gSettings.size() -1);
}





// [[Rcpp::export]]
SEXP trajectory_distances(NumericMatrix S, NumericMatrix T, int hSetting)
{
  
  Settings s = getSettings(hSetting);

  std::vector<std::vector<double>> trajectory,query;
  for (size_t i=0; i < S.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({S(i,0),S(i,1)});
  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    query.push_back({T(i,0),T(i,1)});
  
  
  std::vector<double> distances;
  std::vector<std::string> main;
  
  // for each distance...
  for(int i = 0; i < s.distances.size(); i++){
  
    // if there are no features to encode, dispatch by distance name
    if(s.lengths[i] == 0){
	   distances.push_back(
        distance_dispatchbyname(
          trajectory, query, s.distances[i]));
          
     main.push_back(s.distances[i]);
    
    }
    else{
	
  
      std::vector<std::string> featureNames;
      // encode trajectory according to requested feature(s)
      distances.push_back(string_distance_dispatchbyname(
                  zip_features(trajectory, s.features[i]),
                  zip_features(query, s.features[i]),
                  s.distances[i], s.lengths[i]));
      

      std::string features = boost::algorithm::join(s.features[i], ", ");
      main.push_back(s.distances[i]+": "+features);
    }
  }
    
  return List::create(Named("distances")=distances,
                      Named("main")=main);
}


// /*
//  * Martin
//  * 
//  * 
//  * */
// // [[Rcpp::export]]
// NumericVector trajectory_distance_vector(NumericMatrix S, NumericMatrix T, int hSetting)
// {
//   
//   Settings s = getSettings(hSetting);
//   
//   if (s.distances.size() != 1)
//     throw(std::runtime_error("trajectory_distance_vector can only use one distance at a time"));
//   
//   std::vector<std::vector<double>> trajectory,query;
//   
//   for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
//     query.push_back({T(i,0),T(i,1)});
//   
//   std::vector<double> distances(S.nrow());
//   std::vector<size_t> distributor;
//   
//   for (size_t i=0; i <= S.nrow(); i++)		
//   {
//     if (i == S.nrow()|| std::isnan(S(i,0))) // the last or any nan
//     {
//       //work (else collect)
//       // for each distance... it is only one ;-)
//       double d;
//       for(int i = 0; i < s.distances.size(); i++){
//         
//         // if there are no features to encode, dispatch by distance name
//         if(s.lengths[i] == 0){
//           d = 
//             distance_dispatchbyname(
//               trajectory, query, s.distances[i]);
//           
//           
//           
//         }
//         else{
//           std::vector<std::string> featureNames;
//           // encode trajectory according to requested feature(s)
//           d = string_distance_dispatchbyname(
//             zip_features(trajectory, s.features[i]),
//             zip_features(query, s.features[i]),
//             s.distances[i], s.lengths[i]);
//         }
//       }  
//       if (i < distances.size()) distances[i] = S(i,0); // make it NaN
//       for (auto p:distributor) distances[p] = d;
//       
//       trajectory.clear();
//       distributor.clear();
//       
//     }else{ // collect!
//       trajectory.push_back({S(i,0),S(i,1)});
//       distributor.push_back(i);
//     } 
//     
//   } // foreach trajectory
//   
//   return wrap (distances);
// }

// [[Rcpp::export]]
String encode(NumericMatrix T, std::string feature){


  std::vector<std::vector<double>> trajectory;
  
  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});
  
  // for each feature ...
  std::string ret = feature_dispatchbyname(trajectory, feature);
  

  return ret;
}

/*
 * 
 * returns NumericVector: distance1, distance2, ...
 * 
 * */
// [[Rcpp::export]]
NumericVector distance_vector_ddply(NumericMatrix T, NumericMatrix Q, int hSetting)
{
    
    Settings s = getSettings(hSetting);

    std::vector<std::vector<double>> trajectory,query;
    
    for (size_t i=0; i < Q.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      query.push_back({Q(i,0),Q(i,1)});
    
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    std::vector<double> distances;
    // std::cout << s.distances.size() << std::endl;
    
    
    // for each distance...
    double d;
    for(int i = 0; i < s.distances.size(); i++){
      // if there are no features to encode, dispatch by distance name
      if(s.lengths[i] == 0){
        d = distance_dispatchbyname(
                trajectory, query, s.distances[i]);
        distances.push_back(d);

      }
      else{
        std::vector<std::string> featureNames;
        // encode trajectory according to requested feature(s)

        d = string_distance_dispatchbyname(
              zip_features(trajectory, s.features[i]),
              zip_features(query, s.features[i]),
              s.distances[i], s.lengths[i]);
//         std::cout << d << std::endl;
        distances.push_back(d);
      }
      
    }
//     NumericMatrix res(distances[0].size(), distances.size());
//     for(int i=0; i<distances.size(); i++){
//       NumericVector temp = wrap(distances[i]);
//       res(_,i) = temp;
//     }
//     
    return wrap(distances);
}


/*
 * auxilarity function to get the names of the distances as vector of strings
 * 
 * 
 * */

// [[Rcpp::export]]
CharacterVector getDistanceNames(int hSetting){
  
  Settings s = getSettings(hSetting);
  
  CharacterVector names(s.distances.size());
  
  for(int i=0; i<s.distances.size(); i++){
    std::string distanceName = s.distances[i];
    if(s.lengths[i]==0){
      names[i] = distanceName;
    }
    else{
      names[i] = distanceName + ": " + boost::algorithm::join(s.features[i], ", ");
    }
  }
  
  return names;
  
}


/*
 * Douglas Peucker
 * 
 * 
 * */

// [[Rcpp::export]]
NumericMatrix DouglasPeucker(NumericMatrix S, double epsilon) {
  

  std::vector<std::vector<double>> trajectory;

  
  for (size_t i=0; i < S.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({S(i,0),S(i,1)});
  
  
  std::vector<std::vector<double>> res = trajcomp::douglas_peucker(trajectory, epsilon);
  
  
  NumericMatrix resultMatrix = NumericMatrix(res.size(), res[0].size()) ;
  for(int i = 0; i < res.size(); i++){
    NumericVector temp = wrap(res[i]); 
    resultMatrix(i,_) = temp;
  }

//   std::cout << res[res.size()-1][0] << ", " << res[res.size()-1][1] << std::endl;
  
    
  return resultMatrix;
  
}


// [[Rcpp::export]]
NumericMatrix DouglasPeuckerWeights(NumericMatrix S, double epsilon) {
  

  std::vector<std::vector<double>> trajectory;

  
  for (size_t i=0; i < S.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({S(i,0),S(i,1)});
  
  
  std::vector<double> res = trajcomp::douglas_peucker_weights(trajectory, epsilon);
  
  
  NumericMatrix resultMatrix = NumericMatrix(res.size(), 1) ;
  for(int i = 0; i < res.size(); i++){
    NumericVector temp = wrap(res[i]); 
    resultMatrix(i,_) = temp;
  }

//   std::cout << res[res.size()-1][0] << ", " << res[res.size()-1][1] << std::endl;
  
    
  return resultMatrix;
  
}



/*
* Linear Interpolation
* 
* 
* */
std::pair<int, int> getNeighborIndices(NumericVector tAxis, double t, int nOld){
  int index = 0;
  while(t > tAxis(index)){
    index++;
  }
  return std::make_pair(index-1, index);
}


// [[Rcpp::export]]
NumericMatrix doPLI(NumericMatrix data, int n) {
  

  int nOld = data.nrow();
  

  NumericMatrix rescaledData(n,3);
  
  
  double tau;
  double t;
  
  double tFirst = data(0,2);
  double tLast = data(nOld-1,2);
  double intervall =  (tLast-tFirst)/(n-1);
  
  std::pair<int, int> neighbors;
  int tBefore;
  int tAfter;
  
  
  // first and last element does not change
  rescaledData(0,_) = data(0,_);
  rescaledData(n-1,_) = data(nOld-1,_);
  
  
  for(int i = 1; i < n-1; i++){
    t = tFirst + intervall*i;
    rescaledData(i,2) = t;
    
    std::tie(tBefore,tAfter) = getNeighborIndices(data(_,2), t, n);
    // std::cout << "tBefore: " << tBefore << ", tAfter: " << tAfter << std::endl;
    
    tau = (t-data(tBefore, 2))/(data(tAfter, 2)-data(tBefore, 2));
    // std::cout << "tau: " << tau << std::endl;
    
    for(int j = 0; j < 2; j++){
      rescaledData(i,j) = (1-tau)*data(tBefore, j) + tau*data(tAfter, j);
      // std::cout << "i: " << i << ", j: "<< j << ", "<< rescaledData(i,j) << std::endl;
    }
    
  }
  return rescaledData;
}


//~ 
//~ static void chkIntFn(void *dummy) {
  //~ R_CheckUserInterrupt();
//~ }
//~ 
//~ // this will call the above in a top-level context so it won't longjmp-out of your context
//~ bool checkInterrupt() {
  //~ return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
//~ }
//~ 
//~ struct  R_interrupt  {};




/*Instead of actually providing a progress bar, we abuse this
 * to handle STRG-C from R*/
class traclus_progress: public trajcomp::traclus_progress_visitor
{
	
	
	public:
	//~ traclus_progress():p(10,false){}
	
	void operator()(size_t i, size_t m, std::string phase)
	{
		//~ This is not yet working.
		//~ if ( Progress::check_abort() )
		//~ {
		   //~ std::cout << "Aborting" << std::endl;	
		   //~ 
		   //~ throw(std::runtime_error("Abort"));
		//~ }else{
			//~ std::cout << "Not aborted" << std::endl;
		//~ }
		
	}
	void finish()
	{
		//~ cout << endl;
		//~ fflush(stdout);
	}
};





//' Clusters trajectories according to the TRACLUS framework
//' 
//' @param TrajectoryDB a trajectory database giving 2D trajectories split by NaN or NA
//' @param setting an XML string describing the distance to be used
//' @return An integer handle to pass the (compiled) XML settings to other functions
//' @details
//' Distance functions are described via a specific XML format. Consider reading 
//' the documentation at \url{trajectorycomputing.com/libtrajcomp-xml} 
//' @export
// [[Rcpp::export]]
DataFrame cpp_traclus(NumericMatrix TrajectoryDB, double eps, int minLines)
{
	std::vector<std::vector<std::vector<double>>> DB;
	/* First, copy trajectory database into c++ vectors
	 * This could be optimized by providing a mapping class allowing for 
	 * STL-style access to a trajectory DB given as a Matrix
	 * */
	std::vector<std::vector<double>> trajectory;
    for (size_t i=0; i < TrajectoryDB.nrow(); i++)		//@todo: remove copy by an adapter class
    {
		if (NumericMatrix::is_na (TrajectoryDB(i,0))){ // this is true for NaN or NA
			DB.push_back(trajectory);
			trajectory.clear();
		}else{
			trajectory.push_back({TrajectoryDB(i,0),TrajectoryDB(i,1)});
		}
    }
	Rcout << "Have " << DB.size() << " trajectories" << std::endl;
	
	
	traclus_progress progress;
	std::vector<trajcomp::tLineSegment<std::vector<double>>> result;
	//~ try{
		
	result = trajcomp::traclus(
			DB,
			eps,
			minLines, 
			trajcomp::traclus_partition<std::vector<std::vector<double>>>,  // the partitioning rule from the paper, implemented in trajcomp_traclus.hpp
			trajcomp::calc_dist<std::vector<double>> // the "total" distance from the paper
			, progress);
   // result is a vector of tLineSegment, each with the class. So output here something about it:
   //~ }catch(...)
   //~ {
	   //~ throw(std::runtime_error("Calculation Interrupted"));
   //~ }
   
   NumericVector x1(result.size()),y1(result.size()),x2(result.size()),y2(result.size()),
				trajectory_index(result.size()), cluster(result.size());
   
   
   for (size_t i=0; i < result.size(); i++)
   {
	   auto &r = result[i];
	   x1[i] = r.s[0]; y1[i] = r.s[1];
	   x2[i] = r.e[0]; y2[i] = r.e[1];
	   trajectory_index[i] = r.trajectory_index;
	   cluster[i] = r.cluster;
   }

    return Rcpp::DataFrame::create( 
			Named("x1")= x1, 
			Named("y1") = y1,
			Named("x2") = x2,
			Named("y2") = y2,
			Named("trajectory_index") = trajectory_index,
			Named("cluster") = cluster
			);
}
