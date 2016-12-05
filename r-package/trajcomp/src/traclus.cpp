#include <Rcpp.h>
using namespace Rcpp;
 
 #include<trajcomp/trajcomp.hpp>
 #include<trajcomp/trajcomp_traclus.hpp>
#include "interrupt.h"
 
/*TRACLUS R INTERFACE */


/*Instead of actually providing a progress bar, we abuse this
 * to handle STRG-C from R*/
class traclus_progress: public trajcomp::traclus_progress_visitor
{
	
	
	public:
	//~ traclus_progress():p(10,false){}
	
	void operator()(size_t i, size_t m, std::string phase)
	{
		if (checkInterrupt())
		{
			Rcout << "Should abort" << std::endl;
			stop("Using stop on that");
		}
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

