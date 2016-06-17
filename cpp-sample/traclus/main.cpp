/*
 * This sample demonstrates the TRACLUS algorithm interface from C++
 * 
 * It reads a CSV from a file and outputs a CSV to stdout
 * Input is a file with columns x y i, where x/y are the coordinates
 * and i is the index (see prague.xyi for a sample)
 * 
 * Output is similar with a last column added for the class
 * 
 * This file expects trajcomp headers to be found as trajcomp/, so 
 * you can either download the complete GITHUB repository (then it will grab them 
 * from the R package source) or you can download and install the 
 * libtrajcomp C++ package (essentially copying all headers to /usr/local/include)
 * or otherwise make the used headers available to your compiler
 * 
 * You find the most updated version in the src/trajcomp directory of the
 * R package.
 * 
 * (c) 2016 M. Werner <martin@martinwerner.de>
 * */
 
#include<iostream>
#include<vector>

using namespace std;
#include<trajcomp/trajcomp.hpp>
#include<trajcomp/trajcomp_files.hpp>
#include<trajcomp/trajcomp_traclus.hpp>

using namespace trajcomp;
using namespace trajcomp::tools;


typedef std::vector<std::vector<double>> ttrajectory;
std::vector<ttrajectory> DB;


/*Now let us have a progress bar (for terminal here)*/

/*Implement my progress bar*/
class traclus_progress: public traclus_progress_visitor
{
	public:
	void operator()(size_t i, size_t m, std::string phase)
	{
		progress(i,m,phase);
		fflush(stdout);
	}
	void finish()
	{
		cout << endl;
		fflush(stdout);
	}
};




int main(int argc, char **argv)
{
	/*Parse arguments*/
	if (argc != 4){
		cout << "Usage: traclus [file] [eps] [minLines]" << endl;
		return -1;
	}
	std::string filename(argv[1]);
	double eps = atof(argv[2]);
	double minLines = atoi(argv[3]);
	cout << "Running Traclus on " << filename << " with eps=" << eps << " and " << "minLines=" << minLines <<endl;
		
	
	size_t N = loadXYI(argv[1], DB);
	cout << "Imported " << N << " trajectories" << endl;
	traclus_progress progress;
	auto result = traclus(
			DB,
			eps,
			minLines, 
			traclus_partition<ttrajectory>,  // the partitioning rule from the paper, implemented in trajcomp_traclus.hpp
			calc_dist<ttrajectory::value_type> // the "total" distance from the paper
			, progress);
   // result is a vector of tLineSegment, each with the class. So output here something about it:
   ofstream o("clustering.csv");
   o << "x1 y1 x2 y2 trajectory_index cluster" << endl;
   for (auto s:result)
   {
	   o << s << endl;
   }	
   o.close();
	
	return 0;
	
}

