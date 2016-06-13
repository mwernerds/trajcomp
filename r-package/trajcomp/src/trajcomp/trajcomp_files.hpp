#ifndef TRAJCOMP_FILES_HPP
#define TRAJCOMP_FILES_HPP

#include<string>
#include<fstream>
#include<sstream>

/*
 * Load a Trajectory File in XYI format
 * 
 * Consists of three space-separated columns with the first representing
 * X coordinate, the second representing Y coordinate and the third one
 * representing the index of the trajectory.
 * 
 * Limitations: 
 * 		- Cannot not model time
 * 		- Limited to 2D
 * 
 * Concepts and Types:
 * 		- tcoll must support clear, push_back, value_type
 * 		- tcoll::value_type must support push_back and clear()
 * */


template<typename tcoll>
size_t loadXYI(std::string fname, tcoll &s)
{
	s.clear();
	std::ifstream infile(fname);
	if (!infile)
	{
		throw(std::runtime_error("File not found."));
	}
	std::string line;
	typename tcoll::value_type t;
	int oldk = -1;
	while (std::getline(infile, line))
	{
		std::stringstream iss(line);
		double x,y; double k;
		iss >> x;
		iss >> y;
		iss >> k;
		
		if (oldk != (int) k && oldk != -1)
		{
			s.push_back(t);
			t.clear();
		}else
			t.push_back({x,y});
		oldk = (int)k;
	}
	if (t.size() != 0)
		s.push_back(t);
	return s.size();
}


/*
 * Load a Trajectory Segment File in XYXYI format
 * 
 * Consists of five space-separated columns represnting the 2D 
 * start point, the 2D endpoint and the index
 * 
 * Limitations: 
 * 		- Cannot not model time
 * 		- Limited to 2D
 * 
 * Concepts and Types:
 * 		- tcoll must support clear, push_back, value_type
 * 		- tcoll::value_type must model segments such as tLineSegment in traclus
 * 
 * */

template<typename tcoll>
size_t loadSegmentFile(std::string fname, tcoll &s)
{
	s.clear();
	std::ifstream infile(fname);
	if (!infile)
	{
		throw(std::runtime_error("File not found."));
	}
	std::string line;
	typename tcoll::value_type t;
	int oldk = -1;
	while (std::getline(infile, line))
	{
		std::stringstream iss(line);
		double x,y,x2,y2; double k;
		iss >> x;
		iss >> y;
		iss >> x2;
		iss >> y2;
		iss >> k;
		
		if (oldk != (int) k && oldk != -1)
		{
			s.push_back(t);
			t.clear();
		}else{
			typename tcoll::value_type segment({x,y},{x2,y2},k);
			t.push_back(
			   segment
			);
		}
		oldk = (int)k;
	}
	if (t.size() != 0)
		s.push_back(t);
	return s.size();
}


/*
 * template<typename point>
class tLineSegment
{
	public:
	point s;  /// the start of the segment
	point e;  /// the end of the segment
	size_t trajectory_index;  /// the index of the trajectory, where this segment has been extracted from
	int cluster; /// the cluster ID of the segment
	
	/// The default segment is created as UNCLASSIFIED
	tLineSegment(point &_s, point &_e, size_t index):s(_s),e(_e),trajectory_index(index),cluster(UNCLASSIFIED)
	{
	}
	
	friend ostream& operator<<(ostream& os, const tLineSegment<point>& dt)
	{
		os << trajcomp::tools::make_string(dt.s) << " " << trajcomp::tools::make_string(dt.e) << " " << dt.trajectory_index << " " << dt.cluster;
		return os;
	}
	
};

 * 
 * */


#endif
