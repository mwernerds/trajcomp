#ifndef TRAJCOMP_DIR_READER_HPP
#define TRAJCOMP_DIR_READER_HPP

#include "trajcomp.hpp"
#include<fstream>

#include <sys/types.h>
 #include <sys/stat.h> 
#include <dirent.h>
#include <string.h>

#include<algorithm>
#include<iterator>
#include<string>

namespace trajcomp{

/*directory reader*/

size_t CountWords( std::string s ) 
{ 
	return ((size_t)std::count_if( s.begin(), s.end(), std::ptr_fun <int, int> ( std::isspace ) ) == s.length()) ? 0  : std::count_if( std::find_if( s.begin(), s.end(), std::not1( std::ptr_fun <int, int> ( std::isspace ) ) ), std::unique( s.begin(), s.end() ), std::ptr_fun <int, int> ( std::isspace ) ) + !std::isspace( *s.rbegin() ); 
}


template<class TrajectoryType>
class directory_reader
{
	public:
	std::string directory;
	directory_reader(std::string dir):directory(dir)
	{};
	
	
	TrajectoryType handle_file (const char *dir, const char *fname)
	{
		if (fname[0] == '.') return TrajectoryType();	// ignore ., ..  and hidden files
		char buffer[1024];
	
		strncpy(buffer,dir,1024);
		// Add trailing directory slash when needed
		if (dir[strlen(dir)-1] != '/')
		   strncat(buffer,"/",1024);
		strncat(buffer,fname,1024);
	
		std::ifstream infile;
		infile.open(std::string(buffer));
		
		if (!infile)
		{
			throw(std::runtime_error("File not found."));
		}
		TrajectoryType traj;	// 4-dimensional here
		int dimension = -1;
		
		std::string line;
		
		while (std::getline(infile, line))
		{
			typename TrajectoryType::value_type element;
			typename TrajectoryType::value_type::value_type value;
			if (dimension < 0)
			{
				dimension = CountWords(line);
				traj.dimension = dimension;
			}		
					
		    
			std::stringstream iss(line);
			element.resize(dimension);	
			for (size_t i=0; i < dimension; i++)
			    iss >> element[i];
			//std::cout << tools::make_string(element) << std::endl;
			traj.push_back(element);
		}
    //traj.summary();
	return traj;
}
	
	
	int load(std::vector<TrajectoryType> &db) //std::vector<TimeType> time
	{
		std::cout << "Loading from directory " << directory << std::endl;
		struct dirent **namelist;
            int n;
            n = scandir(directory.c_str(), &namelist, NULL, alphasort);
            if (n < 0)
               perror("scandir");
           else {
			   for (size_t k=0; k < n; k++)
			   {
				   //printf("%s\n", namelist[k]->d_name);
				   TrajectoryType t = handle_file(directory.c_str(),namelist[k]->d_name);
					if (t.size() != 0)
					{
						db.push_back(t);
					}
				   free(namelist[k]);
			   }
                                
               
               free(namelist);
			}
			return n;
	   
            /*DIR *dp;
			struct dirent *ep;
			dp = opendir (directory.c_str());
			if (dp != NULL)
			{
				TrajectoryType t;
				while ( (ep = readdir (dp)) )
				{
					t = handle_file (directory.c_str(),ep->d_name);
					if (t.size() != 0)
					{
						db.push_back(t);
					}
				}
						
				(void) closedir (dp);
			}
			else
				return i;
				//throw(std::runtime_error("geolife: Couldn't open the directory"));
				* */
	}

	
};

}

#endif
