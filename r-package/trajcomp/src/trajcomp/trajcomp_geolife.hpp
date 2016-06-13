#ifndef TRAJCOMP_GEOLIFE_HPP
#define TRAJCOMP_GEOLIFE_HPP

#include "trajcomp.hpp"
#include "trajcomp_files.hpp"
#include<fstream>

#include <sys/types.h>
#include <dirent.h>
#include <string.h>

#include<algorithm>


namespace trajcomp{




/*Roma dataset interface*/
template<class TrajectoryType>
class roma_reader
{
	std::string roma_location;
	public:
	roma_reader(std::string fname): roma_location(fname){};
	
	int load(std::vector<TrajectoryType> &db) //std::vector<TimeType> time
	{
		loadXYI(roma_location,db);
	}

};




	// Dataset reader with time

/*
 * Line 1...6 are useless in this dataset, and can be ignored. Points are described in following lines, one for each line.
Field 1: Latitude in decimal degrees.
Field 2: Longitude in decimal degrees.
Field 3: All set to 0 for this dataset.
Field 4: Altitude in feet (-777 if not valid).
Field 5: Date - number of days (with fractional part) that have passed since 12/30/1899.
Field 6: Date as a string.
Field 7: Time as a string.
Note that field 5 and field 6&7 represent the same date/time in this dataset. You may use either of them.

 * 
 * */
 
 
 

template<class TrajectoryType, size_t Lat=0,size_t  Lon=1,size_t Alt=2, size_t Time=3>
class geolife_reader
{
	public:
	std::string geolife_base;
	size_t MAX_DIR;
	geolife_reader(std::string base, size_t maxdir):geolife_base(base),MAX_DIR(maxdir)
	{};
	geolife_reader(std::string base):geolife_base(base),MAX_DIR(200)
	{}
	
	TrajectoryType handle_file (const char *dir, const char *fname)
	{
		if (fname[0] == '.') 
		{
			throw std::runtime_error("Invisible File");	
		}
		char buffer[1024];
	
		strncpy(buffer,dir,1024);
		strncat(buffer,fname,1024);
	
		std::ifstream infile;
		infile.open(std::string(buffer));
		
		if (!infile)
		{
			throw(std::runtime_error("File not found."));
		}
		TrajectoryType traj;	
		std::string line;
		int skip=6;	// The first 6 lines are useless (see top)
		while (std::getline(infile, line))
		{
			if (skip-->0) continue;
			typename TrajectoryType::value_type element;
			element.resize(4);
			std::replace( line.begin(), line.end(), ',', ' ');
			//std::cout << line << std::endl;
			
			std::stringstream iss(line);
			typename TrajectoryType::value_type::value_type value;
			
			
			if (! (iss >> element[ Lat ]))
			  throw(std::runtime_error("Parsing Error at Lat"));
			if (! (iss >> element[Lon]))
			  throw(std::runtime_error("Parsing Error at Lon"));
			size_t tmp;		
			if (! (iss >> tmp))
			  throw(std::runtime_error("Parsing Error at Reserved"));
			if (! (iss >> element[Alt]))
			  throw(std::runtime_error("Parsing Error at Alt"));
			if (! (iss >> element[Time]))
			  throw(std::runtime_error("Parsing Error at Time"));
				
			
			//std::cout << tools::make_string(element) << std::endl;
			
			 traj.push_back(element);
		}
   
    //traj.summary();
	return traj;
}
	
	
	int load(std::vector<TrajectoryType> &db) //std::vector<TimeType> time
	{
		std::cout << "Loading from Base " << geolife_base << std::endl;
		size_t i;
	    for (i=0; i < MAX_DIR ; i++)  
       {
            char buffer[1024];
            snprintf(buffer,1024,"%s/Data/%03d/Trajectory/",geolife_base.c_str(),i);
#ifdef VERBOSE_GEOLIFE_READER        
            std::cout <<"Going to enumerate directory" << buffer << std::endl;
#endif

            /*DIR *dp;
			struct dirent *ep;*/
			
			struct dirent **namelist;
            int n;
            n = scandir(buffer, &namelist, NULL, alphasort);
            if (n < 0)
               perror("scandir");
           else {
			   for (size_t k=0; k < n; k++)
			   {
				   //intf("%s\n", namelist[k]->d_name);
				   TrajectoryType t;
				   try{
				   t = handle_file(buffer,namelist[k]->d_name);
			      }catch(std::runtime_error &e)
			      {
					  #ifdef VERBOSE_GEOLIFE_READER        
					  cout << "Caught " << e.what() << endl;
					  cout << "Skipping " << namelist[k]->d_name << endl;
					  #endif
					  continue;
				  }
					if (t.size() != 0)
					{
						db.push_back(t);
					}else{
						std::cout << "Warning: " << buffer << "led to empty trajectory." << std::endl;
					}
				   free(namelist[k]);
			   }
                                
               
               free(namelist);
			}
			
			/* Replacing opendir with scandir, which returns sorted entries
			 * 
			 * dp = opendir (buffer);
			if (dp != NULL)
			{
				TrajectoryType t;
				while ( (ep = readdir (dp)) )
				{
					t = handle_file (buffer,ep->d_name);
					if (t.size() != 0)
					{
						db.push_back(t);
					}else{
						std::cout << "Warning: " << buffer << "led to empty trajectory." << std::endl;
					}
				}
						
				(void) closedir (dp);
			}
			else
				return i;
				//throw(std::runtime_error("geolife: Couldn't open the directory"));
				* */
		}
		return i;
	}
	
	
};

}

#endif
