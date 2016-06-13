#ifndef TRAJCOMP_GEO_HPP
#define TRAJCOMP_GEO_HPP

#include <trajcomp/trajcomp.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Constants.hpp>


/* Distance between two WGS84 points
 * Code can be easily extended for other reference systems, if needed
 * 
 * */
 namespace trajcomp{

//////////////   wgs84 distance

template<class sample_type>
class wgs84_element_distance: public element_distance<sample_type,double>
{
	private:
	GeographicLib::Geodesic geod;
		
	public:
	wgs84_element_distance():geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f())
	{
	}

   virtual double operator()(sample_type u, sample_type v )
  {
	  if (u.size() != v.size()) 
		return std::numeric_limits<double>::quiet_NaN();
	  
	  double dist;
	  geod.Inverse(u[0],u[1],v[0],v[1], dist);
	  return dist;

  }
};

///////////////////// wgs84 point-line distance

template<class sample_type>
class wgs84_segment_distance: public element_segment_distance<sample_type,double>
{
	
	private:
	GeographicLib::Geodesic geod;
	double lat,lon;	// holds the nearest point on segment
		
	public:
	wgs84_segment_distance():geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f())
	{
	}
	
   virtual double operator()(sample_type u,sample_type v, sample_type p)
  {
	  // dimensions equal?
	  if (u.size() != v.size() || v.size() != p.size()) 
		return std::numeric_limits<double>::quiet_NaN();
      double a12,s12, azi1, azi2;
	  double d1,d2,d3,d4;
	  double t1=0,t4=1,t3,t2;
	  try{
			// First solve for the azimuth angle azi1
			a12 = geod.Inverse(u[0], u[1], v[0], v[1], s12, azi1, azi2);
			// Construct the Line object to be able to follow the geodesic
			GeographicLib::GeodesicLine line(geod, u[0], u[1], azi1);
				
			// Distance d1,d4
			geod.Inverse(u[0], u[1], p[0],p[1],d1);
			geod.Inverse(v[0], v[1], p[0],p[1],d4);
			
			while((t4-t1)>10E-9)		
			{
				// if middlepoint is smaller than both, 
				t2 = t1 + (t4-t1) / 3;
				t3 = t2 + (t4-t1) / 3;
				line.ArcPosition(a12 * t2, lat, lon);
				geod.Inverse(lat,lon,p[0],p[1], d2);
				line.ArcPosition(a12 * t3, lat, lon);
				geod.Inverse(lat,lon,p[0],p[1],d3);
				// d2,d3 are now the distance of two equally
				// distant inner points of t1,t4
				//cout << d1 << "<#>"<<d2 << "<#>" << d3 << "<#>" << d4 << endl;
				if ((d1 < d2 && d1 < d3 && d1 < d4) || (d2 < d3))
				{
					// refine [t1,t4] => [t1,t3]
					t4 = t3;
					d4 = d3;
				}else{
					// refine [t1,t4] => [t2,t3]
					t1 = t2;
					d1 = d2;
				}
			
			}
		}catch (const exception& e) 
		{
			throw(std::runtime_error(e.what()));
		}
		return std::min(d2,d3);
	};
	
	std::vector<double> getLastNearestPoint()
	{
		return {lat,lon};
	};
      
};




}

#endif
