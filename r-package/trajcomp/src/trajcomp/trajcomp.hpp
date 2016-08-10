#ifndef TRAJCOMP_HPP_INC
#define TRAJCOMP_HPP_INC

#ifndef TRAJCOMP_DISABLE_CPPCHECK
#if __cplusplus < 201103L
  #warning --------------------------------------
  #warning This library needs at least a C++11 compliant compiler
  #warning If you believe that it should work, define TRAJCOMP_DISABLE_CPPCHECK before
  #warning including trajcomp header libraries.
  #warning  --------------------------------------
  #error "Stopping compilation (define TRAJCOMP_DISABLE_CPPCHECK to disable check)"
#endif
#endif



#include<vector>
#include<iostream>
#include<stdexcept>
#include<string>
#include<sstream>
#include<fstream>
#include<limits>
#include<sys/time.h>
#include<stdint.h>

#ifdef TRAJCOMP_SERIALIZATION
	// forward declaration for friendship
	namespace boost {
	namespace serialization {
		class access;
	}
	}
#endif


#include<math.h> // for sqrt
#define COL0 "\033[G"
#define CTLK "\033[K"
//using namespace std;


//#define abs(x) (x<0?-x:x)

namespace trajcomp
{
	#define MIN(x,y) ((x<y)?x:y)
	#define MAX(x,y) ((x>y)?x:y)
	namespace tools{
		
		
		template<class datatype, char delim=' '>
std::string make_string2(datatype data)
{
		std::stringstream ss;
		typename datatype::iterator it;
		ss << std::endl;
		for(it=data.begin(); it != data.end(); it++)
		{
		  typename datatype::value_type::iterator it2;
		  for (it2 = (*it).begin();
			   it2 != (*it).end(); it2++){
				  
				ss << "\t" << *it2 << delim;
			}ss  << std::endl;
		}
			return ss.str();
	}
		
		template<class datatype, char delim=' '>
		std::string make_string(datatype data)
		{
			std::stringstream ss;
			typename datatype::iterator it;
			for(it=data.begin(); it != data.end(); it++)
				ss <<  *it << delim;
			return ss.str();
		}

		void progress(long pos, long size,std::string what="")
		{
			std::cout << COL0 << "Progress("<<what<<"): " << pos << "/" << (double)size << "(" << 100*(double)pos/(double)size << "%)" << CTLK;
		}
		uint64_t ticks(void)
		{	
			struct timeval tv;
			gettimeofday(&tv, 0);
			return uint64_t( tv.tv_sec ) * 1000 + tv.tv_usec / 1000;
		}
		class tictoc
		{
			public:
			uint64_t t;
			void tic() {t = ticks();};
			uint64_t toc(){return ticks() - t; };
		};
			
template<class datatype> 
void matrix_resize(std::vector<std::vector< datatype> > &m, int r, int c)
{
	m.resize(r);
	for(size_t i=0; i <r; i++)
	  m[i].resize(c);	
}
template<class datatype> 
void matrix_resize(std::vector<std::vector<std::vector< datatype > > > &m, int r, int c, int d)
{
	m.resize(r);
	for(size_t i=0; i <r; i++)
	{
	  m[i].resize(c);	
	  for (size_t j=0; j < c; j++)
	    m[i][j].resize(d);
    }
    
}
		


	} // tools
	
template <class ElementType, class DistanceType=double>
class element_distance 
{
public:
   double operator() (void ){return 0;};
   virtual double operator()(ElementType u, ElementType v )=0;
};



template <class ElementType, class DistanceType=double>
class element_segment_distance 
{
public:
   virtual DistanceType operator()(ElementType s1, ElementType s2, ElementType p)=0;
};

/*Default distances for containers of equal length*/

template<class sample_type>
class default_element_distance: public element_distance<sample_type,double>
{
	public:
	double operator() (void ){return 0;};
   virtual double operator()(sample_type u, sample_type v )
  {
	  if (u.size() != v.size()) 
		return std::numeric_limits<double>::quiet_NaN();
      double sum=0;
      for (size_t i=0; i < u.size(); i++)
        sum += (u[i]-v[i])*(u[i]-v[i]);

	return sqrt(sum);
  }
};

template<class sample_type>
class default_element_distance_squared: public element_distance<sample_type,double>
{
	public:
   virtual double operator()(sample_type u, sample_type v )
  {
	  if (u.size() != v.size()) 
		return std::numeric_limits<double>::quiet_NaN();
      double sum=0;
      for (size_t i=0; i < u.size(); i++)
        sum += (u[i]-v[i])*(u[i]-v[i]);

	return sum;
  }
};


template<class sample_type>
class default_segment_distance: public element_segment_distance<sample_type,double>
{
	public:
   virtual double operator()(sample_type u,sample_type v, sample_type p)
  {
	  // dimensions equal?
	  if (u.size() != v.size() || v.size() != p.size()) 
		return std::numeric_limits<double>::quiet_NaN();
      
	  // the length of the segment
	  default_element_distance_squared<sample_type> d2;
	  
	  double l = d2(u,v);
#ifdef DEBUG_default_segment_distance
	  std::cout << "l=" << l << endl;
#endif

	  if (fabs(l) < 1E-12)
			return sqrt(d2(p,v));
		
		double t = 0;
		for (size_t i=0; i <  u.size(); i++)
		   //t += ((u[i]-v[i])*(p[i]-v[i])); @TODO:Remove?
		   t += ((v[i]-u[i])*(p[i]-u[i]));
		t/=l;
#ifdef DEBUG_default_segment_distance
	std::cout << "t= " << t << endl;
#endif		
			
		if (t<0) return sqrt(d2(p,u));
		if (t>1) return sqrt(d2(p,v));
		
		sample_type proj;
		for (size_t i=0; i <  u.size(); i++)
		   //proj.push_back(v[i]+t*(u[i]-v[i]));
		   proj.push_back(u[i]+t*(v[i]-u[i]));
		double ret = sqrt(d2(p,proj));
#ifdef DEBUG_default_segment_distance
		std::cout << "(" << 
				tools::make_string(u) << ";" <<
				tools::make_string(v) << ";" <<
				tools::make_string(p) << ") ==>";
		std::cout << tools::make_string(proj) << "=" << ret <<  std::endl;
#endif
		return ret;
  }
};




	
	// A trajectory is extended from a vector of elements of a dimension
	template <typename ValueType>
	class trajectory : public std::vector<std::vector<ValueType>>
	{
		
		public:
		unsigned int dimension;
		typedef std::vector<std::vector<ValueType>> Base;
		typedef std::vector<ValueType> value_type;
						
		trajectory():dimension(0)
		{};
		trajectory(size_t dim):dimension(dim)
		{};
		void summary()
		{
			std::cout << "Dimension = " << dimension << std::endl;
			std::cout << "Elements  = " << this->size() << std::endl;
		}  
		
		void save(std::string filename)
		{
			std::ofstream of(filename);
			if (of)
			{
				for (typename Base::iterator it=this->begin(); it != this->end(); it++)
				{
					std::string s = tools::make_string(*it);
					of << s << std::endl;   
					
				}
			}else
			{
				throw(std::runtime_error("Error: Can't save to file."));
			}
			
		}
		
		
		void load(std::string filename)
		{
			this->clear();
			std::ifstream infile(filename);
			if (!infile)
			{
				throw(std::runtime_error("File not found."));
			}
			std::string line;
			while (std::getline(infile, line))
			{
				value_type l;
				std::stringstream iss(line);
				ValueType v;
				while (iss >> v) 
					l.push_back(v);
				this->dimension=l.size();
				this->push_back(l);
			}
		}
		
		
		void dump()
		{
			for (typename Base::iterator it=this->begin(); it != this->end(); it++)
		   {
			  std::string s = tools::make_string(*it);
			  std::cout << s << std::endl;   
		   }
		}
		
		void push_back (const typename Base::value_type& val)
		{
			/*if(val.size() != dimension)
			 throw(std::runtime_error("Error: Can't add element of wrong dimension to trajectory."));*/
			Base::push_back(val);
			 
		}
		void push_back (typename Base::value_type && val)
		{
			/*if(val.size() != dimension)
			 throw(std::runtime_error("Error: Can't add element of wrong dimension to trajectory."));*/
			Base::push_back(val);
		}
		
		/*SERIALIZATION*/
		 
#ifdef TRAJCOMP_SERIALIZATION
		 // Allow serialization to access non-public data members.
		friend class boost::serialization::access;
		
		template<typename Archive>
		void serialize(Archive& ar, const unsigned version) {
				ar & d1_ & d2_;  // Simply serialize the data members of Obj
		}
#endif		
		
	};

	
//3.2.1. Uniform Select.
// 	Extracts each k-th point

template<class TrajectoryType>
TrajectoryType uniform_select(TrajectoryType &t,size_t Step=15,bool withEndPoint=true)
{
	TrajectoryType ret(t.dimension);
	typename TrajectoryType::iterator it;
	size_t k=0;
    for (it = t.begin(); it != t.end(); it++)
    {
		if (k % Step == 0)
			ret.push_back(*it);
		k++;
	}
	
	if(withEndPoint)
	  if (((t.size()-1) % Step) != 0)
	  {
	    ret.push_back(t.back());
	}
	return ret;
}

/*
template<class TrajectoryType, typename functional>
TrajectoryType middle_distance_interpolation(TrajectoryType &t,double d_max,
											functional f)
{
	TrajectoryType ret(t.dimension);
	typename TrajectoryType::iterator it;
	
	first = t.begin(); 
	auto second = first;
	second ++;
    for (; second != t.end();)
    {
		double seglen = f(*first, *second);
		
		first ++;
		second ++;
	}
	
	return ret;
}
*/


//3.2.2. Douglas Peucker Algorithm.
// Implemented using a complex class to encapsulate methods for
// recursion and a default function template

template <class TrajectoryType, class DistanceType=double>
class douglas_peucker_impl
{
	public:
	std::vector<unsigned char> booleanMap; // @REMARK: could also use bool specialization
	std::vector<DistanceType> numericMap;  // @REMARK: could be split away if picky about the memory
	element_segment_distance<typename TrajectoryType::value_type, DistanceType> *d;
public:
   douglas_peucker_impl():d(NULL){};
DistanceType largest_contribution(TrajectoryType &u, size_t istart, size_t iend,size_t &the_index)
   {
	// find the point which contributes most
	
	DistanceType the_contribution = -1;
	typename TrajectoryType::value_type s = u[istart];	
	typename TrajectoryType::value_type e = u[iend];
	
	for (size_t i=istart+1; i < iend; i++) 
	{
		
		DistanceType contribution = (*d)(s,e,u[i]);
#ifdef DEBUG_largest_contribution
		std::cout << i << "~" << contribution << std::endl;
#endif		
		if (contribution > the_contribution)
		{
			the_contribution = contribution;
			the_index = i;
		}
	}

	return the_contribution;
   }
   
   void simplify(TrajectoryType &u, size_t istart, size_t iend, DistanceType max_error)
   {
	   // simplify segment and recurse
#ifdef DEBUG_douglas_peucker_simplify	   
	   std::cout << istart << "," << iend << std::endl;
#endif
      
	   
	   if ((istart -iend) < 2) return; // there is nothing I can add.
       

	   size_t the_index=0;
	   DistanceType the_contribution = largest_contribution(u,istart,iend,the_index);
	   if (the_contribution > max_error)
	   {
#ifdef DEBUG_douglas_peucker_simplify
		   std::cout << "DP: Insert " << the_index << std::endl;
#endif
#ifdef TRAJCOMP_DOUGLAS_PEUCKER_TRACK_CONTRIBUTIONS
           numericMap[the_index] = the_contribution;
#endif

		   booleanMap[the_index] = 1;
		   simplify(u,istart,the_index,max_error);
		   simplify(u,the_index,iend,max_error);
	   }
   }

   TrajectoryType operator()(TrajectoryType u, DistanceType max_error,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {
	auto which = indizes(u,max_error,d);

	TrajectoryType ret;
	for (auto w:which)
	  ret.push_back(u[w]);

    return ret;    
  }
  


  std::vector<size_t> indizes(TrajectoryType u, DistanceType max_error,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {
	this->d = &d;
    
	booleanMap.resize(u.size());
	for (size_t i=0; i < u.size(); i++)
		booleanMap[i] = 0;
	booleanMap[0] = booleanMap[u.size()-1] = 1;
    numericMap.resize(u.size());
	simplify(u,0,u.size()-1,max_error);
	//@TODO: call indizes for above.
	// now build the simplified trajectory
	std::vector<size_t> ret;
	for (size_t i=0; i < booleanMap.size(); i++)
	  if (booleanMap[i])
	     ret.push_back(i);
	     
	 return ret;

  }
  
  
  std::vector<DistanceType> weights(TrajectoryType u, DistanceType max_error,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {
	this->d = &d;
    
	booleanMap.resize(u.size());
	numericMap.resize(u.size());
	for (size_t i=0; i < u.size(); i++)
	{
		booleanMap[i] = 0;
		numericMap[i] = 0;
	}
	booleanMap[0] = booleanMap[u.size()-1] = 1;
  
	simplify(u,0,u.size()-1,max_error);
	return numericMap;
  }
  
};



template<class TrajectoryType>
TrajectoryType douglas_peucker(TrajectoryType &t, double epsilon,
element_segment_distance<typename TrajectoryType::value_type,double> &d)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}	


template<class TrajectoryType>
TrajectoryType douglas_peucker(TrajectoryType &t, double epsilon)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}	


template<class TrajectoryType>
std::vector<size_t> douglas_peucker_indizes(TrajectoryType &t, double epsilon)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	std::vector<size_t> ret = dp.indizes(t,epsilon,d);
	return ret;
}	


template<class TrajectoryType>
std::vector<double> douglas_peucker_weights(TrajectoryType &t, double epsilon)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	
	return dp.weights(t,epsilon,d);
}	


//3.2.3. Bellmans Algorithm.
//3.2.4. Bottom up.
//3.2.5. Reservoir Sampling.
//3.2.6. Sliding Window Algorithm.
	
//////////////////// ELEMENTARY DISTANCE ALGORITHMS

//// POINT TO TRAJECTORY DISTANCE

template<class TrajectoryType>
typename TrajectoryType::value_type::value_type point_to_trajectory_point(typename TrajectoryType::value_type &p, TrajectoryType &t)
{
	default_element_distance<typename TrajectoryType::value_type> d;
	typename TrajectoryType::value_type::value_type m = 0;
	bool first = true;
	
	for (typename TrajectoryType::iterator it=t.begin(); it != t.end(); it++)
	{
		typename TrajectoryType::value_type::value_type v=d(p,*it);
		if (first || v<m) 
		{
			   m=v;
			   first=false;
		}
			   
	}
	return m;
}	

template<class TrajectoryType>
typename TrajectoryType::value_type::value_type point_to_trajectory(typename TrajectoryType::value_type &p, TrajectoryType &t)
{
	default_segment_distance<typename TrajectoryType::value_type> d;
	typename TrajectoryType::value_type::value_type m = 0;
	
	bool first = true;
	for (size_t i=0; i < t.size()-1; i++)
	{
		typename TrajectoryType::value_type::value_type v=d(t[i],t[i+1],p);
		if (first || v<m) 
		{
			   m=v;
			   first=false;
		}
		
	}

	return m;
}	




template<class TrajectoryType>
typename TrajectoryType::value_type::value_type closest_pair(TrajectoryType &u, TrajectoryType &v)
{
	default_element_distance<typename TrajectoryType::value_type> d;
	typename TrajectoryType::value_type::value_type m = 0;
	bool first = true;
	
	for (typename TrajectoryType::iterator it1=u.begin(); it1 != u.end(); it1++)
	for (typename TrajectoryType::iterator it2=v.begin(); it2 != v.end(); it2++)
	{
		typename TrajectoryType::value_type::value_type v=d(*it1,*it2);
		if (first || v<m) 
		{
			   m=v;
			   first=false;
		}
			   
	}
	return m;
}


template<typename TrajectoryType>
typename TrajectoryType::value_type::value_type discrete_D(typename TrajectoryType::value_type p, TrajectoryType &v)
{
	typename TrajectoryType::value_type::value_type ret = std::numeric_limits<typename TrajectoryType::value_type::value_type>::infinity();
	default_element_distance<typename TrajectoryType::value_type> d;
	for (auto q:v)
		ret = std::min(ret, d(p,q));
	return ret;
}


template<typename TrajectoryType>
typename TrajectoryType::value_type::value_type discrete_hausdorff(TrajectoryType &u, TrajectoryType &v)
{
	typename TrajectoryType::value_type::value_type ret = -std::numeric_limits<typename TrajectoryType::value_type::value_type>::infinity();
	for (auto p:u)
		ret = std::max(ret, discrete_D(p,v));
	for (auto p:v)
		ret = std::max(ret, discrete_D(p,u));
	return ret;
}






template<class TrajectoryType>
class implicit_distance_matrix
{
	private:
	size_t cols;
	std::vector<typename TrajectoryType::value_type::value_type> D;
	TrajectoryType &u;
	TrajectoryType &v;
	default_element_distance<typename TrajectoryType::value_type> d;
	
	public:
		implicit_distance_matrix(TrajectoryType &theu,TrajectoryType &thev): u(theu), v(thev)
		{
			cols = v.size();
			D.resize(u.size()*v.size());
			for (size_t i=0; i < u.size() * v.size(); i++)
			  D[i] = -1;
		}
		
	typename TrajectoryType::value_type::value_type  operator ()(size_t i, size_t j)
	{
		size_t index = i*cols + j;
		if (D[index] < 0)
		  D[index] = d(u[i],v[j]);
		return D[index];
	}
};




template<class TrajectoryType>
typename TrajectoryType::value_type::value_type sum_of_pairs(TrajectoryType &u, TrajectoryType &v)
{
	if (u.size() != v.size())
	   throw(std::runtime_error("sum_of_pairs distance only defined for trajectories of equal size"));

	default_element_distance<typename TrajectoryType::value_type> d;
	typename TrajectoryType::value_type::value_type m = 0;
	typename TrajectoryType::iterator it1=u.begin();
	typename TrajectoryType::iterator it2=v.begin();
	
	while (it1 != u.end()) // equal length has been tested
	{
		typename TrajectoryType::value_type::value_type v=d(*it1,*it2);
		m += v;
		it1++;
		it2++;
				
	}
	return m;
	 
}





	
//5.1.1 Discrete Frechet Distance
template <class TrajectoryType, class DistanceType=double>
class discrete_frechet_distance_impl
{
	private:
	std::vector<std::vector<DistanceType > > CA;
	element_distance<typename TrajectoryType::value_type,DistanceType> *d;	
	protected:
	
	double c(int i,int j,TrajectoryType &t1, TrajectoryType &t2	)
	{
		if (CA[i][j]>-1)		// don't update
			return CA[i][j];
    if (i==0 && j==0){
        return (CA[i][j] = (*d)(t1[0],t2[0])); 
    }
    if( i>0 && j==0){
        CA[i][j] = MAX( 
			c(i-1,0,t1,t2), 
			(*d)(t1[i],t2[0]));
        return CA[i][j];
    }
    if( i==0 && j>0){
        CA[i][j] = MAX( 
			c(0,j-1,t1,t2),
			(*d)(t1[0],t2[j]));
        return CA[i][j];
    }
    if (i>0 && j>0)
    {
	 CA[i][j] = MAX(MIN(MIN(
			c(i-1,j  ,t1,t2),
		    c(i-1,j-1,t1,t2)),
			c(i  ,j-1,t1,t2)),
			(*d)(t1[i],t2[j]));
        return CA[i][j];
    }
    CA[i][j] = std::numeric_limits<double>::infinity();
    return CA[i][j];
}

	
	
public:
	discrete_frechet_distance_impl():d(NULL){};
    DistanceType operator()(TrajectoryType u, TrajectoryType v,
    element_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {	
	  this->d = &d;
#ifdef DEBUG_DISCRETE_FRECHET
	 std::cout << "Comparing " << u.size() << " and " << v.size() << std::endl;
#endif
	DistanceType distance=0;
    if (u.size() == 0) return 0;	
    for (size_t i=0; i < u.size(); i++) 
    {
		CA.push_back(std::vector<DistanceType>());
       for (size_t j=0; j < v.size(); j++)
          CA[i].push_back(-1);
	  }
	  distance = c(u.size()-1,v.size()-1,u,v); 
	return distance;
  }
   
};      

template<class TrajectoryType>
double discrete_frechet(TrajectoryType &u, TrajectoryType &v)
{
	discrete_frechet_distance_impl<TrajectoryType, double> df;
	default_element_distance<typename TrajectoryType::value_type> d;
	return df(u,v,d);
}	


/*
 * Use a segmentation algorithm returning indizes under the constraint
 * that segments have at least three points (that is one inner point)...
 *  
 * TODO: this copies a vector... could be a filtering interator
 * */

template<class indexarray>
indexarray removeAdjacent(indexarray &ia)
{
	indexarray ret;
	// So this just creates segments. No segment is allowed to
	// end with the entry before the last one.
	for (size_t i=0; i < ia.size()-2; i++)
	{
		if(ia[i] != ia[i+1]-1)
		   ret.push_back(ia[i]);
	}
	// and adding the last one creates another true (non-empty) segment.
	ret.push_back(ia.back());
	return ret;
}




} // NAMESPACE


#endif

