/*
 * TRAJCOMP Opening Window Algorithm
 * 
 * */
 
 
 #ifndef TRAJCOMP_OPENINGWINDOW_INC
 #define TRAJCOMP_OPENINGWINDOW_INC
 
 namespace trajcomp
 {
	 
	
template<typename TrajectoryType, typename DistanceType=double>	 
class opening_window_impl {

	public:
	
	std::vector<unsigned char> booleanMap; // @REMARK: could also use bool specialization
	element_segment_distance<typename TrajectoryType::value_type, DistanceType> *d;

	
		
		
	void calculate(TrajectoryType &u, size_t start, size_t end, DistanceType max_error)
   {
	   // simplify segment and recurse
//#define DEBUG_opening_window_simplify
#ifdef DEBUG_opening_window_simplify	   
	   std::cout << start << "," << end << std::endl;
#endif
		double dmax = -std::numeric_limits<DistanceType>::infinity();
		int index = 0;
		typename TrajectoryType::value_type &s = u[start];	
		typename TrajectoryType::value_type &e = u[end];
	
		for (size_t i=start+1; i < end; i++) 
		{
			DistanceType _d = (*d)(s,e,u[i]);
			if (_d>max_error)
			{
				dmax = _d;
				index = i; 
				break;
			}
		}
			// now i is the index, where the orthogonal error got too large
			if(dmax >= max_error){
				booleanMap[index] = true;
				if(index+2 < u.size()){
					calculate(u,index, index+2,max_error);
				}
			}else{
			if(end+1 < u.size()){
				calculate(u,start, end+1,max_error);
			}
			}
      
		
	}
	
	
	 TrajectoryType operator()(TrajectoryType t, DistanceType epsilon,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
	{
		auto which = indizes(t,epsilon,d);

	TrajectoryType ret;
	for (auto w:which)
	  ret.push_back(t[w]);

    return ret;    
}

	 std::vector<size_t> indizes(TrajectoryType t, DistanceType epsilon,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {
    this->d = &d;
    
	booleanMap.resize(t.size());
	for (size_t i=0; i < t.size(); i++)
		booleanMap[i] = 0;
	booleanMap[0] =booleanMap[t.size()-1]= 1;
	
    calculate(t,0,2,epsilon);
	//simplify(u,0,u.size()-1,max_error);
	// now build the simplified trajectory
	std::vector<size_t> ret;
	for (size_t i=0; i < booleanMap.size(); i++)
	  if (booleanMap[i])
	     ret.push_back(i);

    return ret;    
  }
	
	
	
}; //class

	 
	 //INTERFACE
	 
template<class TrajectoryType>
TrajectoryType opening_window(TrajectoryType &t, double epsilon,
element_segment_distance<typename TrajectoryType::value_type,double> &d)
{
	opening_window_impl<TrajectoryType, double> dp;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}	


template<class TrajectoryType>
TrajectoryType opening_window(TrajectoryType &t, double epsilon)
{
	opening_window_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}	

template<class TrajectoryType>
std::vector<size_t> opening_window_indizes(TrajectoryType &t, double epsilon)
{
	opening_window_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	auto ret = dp.indizes(t,epsilon,d);
	return ret;
}	
	 
	 
	 
 };//namespace
 
 
 #endif
