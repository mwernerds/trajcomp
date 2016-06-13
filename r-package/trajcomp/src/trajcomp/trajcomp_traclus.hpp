#ifndef TRACLUS_HPP_INC
#define TRACLUS_HPP_INC

#include<queue>
#include<set>
#include<unordered_map>
#include<trajcomp/trajcomp.hpp>

#define UNCLASSIFIED -1
#define NOISE -2
#define REMOVED -3


template<typename point>
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


trajcomp::default_element_distance<std::vector<double>> d_euc;

/*Projection point calculation*/
template<typename sample_type>
sample_type projection_point(sample_type &u,sample_type &v, sample_type &p) 
  {
	  
	  // the length of the segment
	  trajcomp::default_element_distance_squared<sample_type> d2;
	  
	  double l = d2(u,v);

	   if (fabs(l) < 1E-12)
	   		return u;
	   
		double t = 0;
		for (size_t i=0; i <  u.size(); i++)
		   t += ((v[i]-u[i])*(p[i]-u[i]));
		t/=l;
		
		sample_type proj;
		for (size_t i=0; i <  u.size(); i++)
		   //proj.push_back(v[i]+t*(u[i]-v[i]));
		   proj.push_back(u[i]+t*(v[i]-u[i]));
		   
		return proj;
  }



template<typename point>
double perpen_dist(point &si,point &ei,point &sj,point &ej)
{
	point ps = projection_point(si,ei,sj);
	point pe = projection_point(si,ei,ej);
	
	double dl1 = d_euc(sj,ps);
	double dl2 = d_euc(ej,pe);
	
	if (fabs(dl1+dl2) < 0.0001) 
	    return 0;
	
	return (dl1*dl1 + dl2*dl2)/(dl1 + dl2); //@TODO: use d2 
};

template<typename point>
double angle_dist(point &si,point &ei,point &sj,point &ej)
{
	double alpha1 = atan2(ei[1] - si[1],ei[0] - si[0]);
	double alpha2 = atan2(ej[1] - sj[1],ej[0] - sj[0]);
	
	double l = d_euc(sj,ej);
	
	return l * fabs(sin(alpha2-alpha1));
}

template<typename point>
double par_dist(point &si,point &ei,point &sj,point &ej)
{
	point ps = projection_point(si,ei,sj);
	point pe = projection_point(si,ei,ej);
	
	double l1 = std::min(d_euc(ps,si),d_euc(ps,ei));
	double l2 = std::min(d_euc(pe,si),d_euc(pe,ei));
	
	return std::min(l1,l2);
}


/*Total distance, which is a weighted sum of above*/
template<typename point>
double total_dist(point &si,point &ei,point &sj,point &ej, 
		double w_perpendicular=0.33, double w_parallel=0.33, double w_angle=0.33)
{
	double td =   w_perpendicular * perpen_dist(si,ei,sj,ej)
			    + w_parallel      *    par_dist(si,ei,sj,ej)
				+ w_angle         *  angle_dist(si,ei,sj,ej);
	return td;
}



template<typename iteratortype>
double MDL_PAR(iteratortype st, iteratortype en)
{
	double d1 = d_euc(*st, *en);
		
	double PER_DIST=0, ANG_DIST=0;
	iteratortype it = st;
	iteratortype it2 = it;
	it2 ++;
	
	while (true)
	{
		double d2 = d_euc(*it, *it2);
		if (d1 >= d2)
		{
			PER_DIST += perpen_dist(*st,*en,*it,*it2);
			ANG_DIST += angle_dist(*st,*en,*it,*it2);
		}else{
			PER_DIST += perpen_dist(*it,*it2,*st,*en);
			ANG_DIST += angle_dist(*it,*it2,*st,*en);
		}
		
		if (it2 == en) 		
			break;
		it ++;
		it2 ++;
	}
	
	
	// information calculation
	double LH = log2(d1);
	double LDH = log2(PER_DIST) + log2(ANG_DIST);
	
	return LH + LDH;
}

// Note: this is never (and must not be) called for an empty trajectory.
template<typename ttraj>
ttraj traclus_partition (ttraj &A)
{
	ttraj CP;
	
	CP.push_back(A[0]);  
	

    
	typename ttraj::iterator it,it2,it2_old;
	it =  A.begin();
	it2_old = it2 = it;
	it2 ++;
	
	while (it2 != A.end())
	{
		double cost = MDL_PAR(it, it2);
		double cost2 = log2(d_euc(*it,*it2));
		//cout << cost << "#" << cost2 << endl;
		//cout << "it2:" << make_string(*it2) << endl;
		if (cost > cost2 && !(fabs(cost) < 0.0001)) // right side: skip over equal points
		{
			//cout << "adding" << make_string(*it2_old);
			CP.push_back(*it2_old);
			it = it2_old;			
			it2++;
			it2_old ++;
		}else{
			it2 ++;
			it2_old ++;
		}
	}
	CP.push_back(A.back());
	return CP;
}


//double calc_dist(point &li1,point &li2,point &lj1,point &lj2, size_t key = 0)
template<typename point>
double calc_dist(std::vector<tLineSegment<point>> &L,size_t i, size_t j)
{
/*	static std::unordered_map<size_t,double> cache; //@REMARK: magic statics are not thread-safe at Microsoft
	std::unordered_map<size_t,double>::const_iterator it;
	it = cache.find(key);
	if (it != cache.end())
	{ 
		return (*it).second;
	}*/
	//@TODO: refactor
	#define li1 L[i].s
	#define li2 L[i].e
	#define lj1 L[j].s
	#define lj2 L[j].e

	double len_i = d_euc(li1,li2);
	double len_j = d_euc(lj1,lj2);
	// total_dist: larger segement first.
	double td = 
	   (len_i >= len_j)?total_dist(li1,li2,lj1,lj2):total_dist(lj1,lj2,li1,li2);

	/*if (cache.size() < 256E3/4)
		cache.insert(std::pair<size_t, double>(key, td));*/
	return td;
}


/* Compute neighborhood containment as vector<bool>
 * 
 * _d can be globally cached functional class and
 *  is called with indizes instead of objects to 
 *  facilitate caching
 * */
 
 
/*template<typename segmentcollection, typename point, typename index_distance_getter>
std::vector<bool> compute_Ne_containement(segmentcollection &D, 
					size_t index, double epsilon, index_distance_getter &_d)
{
	std::vector<bool> ret(D.size());
	
	for(size_t i = 0; i < D.size(); i++)
	{
		double d=0;
		d = _d(i,index);
		if (d <= epsilon)
			ret[index] = true;
	}
	return ret;	
};*/

/* Realize the neighborhood in memory as a segment collection
 * 
 * 
 * */

/*template<typename segmentcollection, typename point, typename index_distance_getter>
segmentcollection compute_Ne(segmentcollection &D, 
		size_t index, double epsilon, index_distance_getter &_d)
{
	std::vector<bool> result = compute_Ne_containment(D,index,epsilon,_d);

	segmentcollection ret;
	for(size_t i = 0; i < result.size(); i++)
	  if(result[i])
		ret.push_back(D[i]);
	return ret;	
};
*/

/*
 * Traclus Implementation of libTrajcomp
 * 
 * The Traclus algorithm is a classical algorithm for trajectory clustering.
 * In this C++ header, the two most important steps are implemented:
 * 		1) Creating a trajectory segment set
 * 		2) Clustering of these segments
 *      3) Removement of cluster segments, which have enoguh different 
 *         input segments but not enough different trajectories
 * 
 * Basically, traclus has two parameters and follows the density clustering
 * approach. 
 * 	 epsilon:  threshold defining nearness of trajectories
 *             Note that this epsilon is relative to a non-straightforward definition
 *             of similarity given in total_distance and calc_dist
 *
 *   minLines: Minimum number of distinct lines to make a segment a core
 *             line segment 
 * 
 * */


/*CONSTANTS FOR SPECIAL CLUSTERS*/



/**
 * Compute the epsilon-neighborhood as a vector of indizes given a line segment
 * as an index into the collection of line segments L
 * 
 * Note that the line segment idx is *not* added 
 * */
 //std::vector<tLineSegment<point>> &L,size_t i, size_t j
 
template<typename point, typename distance_functor>
std::vector<size_t> compute_NeIndizes(std::vector<tLineSegment<point>> &L,size_t idx,
		double eps, distance_functor &_d)
{
	std::vector<size_t> ret;
	#pragma omp parallel for
	for (size_t i=0; i < L.size(); i++)
	  if (idx != i)
//		if (_d(L[i].s,L[i].e,L[idx].s,L[idx].e,i*L.size() + idx) <= eps)
		if (_d(L,i,idx) <= eps)
		{
		  #pragma omp critical
	      ret.push_back(i);
	    }
	
	return ret;	
}
/**
 * Expand the epsilon-connected neighborhood of a cluster
 * using the line segment collection L, the queue of unprocessed elements
 * to be looked at, epsilon, minlines, the clusterID to be assigned and 
 * the distance functional.
 * 
 */

template<typename point,typename distance_functor>
void expandCluster(std::vector<tLineSegment<point>> &L,
				   std::queue<size_t> &Q,
					double eps,	size_t minLines, size_t clusterID,distance_functor &_d)
{
     while (!Q.empty())
     {
		  size_t m = Q.front();
		 auto Ne = compute_NeIndizes(L,m,eps,_d);
		 Ne.push_back(m);		 		 
		 if (Ne.size() >= minLines)
		 {
			for (auto it = Ne.begin(); it != Ne.end(); it++)
			{
		
				if (L[*it].cluster == UNCLASSIFIED)	
					Q.push(*it);
				if (L[*it].cluster < 0)
				    	L[*it].cluster = clusterID;
			}
		 }
		 // WAR IN l. 310, hat lorenz wieder anders gemacht
//		 if (L[m].cluster < 0)
	//					L[m].cluster = clusterID;

		 Q.pop();
	 }
}

/**
 * Perform the grouping. This is the central clustering work:
 * 
 * Finds the next UNCLASSIFIED segment, calculates the epsilon-neighborhood,
 * if it is large enough (minLines) expands the Cluster using density-connectedness and
 * the current clusterID, which is then incremented, 
 * otherwise it marks the eps-neighborhood as NOISE segments
 * 
 * @Question: something, that is in eps-distance to some segment which
 * is not a core segment will never become a core segment?
 * */
template<typename point,typename distance_functor, typename progress_visitor>
void grouping(std::vector<tLineSegment<point>> &L, double eps, size_t minLines,
				distance_functor &_d, progress_visitor &pvis)
{
	size_t clusterID=0;
	std::queue<size_t> Q;
	pvis.init(L.size());
	for (size_t i=0; i < L.size(); i++)
	{
		if (i%100==0)
		pvis(i,L.size(),"Phase 2: Clustering Segments");
		if (L[i].cluster == UNCLASSIFIED)
		{
			std::vector<size_t> Ne = compute_NeIndizes<point>(L,i,eps,_d);
			//cout << "Ne[" << i << "]="<< Ne.size() << endl;
			
			if (Ne.size()+1 >= minLines)
			{
				L[i].cluster = clusterID;
				for (auto it = Ne.begin(); it != Ne.end(); it++)
				{
					L[*it].cluster = clusterID;
					Q.push(*it);
				}
				expandCluster(L,Q, eps,minLines,clusterID,_d);
				clusterID++;
			}else  // not minLines
			{
				L[i].cluster = NOISE; 
			}
		}		
		
	}
	/*Step 3: check that clusters come from more than minLines 
	 * 		  different trajectories*/
	for (int i=0; i < (int)clusterID; i++)
	{
		std::set<size_t> sources;
		for (size_t j = 0; j < L.size(); j++)
		  if (L[j].cluster == i)
		    sources.insert (L[j].trajectory_index);
		    
		if (sources.size() < minLines)
		{
		   for (size_t j = 0; j < L.size(); j++)
		      if (L[j].cluster == i)
			     L[j].cluster = REMOVED;
		}
	}
	/*Step 3a: ClusterID compression into a range*/
	
	
	
    /* Step 4: Representation*/	
    // left out, will do that later. This step creates a representative
    // trajectory for a cluster.
    	
};


class traclus_progress_visitor
{
	public:
	void init(size_t count)
	{
		(void)count;
	}
	void operator() (size_t step, size_t count, std::string phase)
	{
	    (void)step;(void)count;(void)phase;	// suppresses unused warning
	}
	void finish()
	{
	}
	
};


/*Traclus Main Function with progress_visitor */

template<typename TrajectoryCollection, typename partitioning_functional,
typename distance_functional, typename progress_visitor>
std::vector<tLineSegment<typename TrajectoryCollection::value_type::value_type>> 
traclus(TrajectoryCollection &A, double eps, size_t minLines, partitioning_functional &part,
distance_functional &_d, progress_visitor  &pvis)
{
	
	typedef typename TrajectoryCollection::value_type TrajectoryType;
	std::vector<size_t> index_map;
	
	
	std::vector<tLineSegment<typename TrajectoryType::value_type>> segments;
	TrajectoryType L;
	size_t i=0;
	pvis.init(A.size());
	auto it = A.begin();
	while (it != A.end())
	{
		if((*it).size() == 0)		// ignore empty
			continue;
		TrajectoryType CPi;
		 CPi = part(*it);
		 for (size_t k=0; k < CPi.size() -1; k++)
		 {
			 segments.push_back(tLineSegment<typename TrajectoryType::value_type>(CPi[k],CPi[k+1],i));
		 }
		it++;i++;
		pvis(i,A.size(),"Phase 1: Segment Creation");
	}

	 
	 grouping(segments,eps,minLines,_d, pvis);
     pvis.finish();
	 return segments;
};


// and without a progress bar (i.e. an empty one)
template<typename TrajectoryCollection, typename partitioning_functional,
typename distance_functional>
std::vector<tLineSegment<typename TrajectoryCollection::value_type::value_type>>  traclus(TrajectoryCollection &A, double eps, size_t minLines, partitioning_functional &part,
distance_functional &_d)
{
   traclus_progress_visitor v;
   return traclus(A,eps,minLines,part,_d,v);	
};






#endif
