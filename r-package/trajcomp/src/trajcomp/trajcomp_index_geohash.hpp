#ifndef GEOHASH_BLOOM_INDEX
#define GEOHASH_BLOOM_INDEX

#include<trajcomp/trajcomp.hpp>
#include<iostream>
#include<vector>
#include <trajcomp/trajcomp_geohash.hpp>
#include <trajcomp/bloomfilter.hpp>

#include<algorithm>
#include<string>

using namespace std;


namespace trajcomp
{
	
	// Forward declared for friendship
template<class Trajectory, class BF>
class CorrelationTree;
	

		
/*Geohash Index of Trajectory*/		
template<class Trajectory, class BF, typename bloom_d_type = size_t>
class geohash_bloom_index
{
	
	
	private:
	std::vector<BF> ftracks; // filter tracks
	std::vector<size_t> N;   // the number of elements of each filter
	bloom_d_type cfg_bloom_d;
	size_t cfg_bloom_size;
	size_t cfg_seed;
	size_t cfg_geohash_length;
	
	std::vector<std::vector<double>> foz_histograms;
	size_t foz_k;
	
	public:
	typedef BF query_type;
	geohash_bloom_index():cfg_bloom_d(5),cfg_bloom_size(1024),cfg_seed(123), cfg_geohash_length(8)
	{
	}
	geohash_bloom_index(bloom_d_type d,size_t size,size_t geohash_length, size_t seed=0):cfg_bloom_d(d),cfg_bloom_size(size),cfg_seed(seed), cfg_geohash_length(geohash_length)
	{
		
	};
	
	double foz(size_t k)
    {
		return ftracks[k].foz();
	}	

	void generate(std::vector<Trajectory> &data)
	{
		cout << "Generate using " << cfg_geohash_length <<"on" << data.size() <<endl;
		geohash gh;
	
		ftracks.resize(data.size());
		N.resize(data.size());
		auto ft = ftracks.begin();
		auto itN = N.begin();
		size_t i=0; 
		for (auto it = data.begin(); it != data.end(); it++)
		{
			BF bf(cfg_bloom_d,cfg_bloom_size); // Remark: there is no seed anymore in libtrajcomp ,cfg_seed);
			tools::progress(i++,data.size(),"Creating Index");

			std::vector<std::string> set = trajectory2geohashset<Trajectory>(*it,cfg_geohash_length);
			for(size_t i=0; i< set.size(); i++)
			   bf.add(set[i]);
			*ft = bf;
			*itN = set.size();
			itN ++;
			ft ++;
		}
		cout << endl;
	}
	
	query_type empty_query()
	{
		query_type bf(cfg_bloom_d,cfg_bloom_size,cfg_seed);
		return bf;
	}
	
	query_type trajectory_query(Trajectory &t)
	{
		query_type bf(cfg_bloom_d,cfg_bloom_size);//,cfg_seed);
		auto set = trajectory2geohashset(t,cfg_geohash_length);
		for (auto it = set.begin(); it != set.end(); it++)
		  bf.add(*it);
		return bf;
	}
	
	
	
	double LB(std::vector<string> &S1, BF &f2, size_t N2)
	{
		//cout << "DEBUG: N2=" << N2 << endl;
		 // phi1 = #Elem of A, die auf B passen
		double phi1 = 0;
		for (size_t i=0; i < S1.size(); i++)
			if (f2.contains(S1[i]))
				phi1 ++;
		// phi2 = #Elem of A, die nicht auf B passen + #B
		double phi2 = N2;
		for (size_t i=0; i < S1.size(); i++)
			if (!f2.contains(S1[i]))
				phi2 ++;
		double LB = 1 - phi1/phi2;
		return LB;
	}
	double LB(std::vector<string> &S1, size_t k)
	{
		return LB(S1,ftracks[k],N[k]);
	}
	double LB(Trajectory t, size_t k)
	{
		auto set = trajectory2geohashset(t,cfg_geohash_length);
		return LB(set,k);		
	}
	

	
	double jaccard(query_type &q, size_t k)
	{
		double jac = ( 1 - ( ftracks[k]._eIntersect(q)/ ftracks[k]._eUnion(q)));
		if (std::isnan(jac))
		{
			jac = 1;
		}
		return jac;
	}
	
	double intersect(query_type &q, size_t k)
	{
		return ftracks[k]._eIntersect(q);
	}
	
	double relative_subset_contradicts(query_type &q, size_t k)
	{
		double query_ones=0;	
		double contradictions=0;
		//cout << "query" << q.filter.size() <<";"<< q.d <<"," << q.gSeed << endl;
		//cout << "k" <<  k<<" <ftracksize" << ftracks.size() << endl;
		//cout << "match" << ftracks[k].filter.size() <<"," << ftracks[k].d <<"," << ftracks[k].gSeed << endl;
		for (size_t i=0; i<q.filter.size(); i++)
		{
			 if (q.filter[i] == 1)
			 {
				 query_ones ++;
				 if (ftracks[k].filter[i] != 1)
					contradictions ++;
			 }
		}
		//if (contradictions != query_ones)
			//std::cout << "contrad/query: " << contradictions << "/"<< query_ones << std::endl;
		return contradictions / query_ones;
	}
	
	
	
	
	
	size_t bytes()
	{
		size_t ret = 0;
		for (size_t i=0; i < ftracks.size(); i++)
		{
			ret += cfg_bloom_size / 8;
		}
		return ret;
	}
		
	// fast histogram query
	std::vector<double> createFozHistogram(query_type &q, size_t k)
	{
		if (q.filter.size() % k != 0)
		  throw(std::runtime_error("createFozHistogram: Slot number must divide the hash size."));
		std::vector<double> foz_histogram(k);
		size_t slice = q.filter.size() / k;
		for (size_t i=0; i < k; i++)
		{
			size_t low = i*slice;
			size_t high = (i+1)*slice;
			double sum = 0;
			for (size_t j=low; j < high; j++)
				sum += q.filter[j]?1:0;
			foz_histogram[i] = sum / slice;
		}
		return foz_histogram;
	}
	
	void createFozHistograms(size_t k)
	{
		cout << "Generating Histograms: " << k << endl;
		foz_k = k;
		foz_histograms.clear();
		auto ft = ftracks.begin();
		size_t i=0; 
		for (auto it = ftracks.begin(); it != ftracks.end(); it++)
		{
			foz_histograms.push_back(createFozHistogram(*it,k));
			//tools::progress(i++,ft.size(),"Creating Index");
		}
		cout << "Done." << endl;
	}
	
	std::vector<size_t> fast_subset(query_type &q)
	{
		std::vector<size_t> candidates;
		std::vector<double> qh = createFozHistogram(q, foz_k);
		for (size_t i=0; i<  foz_histograms.size(); i++)
		{
			bool cand = true;
			for (size_t j=0; j < foz_k; j++)
			{
				
				if (qh[j] > foz_histograms[i][j])
				{
					cand = false;
					break;	// conclude that i is not a candidate
				}
			}
			if (cand) 
				candidates.push_back(i);
		}
		return candidates;
	}
	
	friend class CorrelationTree<Trajectory,BF>;
	
};



template<class Trajectory>
class geohash_reference_index
{
	typedef std::vector<std::string> query_type;
	
	public:
	std::vector<std::vector<std::string>> ftracks; // filter tracks
	size_t cfg_geohash_length;
	
	
	
	
	public:
	geohash_reference_index(size_t geohash_length): cfg_geohash_length(geohash_length)
	{
		
	};
	
	
	double getsize(size_t k)
	{
		return ftracks[k].size();
	}
	
	void generate(std::vector<Trajectory> &data)
	{
		ftracks.resize(data.size());
		auto ft = ftracks.begin();
		size_t i=0; 
		for (auto it = data.begin(); it != data.end(); it++)
		{
			query_type s;
			tools::progress(i++,data.size(),"Creating Reference Index");
			*ft = trajectory2geohashset(*it,cfg_geohash_length);
			ft ++;
		}
		cout << endl;
	}
	
	query_type trajectory_query(Trajectory &t)
	{
		return trajectory2geohashset(t,cfg_geohash_length);
	}
	
	double jaccard(query_type &q, query_type &b)
	{
		
		/*double num_equal = _intersect(q,b);
		double num_total = (q.size() + b.size() - num_equal);
	
		return (1 - (num_equal / num_total));*/
		double U,I;
		auto A = q;
		auto B = b;
			sort(A.begin(), A.end());
			sort(B.begin(), B.end());
			auto alast = std::unique(A.begin(),A.end());
			auto blast = std::unique(B.begin(),B.end());
			A.erase(alast, A.end());
			B.erase(blast, B.end());
			// count elements of the union
			std::vector<std::string> u =  A;
			u.insert(u.end(),B.begin(), B.end());
			sort(u.begin(),u.end());
			auto ulast = std::unique(u.begin(),u.end());
			u.erase(ulast,u.end());
			U = u.size();
			I = A.size() + B.size() - U;
			return 1-(I/U);
	}
	
	double jaccard(query_type &q, size_t k)
	{
		return jaccard(q,ftracks[k]);
	}
	
	double intersect(query_type &q, size_t k)
	{
		return _intersect(q,ftracks[k]);
	}
	
	inline double _intersect(query_type &q, query_type &b)
	{
		
		double num_equal=0;
		for (size_t i=0; i < q.size(); i++)
			for (size_t j=0; j < b.size(); j++)
				if (q[i] == b[j])
					num_equal ++;
		return (num_equal);
	}

	size_t bytes()
	{
		size_t ret = 0;
		for (size_t i=0; i < ftracks.size(); i++)
		{
			for (size_t j=0; j < ftracks[i].size(); j++)
			   ret += ftracks[i][j].size();
		}
		return ret;
	}
	
	
	
};

	/*Correlation Tree*/
	
template<class Trajectory, class BF>
class CorrelationTree
{
	public:
	std::vector<CorrelationTree> children;
		
	CorrelationTree()
	{
		cout << "Building Correlation Tree on GHindex" << endl;
	}
	
	size_t HD(BF a, BF b)
	{
		if (a.filter.size() != b.filter.size())
		   throw(std::runtime_error("Hamming Distance only defined between equal length vectors"));
		size_t ret = 0;
		for (size_t i = 0; i < a.filter.size(); i++)
		  if (a.filter[i] != b.filter[i])
		     ret ++;
		return ret;
	}
	
	int _buildTree(std::vector<BF> &bfs, size_t k)
	{
		/*Build all children and call each */
		std::vector<BF> classes;
		std::vector<int> class_assignments;
		class_assignments.resize(bfs.size());
		for (size_t i=0; i < class_assignments.size(); i++)
		   class_assignments[i] = -1;
		for (size_t i=0; i< bfs.size(); i++)
		{
			// find class of bfs[i]
			for (size_t j=0; j < classes.size(); j++)
			{
				if (HD(classes[j],bfs[i]) < k)
				{
					class_assignments[i] = j; 
					break; 
				}
			}
			if (class_assignments[i] < 0)
			{
				// create new class
				classes.push_back(bfs[i]);
				class_assignments[i] = classes.size() -1;
			}
		}
		cout << tools::make_string(class_assignments);
		
		
		return 0;
	}
	
	int build(geohash_bloom_index<Trajectory,BF> &idx)
	{
		return _buildTree(idx.ftracks,25);
	}
	
	
	
};



}

#endif
