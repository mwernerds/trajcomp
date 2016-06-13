#ifndef TRAJCOMP_INVERTEDINDEX_HPP
#define TRAJCOMP_INVERTEDINDEX_HPP

/*
 * Provides a trajcomp inverted index using geohash (or similar stringify functional)
 * in order to index trajectories
 * 
 * Can also generate execution plans for similarity search at least for functionals (lambdas).
 * A specialization for back_inserter is available, such that the following calls are typically good:
 * 
 * 	IndexInvertedIndex i;
	trajcomp::geohash gh;
	i.add(trajectories, [&] (std::vector<double> &d) {return gh(d[0],d[1]); } );
	
	i.dump();
	cout << "--------------------" << endl;
	i.execute_plan("wx4gp",trajectories.size(),[](size_t i) { cout << i << "-" ;});
	cout << "--------------------" << endl;
	std::vector<size_t> plan;
	i.execute_plan("wx4gp",trajectories.size(),std::back_inserter(plan));
	cout << "Stored plan: " << tools::make_string(plan) << endl;
  	cout << "--------------------" << endl;
	return 0;	

 * 
 * 
 * */

#include<map>
#include<string>
#include<iterator>
#include "trajcomp.hpp"


namespace trajcomp{
class IndexInvertedIndex
{
		public:
		std::map<std::string,std::vector<size_t>> map;
		template<typename collection, typename stringcoder>
		void add(collection &c, stringcoder stc)
		{
			for (size_t i=0; i < c.size(); i++)
			{
				for (auto it=c[i].begin(); it != c[i].end(); it++)
				{
					auto s = stc(*it);
					//map[s].insert(i);
					// @NOTE that map[s] is sorted due to linear i. This is not parallelizable 
					// without a sort afterwards
					// you can also use the snippet 1 at the end of this file.
					if (std::lower_bound(map[s].begin(),map[s].end(), i) == map[s].end())
					{
						map[s].push_back(i);
					}
				}
			}
		}		
		
		void dump()
		{
			for (auto &m:map)
			{
				std::cout << "Item: " << m.first << "==>" << tools::make_string(m.second) << endl;
			}
		}
		
		
		
		

		template<typename callable>
		void execute_plan(std::string query, size_t max, const callable &f)
		{

			auto s = map[query];
			// First run s;
			auto it = s.begin();
			for (; it != s.end(); it++)
			{
				f(*it);
			}
			auto the_end = s.end();
			for (size_t i=0; i < max; i++)
			{
				if (it != the_end)
				if (i == *it)
				{
					it++;
					continue; // skip this
				}
				f(i);
			}
		}
		
		// while this catches back_inserter
		template<class bitype>
		void execute_plan(std::string query, size_t max, std::back_insert_iterator<bitype> f)
		{
			
			execute_plan(query,max,[&](size_t i) {*f=i; f++;});
		}
		
};


};//namespace
/*
 * 
 * THIS SNIPPET CAN BE USED, IF INSERT IS NON-MONOTONOUS...
 * 
template <class Vector, class T>
void insert_into_vector(Vector& v, const T& t) {
typename
Vector::iterator i
= std::lower_bound(v.begin(), v.end(), t);
if (i == v.end() || t < *i)
v.insert(i, t);

}

 * */
#endif
