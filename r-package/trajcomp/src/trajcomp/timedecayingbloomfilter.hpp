#ifndef TIMEDECAYINGBLOOMFILTER
#define TIMEDECAYINGBLOOMFILTER

#include "bloomfilter.hpp"

namespace trajcomp{
	namespace bloom {

class TimeDecayingBloomFilter : public BloomFilter<int>
{
	public:
		
		TimeDecayingBloomFilter() : BloomFilter<int>()
		{
		};
  
		TimeDecayingBloomFilter(int d, int size, int seed =-1) 
				:BloomFilter<int>(d,size,seed)
		{ 
			
		};
		
		void add(std::string elem, unsigned int slots=1)
		{
#ifdef LOG_OPS
		std::cout << "TDBF::add(" << elem<< ","<<slots<<");" << std::endl;
#endif
			for (size_t i=0; i < d; i++)
			{
				filter[hash(i,elem)] += slots;
			}
			numElements ++;
		}

		bool contains(std::string elem)
		{
#ifdef LOG_OPS
			std::cout << "TDBF::contains(" << elem <<");" << std::endl;
#endif
			for (size_t i=0; i < d; i++)
			{
				if (filter[hash(i,elem)] <= 0 ) 
					return false;
			}
			return true;
		}
		
		void decay(unsigned int slots)
		{
			for (size_t i=0; i < filter.size(); i++)
			{
				filter[i] -= slots;
				if (filter[i] < 0)
					filter[i] =0;
			}
		}
		

		
		
		
	
}; // TimeDecayingFilter

} // bloom
} // trajcomp


#endif
