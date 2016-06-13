#ifndef COUNTINGBLOOMFILTER
#define COUNTINGBLOOMFILTER

#include "trajcomp/bloomfilter.hpp"

namespace trajcomp{
	namespace bloom {

class CountingBloomFilter : public BloomFilter<int>
{
	public:
		
		CountingBloomFilter() : BloomFilter<int>()
		{
		};
  
		CountingBloomFilter(int d, int size, int seed =-1) 
				:BloomFilter<int>(d,size,seed)
		{ 
			
		};
		
		void add(std::string elem, unsigned int num=1)
		{
#ifdef LOG_OPS
		std::cout << "CBF::add(" << elem<< ","<<num<<");" << std::endl;
#endif
			for (size_t i=0; i < d; i++)
			{
				filter[hash(i,elem)] += num;
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
		
		unsigned int count(std::string elem)
		{
#ifdef LOG_OPS
			std::cout << "TDBF::count(" << elem <<");" << std::endl;
#endif
			int num = std::numeric_limits<int>::max();
			for (size_t i=0; i < d; i++)
			{
				num = std::min(filter[hash(i,elem)],num) ;
			}
			return num;
		}

		bool remove(std::string elem, unsigned int num=1)
		{
#ifdef LOG_OPS
		std::cout << "CBF::add(" << elem<< ","<<num<<");" << std::endl;
#endif
			bool wasable = true;
			for (size_t i=0; i < d; i++)
			{
				size_t slot = hash(i,elem);
				filter[slot] -= num;
				if (filter[slot] < 0)
				{
					wasable = false;	// Removing more than inside 
					filter[slot] =0;
				}
			}
			return wasable;
		}



		
	
}; // TimeDecayingFilter

} // bloom
} // trajcomp


#endif
