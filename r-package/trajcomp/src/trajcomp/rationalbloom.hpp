#ifndef RATIONAL_BLOOM_HPP_INC
#define RATIONAL_BLOOM_HPP_INC

#include<trajcomp/trajcomp.hpp> 	// tools::make_string
#include<trajcomp/murmur.hpp>		// murmur
#include<trajcomp/bloomfilter.hpp>	// BloomFilter
#include<trajcomp/trajcomp_stringset.hpp>  // String sets


using namespace std;
using namespace trajcomp::bloom;	// for BloomFilter
using namespace trajcomp;	// for murmur



template <class element_type> 
class RationalBloomFilter{
public:
 
  std::vector<element_type> filter;
  double d;
  double totalhashes=0;
  double totalqueries=0;
   
  RationalBloomFilter()
  {
	  configure(0,0);
  }
  
  RationalBloomFilter(double d, int size)
  { 
	resize(size); 
	this->d = d;
	
  }
  
  void configure(double d, int size)
  {
	resize(size); 
	this->d = d;
  }
  double seeded_rand01(std::string s)
  {
	  const size_t max = 1E4;
	  auto h =trajcomp::murmur::murmur(s,0)[0] % max;
	  double r = h / (double) (max-1); 
	  //~ cout << "@seededrand: " << r << endl;
	  // this is good enough.
	  return r;
  }
  
  
  
  
  void summary()
  {
	  cout << "Filter Summary:\t"<< endl;
	  cout << "Size: \t" << filter.size() << endl;
	  int ones = 0;
	  for(size_t i=0; i < filter.size(); i++)
		if (filter[i] == 1) 
			ones ++;
	  cout << "Ones: \t" << ones << endl;
  }
  void resize(int size)
  {
	  filter = std::vector<element_type>(size);
  }
  
  double foz()
  {
	  size_t o=0;
	  for (size_t i=0; i < filter.size(); i++)
	    if (filter[i] == 0)
	      o++;
	  return  (double) o / (double) filter.size();
  }

   double _esize()
  {
	  return -log(foz())*(double)filter.size()/(double) d;
  }

  unsigned int esize()
  {
	  double num = _esize();
	  return (unsigned int) num;
  }
   double _eUnion(RationalBloomFilter &b)
  {
	  //sets. Swamidass & Baldi (2007) show
	  double dot = 0;
	  for (size_t i =0; i < filter.size(); i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  // dot 
	  return -log(1- dot / (double) filter.size())*filter.size()/d;
	  
  }
  unsigned int eUnion(RationalBloomFilter &b)
  {
	    
	  return (unsigned int) _eUnion(b);
  }
  double _eIntersect(RationalBloomFilter &b)
  {
	  double dot = 0;
	  //cout << filter.size() << "==" << filter.size() << endl;
	  for (size_t i =0; i < filter.size(); i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  double union_size = _eUnion(b);
	  double e = _esize() + b._esize() - union_size;
	  
	  return e;
	  
  }
  unsigned int eIntersect (RationalBloomFilter &b)
  {
	    
	  return (unsigned int) _eIntersect(b);
  }



friend std::ostream& operator<< (std::ostream& lhs, const RationalBloomFilter & p)
{
   for (size_t i=0; i < p.filter.size(); i++)
     lhs << (unsigned int) p.filter[i] << " ";
   return lhs;
}

   void dump()
   {
	for (size_t i=0; i < filter.size(); i++) cout << filter[i] <<" ";
   cout << endl;   
   }


void add(std::string &elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::add(" << code.size() <<"," << elem<< ");" << std::endl;
	#endif
	size_t hashes_to_use = (size_t) d;
	double p = seeded_rand01(elem);
	if (p <= (d - (size_t) d))		// in case, one more
		hashes_to_use ++;
		
	totalqueries++;
	totalhashes += hashes_to_use;
	/*hashes_to_use is now an integer number with a distribution such that the mean of hashes_to_use approaches
	 * the rational number d over iterations.
	 * 
	 * */
	
    for (size_t i=0; i < hashes_to_use; i++)
    {
		size_t k = 	hash(i,elem,filter.size());
		filter[k] = 1;
    }
	

}

bool contains(std::string elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::contains(" << code.size() <<"," << elem<< ")==";
	#endif
	size_t hashes_to_use = (size_t) d;
	double p = seeded_rand01(elem);
	if (p <= (d - (size_t) d))		// in case, one more
		hashes_to_use ++;

	for (size_t i=0; i < hashes_to_use; i++)
    {
		if (filter[hash(i,elem,filter.size())] ==0) 
	         return false;
    }
    return true;
}

unsigned int hash(size_t i,  std::string &elem, size_t max)
{
	// static
	trajcomp::murmur::NonCachingMurmur hasher;
	unsigned int ret;
	ret = hasher(i,elem,max);
	return ret;
}


		bool empty() 
		{
			for (size_t i=0; i < filter.size(); i++)
			  if (filter[i] != 0)
				return false;
			return true;
		}

  // relations between BF
  bool subset(RationalBloomFilter &b)
  {
	  if (b.filter.size() != filter.size())
	     throw(std::runtime_error("Subset relation only on equal sized arrays"));
	  for (size_t i=0; i < filter.size(); i++)
	    if (filter[i] != 0)
	      if (b.filter[i] == 0)
	         return false;
	  return true;
  }



}; // class RationalBloomFilter


#endif
