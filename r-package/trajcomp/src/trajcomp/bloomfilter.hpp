#ifndef BF_HPP_INC
#define BF_HPP_INC
#include<vector>
#include<iostream>
#include<sstream>
#include<string>
#include<unordered_map>
#include<math.h> 

#include "murmur.hpp"


using namespace std;


namespace trajcomp{
	namespace bloom {




template <class element_type> 
class BloomFilter{
public:
  std::vector<element_type> filter;
  int d;
  size_t fsize=0;
  
  unsigned int numElements;
  
  
  BloomFilter()
  {
	  configure(0,0);
  }
  
  BloomFilter(int d, int size)
  { 
	resize(size); 
	this->d = d;
	numElements = 0;
  }
  
  void configure(int d, int size)
  {
	resize(size); 
	this->d = d;
	this->fsize = size;
	numElements = 0;
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
	  numElements = 0;
	  this->fsize = size;
	  filter = std::vector<element_type>(size);
      //filter.clear(); for (size_t i =0; i < size; i++) filter.push_back(0);
  }
  
  
  
  double foz()
  {
	  size_t o=0;
	  for (size_t i=0; i < fsize; i++)
	    if (filter[i] == 0)
	      o++;
	  return  (double) o / (double) fsize;
  }

   double _esize()
  {
	  return -log(foz())*(double)fsize/(double) d;
  }

  unsigned int esize()
  {
	  double num = _esize();
	  return (unsigned int) num;
  }
   double _eUnion(BloomFilter &b)
  {
	  //sets. Swamidass & Baldi (2007) show
	  double dot = 0;
	  for (size_t i =0; i < filter.size(); i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  // dot 
	  return -log(1- dot / (double) fsize)*fsize/d;
	  
  }
  unsigned int eUnion(BloomFilter &b)
  {
	    
	  return (unsigned int) _eUnion(b);
  }
  double _eIntersect(BloomFilter &b)
  {
	  double dot = 0;
	  //cout << fsize << "==" << filter.size() << endl;
	  for (size_t i =0; i < fsize; i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  double union_size = _eUnion(b);
	  double e = _esize() + b._esize() - union_size;
	  
	  return e;
	  
  }
  unsigned int eIntersect (BloomFilter &b)
  {
	    
	  return (unsigned int) _eIntersect(b);
  }



friend std::ostream& operator<< (std::ostream& lhs, const BloomFilter & p)
{
   for (size_t i=0; i < p.filter.size(); i++)
     lhs << (unsigned int) p.filter[i] << " ";
   return lhs;
}

   void dump()
   {
	for (size_t i=0; i < filter.size(); i++) cout << (unsigned int) filter[i] <<" ";
   cout << endl;   
   }


/*unsigned int hash(int index,const char *s, char orient='+')
{
    std::string st(s);
    return hash(index,st, orient);
}

unsigned int hash(int index,const std::string &elem)
{
	unsigned int ret;
  std::stringstream stream;
  stream << index << elem << orient;
#ifdef TRAJCOMP_HASHCACHE
   std::unordered_map<std::string,unsigned int>::const_iterator 
		it = hashcache.find(stream.str());
	if (it == hashcache.end())
	{ // calculate hash and store
#endif
  std::vector<uint32_t> hash = 
	//trajcomp::murmur::murmur(stream.str(),this->gSeed);
	trajcomp::murmur::murmur(elem,index);
	ret = hash[1] % filter.size(); 
	
#ifdef TRAJCOMP_HASHCACHE
	hashcache.insert(std::pair<std::string, unsigned int>(stream.str(), ret));
	}else{
		ret = (*it).second;
	}
#endif	
  return ret;
}*/



void add_byval(std::string elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::add(" << code.size() <<"," << elem<< ");" << std::endl;
	#endif
	  for (size_t i=0; i < d; i++)
    {
//		cout << "Call: " << i << elem << filter.size() << endl;
		size_t k = 	hash(i,elem,filter.size());
		filter[k] = 1;
    }
    numElements ++;
}


void add(std::string &elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::add(" << code.size() <<"," << elem<< ");" << std::endl;
	#endif
	  for (size_t i=0; i < d; i++)
    {
//		cout << "Call: " << i << elem << filter.size() << endl;
		size_t k = 	hash(i,elem,filter.size());
		filter[k] = 1;
    }
    numElements ++;
}

bool contains(std::string &elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::contains(" << code.size() <<"," << elem<< ")==";
	#endif
	for (size_t i=0; i < d; i++)
    {
		if (filter[hash(i,elem,filter.size())] ==0) 
	         return false;
    }
    return true;
}

bool contains_byval(std::string elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::contains(" << code.size() <<"," << elem<< ")==";
	#endif
	for (size_t i=0; i < d; i++)
    {
		if (filter[hash(i,elem,filter.size())] ==0) 
	         return false;
    }
    return true;
}


unsigned int hash(size_t i,  const std::string elem, size_t max)
{
	trajcomp::murmur::NonCachingMurmur hasher;
	unsigned int ret = hasher(i,elem,max);
	return ret;
}

bool contains(const std::string &elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::contains(" << code.size() <<"," << elem<< ")==";
	#endif
	for (size_t i=0; i < d; i++)
    {
		if (filter[hash(i,elem,filter.size())] ==0) 
	         return false;
    }
    return true;
}


		bool empty() 
		{
			for (size_t i=0; i < filter.size(); i++)
			  if (filter[i] != 0)
				return false;
			return true;
		}

  // relations between BF
  bool subset(BloomFilter &b)
  {
	  if (b.filter.size() != filter.size())
	     throw(std::runtime_error("Subset relation only on equal sized arrays"));
	  for (size_t i=0; i < filter.size(); i++)
	    if (filter[i] != 0)
	      if (b.filter[i] == 0)
	         return false;
	  return true;
  }



}; // class BloomFilter



} // bloom
} // trajcomp
#endif
