#ifndef STRINGSET_HPP_INC
#define STRINGSET_HPP_INC
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<stdlib.h>
namespace trajcomp
{
	typedef std::vector<std::string> tss;
	
	
	
	
	void ss_unique(std::vector<std::string> &S)
	{
		sort(S.begin(), S.end());
		S.erase(unique(S.begin(),S.end()),S.end());
	}
	
	std::vector<std::string> ss_union(tss &A, tss &B)
	{
		tss ret = A;
		ret.insert(ret.end(),B.begin(), B.end());
		ss_unique(ret);
		return ret;
	}
	
	std::vector<std::string> ss_intersect(tss &A, tss &B)
	{
		
		// loop through A and find in B
		tss ret;
		for (auto it = A.begin(); it != A.end(); it++)
		{
			if (std::find(B.begin(),B.end(),*it)!=B.end())
				ret.push_back(*it);
		}
		return ret;
	}
	
	
	
	
	std::vector<std::string> ss_random_nonunique(int size, int length=3, int chars=3)
	{
		std::stringstream ss;
		std::vector<std::string> ret;
		for (size_t s=0; s < size; s++)
		{
			ss.str("");
			for (size_t i=0; i < length; i++){
				unsigned char ch = (unsigned char) 'a' + (unsigned char) (rand() % chars);
			
				ss << ch;
			}
			ret.push_back(ss.str());
		}
		return ret;
	};
	std::vector<std::string> ss_random(int size, int length=3, int chars=3)
	{
		std::vector<std::string> s = ss_random_nonunique(size,length,chars);
		ss_unique(s);
		return s;
	}


}

#endif
