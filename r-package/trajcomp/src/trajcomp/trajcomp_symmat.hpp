#ifndef TRAJCOMP_SYMMAT_INC
#define TRAJCOMP_SYMMAT_INC


#include<vector>
/*a nice trick to make rvalues to lvalues*/

template <typename T>
T& as_lvalue(T&& x)
{
    return x;
}

/*Symmetric Square Matrix storing only Upper Triangular, diagonal is fixed zero
 * 
 * Especially useful for distance matrices. * 
 * */

template<class tvalue> 
class UTSquareMatrix
{
	public:
	size_t _size;
	std::vector<tvalue> data;
	
	size_t size(void)
	{
		
		return _size;
	}
	
	void resize(size_t m)
	{
		_size = m;
		size_t size = m*(m-1)/2;
		data.resize(size);
	}

	struct tsymmetricgetter{
		tvalue temp; // for making a zero reference-accessible
		std::vector<tvalue> &base;
		size_t i;
		tsymmetricgetter (std::vector<tvalue> &_base, size_t _i): base(_base),i(_i)
		{};
		tvalue &operator[] (size_t j)
		{
			if (i==j)
			{
				temp =0;
			   return  temp;
		    }
		    if (i < j )
		      swap(i,j);
		    		    
			return as_lvalue (base[((i*(i-1))/2) + (j)]);
		}
	};

	tsymmetricgetter operator[](size_t i)
	{
		return tsymmetricgetter(data, i);
	}
	
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & _size;
		ar & data;
    }
	
	
};


/* Test Main <=> Documentation ;-)
 * 
 * // #include<iostream>
#include <fstream>
//include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace std;
int main(void)
{
	UTSquareMatrix<int> m;
	const int s = 6;
	m.resize(s);
	// Fill in matrix
	
	for (size_t i=0; i < s; i++)
	{
		for (size_t j=0; j < s; j++)
		{
			m[i][j] = i+j;
			
    	}
	}
	// Dump before save
	for (size_t i=0; i < s; i++)
	{
		for (size_t j=0; j < s; j++)
		{
			cout << m[i][j]<< "\t";
		}
		cout << endl;
	}
	// Serialize
	 std::ofstream ofs("serial.dat");

    // save data to archive
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << m;
    	// archive and stream closed when destructors are called
    }
	for (size_t i=0; i < s; i++)
	{
		for (size_t j=0; j < s; j++)
		{
			m[i][j] = 4;
			
    	}
	}

for (size_t i=0; i < s; i++)
	{
		for (size_t j=0; j < s; j++)
		{
			cout << m[i][j]<< "\t";
		}
		cout << endl;
	}
	

    UTSquareMatrix<int> n;
    
    {
        // create and open an archive for input
        std::ifstream ifs("serial.dat");
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> n;
        // archive and stream closed when destructors are called
    }
    cout << "--" << endl;
	for (size_t i=0; i < s; i++)
	{
		for (size_t j=0; j < s; j++)
		{
			cout << n[i][j]<< "\t";
		}
		cout << endl;
	}
	



	
}

*/ 
#endif
