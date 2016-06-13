#ifndef EDIT_DISTANCE_HPP
#define EDIT_DISTANCE_HPP

#include<assert.h>

namespace edit_distance {

	
	std::vector<std::string> explode_const(const std::string &a, size_t k) 
	{
		std::vector<std::string> ret;
		assert(a.size() % k == 0);
		for (size_t i=0; i < a.size() / k; i++)
		   ret.push_back(a.substr(i*k,k));
		return ret;
	}
	
	// ed(explode_const(A,3),explode_const(B,3))


	// uses concept [], != , size()
	template <typename in1, typename in2>
	int ed(const in1 &t1, const in2 &t2) {
		int rows = t1.size() + 1;
		int cols = t2.size() + 1;
		int matrix[rows][cols];

		// initialize first row (0,1,2,...)
		for (int i = 0; i < cols; i++) {
			matrix[0][i] = i;
		}

		// initialize first column (0,1,2,...)
		for (int i = 1; i < rows; i++) {
			matrix[i][0] = i;
		}

		// initialize rest of the matrix (0....0)
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				matrix[i][j] = 0;
			}
		}

		int delt;
		// fill the matrix with "distances"
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				if (t1[i - 1] != t2[j - 1]) {
					delt = 1;
				}
				else delt = 0;
				matrix[i][j] = std::min(matrix[i - 1][j - 1] + delt, 
								std::min(matrix[i - 1][ j] + 1, matrix[i][ j - 1] + 1));
			}
		}
		return matrix[rows-1][cols-1];
}
	
	
	

	// edit distance for strings (dynamic programming)
	/*int ed(const std::string &t1, const std::string &t2) {
		int rows = t1.size() + 1;
		int cols = t2.size() + 1;
		int matrix[rows][cols];

		// initialize first row (0,1,2,...)
		for (int i = 0; i < cols; i++) {
			matrix[0][i] = i;
		}

		// initialize first column (0,1,2,...)
		for (int i = 1; i < rows; i++) {
			matrix[i][0] = i;
		}

		// initialize rest of the matrix (0....0)
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				matrix[i][j] = 0;
			}
		}

		int delt;
		// fill the matrix with "distances"
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				if (t1[i - 1] != t2[j - 1]) {
					delt = 1;
				}
				else delt = 0;
				matrix[i][j] = std::min(matrix[i - 1][j - 1] + delt, 
								std::min(matrix[i - 1][ j] + 1, matrix[i][ j - 1] + 1));
			}
		}
		return matrix[rows-1][cols-1];
}*/

}

#endif
