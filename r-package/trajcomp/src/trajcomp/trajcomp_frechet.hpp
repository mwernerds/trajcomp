#ifndef FULL_FRECHET
#define FULL_FRECHET


#include<algorithm>






/*This one is a "natural" string version of a NxMx2 array.
 * Below, however, octave style output is given
 * 
 * template<class datatype, char delim=' '>
std::string make_string3(datatype data)
{
		std::stringstream ss;
		typename datatype::iterator it;
		ss << std::endl;
		for(it=data.begin(); it != data.end(); it++)
		{
		  typename datatype::value_type::iterator it2;
		  for (it2 = (*it).begin();
			   it2 != (*it).end(); it2++){
				  
				ss << "\t(" << (*it2)[0]<<";"<<(*it2)[1]<<")" << delim;
			}ss  << std::endl;
		}
			return ss.str();
	}

*/

// the ugly octave way of viewing NxMx2
template<class datatype, char delim=' '>
std::string make_string3(datatype data)
{
		std::stringstream ss;
		typename datatype::iterator it;
		ss << "[:,:,1]" << std::endl;
		for(it=data.begin(); it != data.end(); it++)
		{
		  typename datatype::value_type::iterator it2;
		  for (it2 = (*it).begin();
			   it2 != (*it).end(); it2++){
				  
				ss << "\t" << (*it2)[0] << delim;
			}ss  << std::endl;
		}
		ss << "[:,:,2]" << std::endl;
		
		for(it=data.begin(); it != data.end(); it++)
		{
		  typename datatype::value_type::iterator it2;
		  for (it2 = (*it).begin();
			   it2 != (*it).end(); it2++){
				  
				  ss << "\t" << (*it2)[1] << delim;
				
			}ss  << std::endl;
		}
		return ss.str();
	}


/*
 * Following http://www.mathworks.com/matlabcentral/fileexchange/38714-frechet-distance-calculation
 * by
 * Copyright (c) 2012, Richard
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

 * */


namespace trajcomp{

//5.1.X Frechet Distance

template<class datatype> 
void matrix_resize(std::vector<std::vector< datatype> > &m, int r, int c)
{
	m.resize(r);
	for(size_t i=0; i <r; i++)
	  m[i].resize(c);	
}
template<class datatype> 
void matrix_resize(std::vector<std::vector<std::vector< datatype > > > &m, int r, int c, int d)
{
	m.resize(r);
	for(size_t i=0; i <r; i++)
	{
	  m[i].resize(c);	
	  for (size_t j=0; j < c; j++)
	    m[i][j].resize(d);
    }
    
}



template <class TrajectoryType, class DistanceType=double>
class frechet_distance_impl
{
	private:
	element_distance<typename TrajectoryType::value_type> *d2;

	std::vector<double> lP, lQ;
	std::vector<std::vector<double>> lPQ;
	double max_lPQ, min_lPQ;
	std::vector<std::vector<double>> bP;
	std::vector<std::vector<double>> bQ;
	
	std::vector<std::vector<double> > A,B,C,D;
	std::vector<std::vector<std::vector<double> > > LF,BF,LR,BR;
	
			
	public:
		frechet_distance_impl(){};

		double getMaxPointDistance(){return max_lPQ;};
		bool decide(TrajectoryType &P,TrajectoryType &Q, double epsilon, bool dump_cells = false)
		{
			size_t i,j;	
			const double NA = -1; // instead of NA
			
			

		
		
		// 2. Step Solve line segment circle intersections
		
		double xp,yp,xp1,yp1;
		double xq,yq,xq1,yq1;
		double ap, aq,bp,bq,c,dp,dq, cp ,cq;
		double up, um;
		//cout << P.size() << "?" << Q.size() << endl;
		
		for ( i=0; i < P.size()-1; i++)
		{
				xp = P[i][0];yp = P[i][1];
				xp1 = P[i+1][0];yp1 = P[i+1][1];
				for( j=0; j < Q.size()-1; j++)
				{
					xq = Q[j][0];yq = Q[j][1];
					xq1 = Q[j+1][0];yq1 = Q[j+1][1];
					// now solve intersection
					ap = lP[i];  aq = lQ[j];
					bp = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
					bq = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
					c = lPQ[i][j] - epsilon*epsilon;
					dp = bp*bp - 4*ap*c;
					dq = bq*bq - 4*aq*c;
			//		cout << dp << "#" << dq << endl;
					if (dp < 0)  //%--  a_ij, b_ij, LF_ij
					{
						A[i][j] = NA; B[i][j] = NA; 
						LF[i][j][0] = NA; LF[i][j][1] = NA;
					}else{
						up = (-bp+sqrt(dp))/(2*ap); um = (-bp-sqrt(dp))/(2*ap);
						//cout << up << ";" << um << endl;
						if (((up<0)&&(um<0))||((up>1)&&(um>1)))
						{
							//%--line segment outside circle
							A[i][j] = NA; B[i][j] = NA;
							LF[i][j][0] = NA; LF[i][j][1] = NA;
						}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
							//%--line segment is interior to circle
							A[i][j] = 0; B[i][j] = 1;
							LF[i][j][0] = 0; LF[i][j][1] = 1;
						}
						else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
						//%--one intersection (b_i,j)
						A[i][j] = 0; B[i][j] = std::max(um,up);
						LF[i][j][0] = 0; LF[i][j][1] = B[i][j];
						}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
							//%--one intersection (a_i,j)
							A[i][j] = std::min(um,up); B[i][j] = 1;
							LF[i][j][0] = A[i][j]; LF[i][j][1] = 1;
						}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
							//%--two intersections
						A[i][j] = std::min(um,up); B[i][j] = std::max(um,up);
						LF[i][j][0] = A[i][j]; LF[i][j][1] = B[i][j];
						}else
						{
							throw std::runtime_error("Unexpected case in frechet_decide at LF");
						}
				
					}// dp < 0
					if (dq < 0)  //%--   %--  c_ij, d_ij, BF_ij
					{
						C[i][j] = NA; D[i][j] = NA; 
						BF[i][j][0] = NA; BF[i][j][1] = NA;
					}else{
						up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
						if (((up<0)&&(um<0))||((up>1)&&(um>1)))
						{
							//%--line segment outside circle
							C[i][j] = NA; D[i][j] = NA;
							BF[i][j][0] = NA; BF[i][j][1] = NA;
						}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
							//%--line segment is interior to circle
							C[i][j] = 0; D[i][j] = 1;
							BF[i][j][0] = 0; BF[i][j][1] = 1;
						}
						else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
						//%--one intersection d_i,j)
						C[i][j] = 0; D[i][j] = std::max(um,up);
						BF[i][j][0] = 0; BF[i][j][1] = D[i][j];
						}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
							//%--one intersection (a_i,j)
							C[i][j] = std::min(um,up); D[i][j] = 1;
							BF[i][j][0] = C[i][j]; BF[i][j][1] = 1;
						}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
							//%--two intersections
						C[i][j] = std::min(um,up); D[i][j] = std::max(um,up);
						BF[i][j][0] = C[i][j]; BF[i][j][1] = D[i][j];
						}else
						{
							throw std::runtime_error("Unexpected case in frechet_decide at LF");
						}
				
					}// dq < 0
				}
		}// for
		
		/*std::cout << "A" << make_string2(A) << std::endl;
		std::cout << "B" << make_string2(B) << std::endl;

		std::cout << "LF" << make_string3(LF) << std::endl;
		std::cout << "BF" << make_string3(LR) << std::endl;
		*/
		
		// 3. Step: Top Row Analysis
		
		i = P.size() -1; // to have the if then else equal to above
		xp = P[P.size()-1][0];yp = P[P.size()-1][1];
		
		for ( j=0;j < Q.size()-1; j++)
		{
			xq  = Q[j][0];yq = Q[j][1];
			xq1 = Q[j+1][0];yq1 = Q[j+1][1];
			
			aq = lQ[j];
			bq = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
			c = lPQ[P.size()-1][j] - epsilon*epsilon;
			dq = bq*bq - 4*aq*c;
			
			if (dq < 0) //c_ij, d_ij, BF_ij
			{
				C[P.size()-1][j] = NA; D[P.size()-1][j] = NA;
				BF[P.size()-1][j][0] = NA; BF[P.size()-1][j][1] = NA;
			}
			else{
				
				up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
				//cout << up << "$" << um << endl;
				if (((up<0)&&(um<0))||((up>1)&&(um>1)))
				{
					//%--line segment outside circle
					C[i][j] = NA; D[i][j] = NA;
					BF[i][j][0] = NA; BF[i][j][1] = NA;
				}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
					C[i][j] = 0; D[i][j] = 1;
					BF[i][j][0] = 0; BF[i][j][1] = 1;
					}
				else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
				//%--one intersection d_i,j)
				C[i][j] = 0; D[i][j] = std::max(um,up);
				BF[i][j][0] = 0; BF[i][j][1] = D[i][j];
				}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
					//%--one intersection (a_i,j)
					C[i][j] = std::min(um,up); D[i][j] = 1;
					BF[i][j][0] = C[i][j]; BF[i][j][1] = 1;
				}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
					//%--two intersections
				C[i][j] = std::min(um,up); D[i][j] = std::max(um,up);
				BF[i][j][0] = C[i][j]; BF[i][j][1] = D[i][j];
				}else
				{
					throw std::runtime_error("Unexpected case in frechet_decide at LF");
				}
			}//else
    		
    		
		}// for j
	//	std::cout << "C" << make_string2(C) << std::endl;
	//	std::cout << "D" << make_string2(D) << std::endl;
	//	std::cout << "BF" << make_string3(BF) << std::endl;

			// 4. Step: Right Column Analysis	
			
		
		j = Q.size() -1;
		xq = Q[j][0];yq = Q[j][1];
		for (i=0; i < P.size()-1; i++) 
		{
			xp = P[i][0];yp = P[i][1];
			xp1 = P[i+1][0];yp1 = P[i+1][1];
			// and again line segment circle intersection
			ap = lP[i]; 
			bp = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
			c = lPQ[i][j] - epsilon*epsilon;
			dp = bp*bp - 4*ap*c;
			//cout << "right col dp " << dp << endl;
    		if (dp < 0)  //%--  a_ij, b_ij, LF_ij
			{
				A[i][j] = NA; B[i][j] = NA; 
				LF[i][j][0] = NA; LF[i][j][1] = NA;
			}else{
				up = (-bp+sqrt(dp))/(2*ap); um = (-bp-sqrt(dp))/(2*ap);
				//cout << up << "&%" << um << endl;
				if (((up<0)&&(um<0))||((up>1)&&(um>1)))
				{
					//%--line segment outside circle
					A[i][j] = NA; B[i][j] = NA;
					LF[i][j][0] = NA; LF[i][j][1] = NA;
				}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
					//%--line segment is interior to circle
					A[i][j] = 0; B[i][j] = 1;
					LF[i][j][0] = 0; LF[i][j][1] = 1;
				}else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
					//%--one intersection (b_i,j)
					A[i][j] = 0; B[i][j] = std::max(um,up);
					LF[i][j][0] = 0; LF[i][j][1] = B[i][j];
				}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
					//%--one intersection (a_i,j)
					A[i][j] = std::min(um,up); B[i][j] = 1;
					LF[i][j][0] = A[i][j]; LF[i][j][1] = 1;
				}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
					//%--two intersections
					A[i][j] = std::min(um,up); B[i][j] = std::max(um,up);
					LF[i][j][0] = A[i][j]; LF[i][j][1] = B[i][j];
				}else
				{
					throw std::runtime_error("Unexpected case in frechet_decide at LF");
				}
				
			}// dp < 0
			
		}
		
		
		// 5. Step: Top Right Cell
		i = P.size()-1;
		j= Q.size() -1;
		xp = P[P.size()-2][0];yp = P[P.size()-2][1];
		xp1 = P[P.size()-1][0];yp1 = P[P.size()-1][1];
		xq = Q[Q.size()-2][0];yq = Q[Q.size()-2][1];
		xq1 = Q[Q.size()-1][0];yq1 = Q[Q.size()-1][1];
		ap = lP[P.size()-2];aq = lQ[Q.size()-2];
		bp = 2*((xp-xq1)*(xp1-xp)+(yp-yq1)*(yp1-yp));
		bq = 2*((xq-xp1)*(xq1-xq)+(yq-yp1)*(yq1-yq));
		cp = lPQ[P.size()-2][Q.size()-1] - epsilon*epsilon;
		cq = lPQ[P.size()-1][Q.size()-1] - epsilon*epsilon;
		dp = bp*bp - 4*ap*cp;
		dq = bq*bq - 4*aq*cq;
		
		
		if (dp < 0)  //%--  a_ij, b_ij, LF_ij
		{
			A[i][j] = NA; B[i][j] = NA; 
			LF[i][j][0] = NA; LF[i][j][1] = NA;
		}else{
			up = (-bp+sqrt(dp))/(2*ap); um = (-bp-sqrt(dp))/(2*ap);
			if (((up<0)&&(um<0))||((up>1)&&(um>1)))
			{
				//%--line segment outside circle
				A[i][j] = NA; B[i][j] = NA;
				LF[i][j][0] = NA; LF[i][j][1] = NA;
			}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
				//%--line segment is interior to circle
				A[i][j] = 0; B[i][j] = 1;
				LF[i][j][0] = 0; LF[i][j][1] = 1;
			}else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
				//%--one intersection (b_i,j)
				A[i][j] = 0; B[i][j] = std::max(um,up);
				LF[i][j][0] = 0; LF[i][j][1] = B[i][j];
			}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
				//%--one intersection (a_i,j)
				A[i][j] = std::min(um,up); B[i][j] = 1;
				LF[i][j][0] = A[i][j]; LF[i][j][1] = 1;
			}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
				//%--two intersections
				A[i][j] = std::min(um,up); B[i][j] = std::max(um,up);
				LF[i][j][0] = A[i][j]; LF[i][j][1] = B[i][j];
			}else
			{
				throw std::runtime_error("Unexpected case in frechet_decide at LF");
			}
				
		}// dp < 0
		
		
		if (dq < 0) //c_ij, d_ij, BF_ij
		{
			C[i][j] = NA; D[i][j] = NA;
			BF[i][j][0] = NA; BF[i][j][1] = NA;
		}
		else{
			up = (-bq+sqrt(dq))/(2*aq); um = (-bq-sqrt(dq))/(2*aq);
			if (((up<0)&&(um<0))||((up>1)&&(um>1)))
			{
				//%--line segment outside circle
				C[i][j] = NA; D[i][j] = NA;
				BF[i][j][0] = NA; BF[i][j][1] = NA;
			}else if ((std::min(um,up)<0)&&(std::max(um,up)>1)){
				//%--line segment is interior to circle
				C[i][j] = 0; D[i][j] = 1;
				BF[i][j][0] = 0; BF[i][j][1] = 1;
			}
			else if ((((std::min(um,up))<=0)&&(std::max(um,up))<=1)){
				//%--one intersection d_i,j)
				C[i][j] = 0; D[i][j] = std::max(um,up);
				BF[i][j][0] = 0; BF[i][j][1] = D[i][j];
			}else if((std::min(um,up)>=0)&&(std::max(um,up)>1)){
				//%--one intersection (a_i,j)
				C[i][j] = std::min(um,up); D[i][j] = 1;
				BF[i][j][0] = C[i][j]; BF[i][j][1] = 1;
			}else if ((std::min(um,up)>=0)&&(std::max(um,up)<=1)){
				//%--two intersections
				C[i][j] = std::min(um,up); D[i][j] = std::max(um,up);
				BF[i][j][0] = C[i][j]; BF[i][j][1] = D[i][j];
			}else
			{
				throw std::runtime_error("Unexpected case in frechet_decide at LF");
			}
		}//dq<0
		
		/*cout << "##############'" << endl;
		cout << make_string2(B) << endl;
		cout << make_string3(LF) << endl;
		cout << make_string3(BF) << endl;
		cout << "#################" << endl;
		*/
		
		
		/*std::cout << "A" << make_string2(A) << std::endl;
		std::cout << "B" << make_string2(B) << std::endl;
		std::cout << "C" << make_string2(C) << std::endl;
		std::cout << "D" << make_string2(D) << std::endl;
		
		std::cout << "LF" << make_string3(LF) << std::endl;
		std::cout << "BF" << make_string3(BF) << std::endl;
		*/
				
		
	// Now the free space diagram is complete
	
	// Step 6: Compute reachable sets for each cell
	
	for (i=1; i < P.size(); i++)
		LR[i][0][0] = LR[i][0][1] = NA; 
	for (j=1; j < Q.size(); j++)
		BR[0][j][0] = BR[0][j][1] = NA; 

	for (i=0; i < P.size()-1; i++)
	  for(j=0; j < Q.size()-1; j++)
	  {
	    if (i==0 && j==0)
	    {
           if ((LF[i][j][0]==0)&&(BF[i][j][0])==0){ //   %--start at the origin
                LR[i][j+1][0] = LF[i][j+1][0]; LR[i][j+1][1] = LF[i][j+1][1];
                BR[j+1][i][0] = BF[j+1][i][0]; BR[j+1][i][1] = BF[j+1][i][1];
           }else{
                LR[i][j+1][0] = NA; LR[i][j+1][1] = NA;
                BR[j+1][i][0] = NA; BR[j+1][i][1] = NA;
			}
			
		}else // if i==0 && j ==0
		{
			if (LR[i][j][0]==NA && BR[i][j][0]==NA)
			{
				LR[i][j+1][0]=NA; LR[i][j+1][1] = NA;
				BR[i+1][j][0]=NA; BR[i+1][j][1] = NA;
			}else if (LR[i][j][0]!=NA && BR[i][j][0]==NA){
				BR[i+1][j][0] = BF[i+1][j][0];BR[i+1][j][1] = BF[i+1][j][1];
                if ((LF[i][j+1][1]<LR[i][j][0])||(LF[i][j+1][1]==NA)){
                    LR[i][j+1][0] = NA; LR[i][j+1][1] = NA;
				}else if (LF[i][j+1][0]>LR[i][j][0]){
                    LR[i][j+1][0] = LF[i][j+1][0];
                    LR[i][j+1][1] = LF[i][j+1][1];
				}
                else{
					LR[i][j+1][0] = LR[i][j][0];
                    LR[i][j+1][1] = LF[i][j+1][1];
				}
                
			}else if (LR[i][j][0] == NA && BR[i][j][0] != NA){
                 LR[i][j+1][0] = LF[i][j+1][0]; LR[i][j+1][1] = LF[i][j+1][1];
                if (BF[i+1][j][1]<BR[i][j][0]|| BF[i+1][j][1]==NA){
                    BR[i+1][j][0] = NA; BR[i+1][j][1] = NA;
                }else if (BF[i+1][j][0]>BR[i][j][0]){
                    BR[i+1][j][0] = BF[i+1][j][0];
                    BR[i+1][j][1] = BF[i+1][j][1];
				}else{
                    BR[i+1][j][0] = BR[i][j][0];
                    BR[i+1][j][1] = BF[i+1][j][1];
                }
			}else{  
                LR[i][j+1][0] = LF[i][j+1][0]; LR[i][j+1][1] = LF[i][j+1][1];
                BR[i+1][j][0] = BF[i+1][j][0]; BR[i+1][j][1] = BF[i+1][j][1];
            }
        
		}
	  } // for ij	

	
	  
	/*cout << "BR" << make_string3(BR) << endl;
	cout << "LR" << make_string3(LR) << endl;
	*/
	/*Dump Cells*/
	if (dump_cells)
	for (i=0; i < P.size()-1; i++)
	for (j=0; j < Q.size() -1; j++)
	{
		// octave indices!
            printf("cell(%d,%d):\n",i+1,j+1); // octave indices
            printf("\tLF(i,j)=[%f,%f] ",    LF[i][j][0],LF[i][j][1]);
            printf("\tBF(i,j)=[%f,%f] \n",  BF[i][j][0],BF[i][j][1]);
            printf("\tLR(i,j)=[%f,%f] ",    LR[i][j][0],LR[i][j][1]);
            printf("\tBR(i,j)=[%f,%f] \n",  BR[i][j][0],BR[i][j][1]);
            printf("\tLR(i,j+1)=[%f,%f] ",  LR[i][j+1][0],LR[i][j+1][1]);
            printf("\tBR(i+1,j)=[%f,%f] \n",BR[i+1][j][0],BR[i+1][j][1]);
		
	}
	
	
	bool decision = (BR[P.size()-1][Q.size()-2][1]==1)||(LR[P.size()-2][Q.size()-1][1]==1);
		return decision;

	}
	
	void bootstrap(TrajectoryType &P, TrajectoryType &Q)
	{
		matrix_resize(A,P.size(),Q.size());
		matrix_resize(B,P.size(),Q.size());
		matrix_resize(C,P.size(),Q.size());
		matrix_resize(D,P.size(),Q.size());

		matrix_resize(LF,P.size(),Q.size(),2);
		matrix_resize(BF,P.size(),Q.size(),2);
		matrix_resize(LR,P.size(),Q.size(),2);
		matrix_resize(BR,P.size(),Q.size(),2);
	
		
			
			for (size_t i=1; i < P.size(); i++) lP.push_back((*d2)(P[i-1],P[i]));
			for (size_t j=1; j < Q.size(); j++) lQ.push_back((*d2)(Q[j-1],Q[j]));
			matrix_resize(lPQ, P.size(),Q.size());
			max_lPQ = 0; min_lPQ = std::numeric_limits<double>::infinity();
			for (size_t i=0; i < P.size(); i++)
				for (size_t j=0; j < Q.size(); j++) 
				{
					lPQ[i][j] = (*d2)(P[i],Q[j]);
					if (lPQ[i][j] > max_lPQ)
						max_lPQ = lPQ[i][j];
					if (lPQ[i][j] < min_lPQ)
						min_lPQ = lPQ[i][j];

				}
		max_lPQ =sqrt(max_lPQ);
		min_lPQ = sqrt(min_lPQ);
		//std::cout << "lP " << tools::make_string(lP) << std::endl;
		//std::cout << "lQ " << tools::make_string(lQ) << std::endl;
		//std::cout << "lPQ " << make_string2(lPQ) << std::endl;
		matrix_resize(bP,P.size(),Q.size());
		matrix_resize(bQ,P.size(),Q.size());
		size_t I=P.size() -1;
		size_t J=Q.size() -1;
		
		double xp,xp1,yp,yp1,xq,xq1,yq,yq1;
		for (size_t i=0; i < P.size() -1; i++)
		{
			xp  = P[i][0];   
			yp = P[i][1];
			xp1 = P[i+1][0]; 
			yp1 = P[i+1][1];
			for (size_t j=0; j < Q.size() -1; j++)
			{
				xq = Q[j][0]; yq = Q[j][1];
				xq1 = Q[j+1][0]; yq1 = Q[j+1][1];
				bP[i][j] = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
				bQ[i][j] = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
			}
			
		}
		xp = P[I][0]; yp=P[I][1];
		for (size_t j=0; j < Q.size()-1; j++)
		{
			xq = Q[j][0]; yq = Q[j][1];
			xq1 = Q[j+1][0]; yq1 = Q[j+1][1];
			//bP[i][j] = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
			bQ[I][j] = 2*((xq-xp)*(xq1-xq)+(yq-yp)*(yq1-yq));
		}
		xq = Q[J][0]; yq=Q[J][1];
		for (size_t i=0; i < P.size()-1; i++)
		{
			xp  = P[i][0];   yp = P[i][1];
			xp1 = P[i+1][0]; yp1 = P[i+1][1];
			bP[i][J] = 2*((xp-xq)*(xp1-xp)+(yp-yq)*(yp1-yp));
		}
		xp = P[I-1][0]; yp = P[I-1][1];
		xp1 = P[I][0]; yp1= P[I][1];
		
		xq = Q[J-1][0]; yq = Q[J-1][1];
		xq1 = Q[J][0]; yq1= Q[J][1];
		bP[I][J] = 2*((xp-xq1)*(xp1-xp)+(yp-yq1)*(yp1-yp));
		bQ[I][J] = 2*((xq-xp1)*(xq1-xq)+(yq-yp1)*(yq1-yq));
			
	}
	
	
	 // Frechet compute by interval checking
	 DistanceType frechet_compute(TrajectoryType &P, TrajectoryType &Q)
	 {
		 // Calculate low and high
		 
		// Type A critical Values
		// all length between two points
	/* Just hack away all criticals
		std::vector<DistanceType> E;
		for (size_t i=0; i < P.size(); i++)
		for (size_t j=0; j < Q.size(); j++)
			E.push_back(sqrt(lPQ[i][j]));
		// Type B critical values
		for (size_t i=0; i < P.size()-1; i++)
		{
			double ap = lP[i];
			for (size_t j=0; j<Q.size() -1; j++)
			{
				double aq = lQ[j];
				double bp = bP[i][j];
				double tst = -bp / (2*ap);
				if (tst >= 0 && tst <= 1)
				  E.push_back(sqrt(-(((.5*bp)*(.5*bp))/ap)+lPQ[i][j]));
				double bq = bQ[i][j];
				tst = -bq / (2*aq);
				if (tst >= 0 && tst <= 1)
				  E.push_back(sqrt(-(((.5*bq)*(.5*bq))/aq)+lPQ[i][j]));
			}
		}
		
		
		std::sort(E.begin(),E.end());
		// now linear search, binary would be better
		double low = 0, high = max_lPQ;
		for (size_t i=0; i <E.size(); i++)
		{
			if (!decide(P,Q,E[i]))
			{
				low = E[i];
			}else{
				high = E[i];
				break;
			}
		}*/
		double low = min_lPQ, high = max_lPQ;
		 // Refinement
		//cout << low << "<= d <=" <<high << endl;
		while(fabs(low-high) > 0.01)
		{
			double epsilon = (low + high)/2;
			if (!decide(P,Q,epsilon))
			{
				low = epsilon;
			}else{
				high = epsilon;
			}			
		}
		
		//cout << low << "<= d <=" <<high << endl;
		return high;
		 
		 
		 
	 }
	 
	 
	// Ciritical Value Analysis follows
	DistanceType frechet_compute_critval2(TrajectoryType &P, TrajectoryType &Q)
	{
		bootstrap(P,Q);
		std::vector<DistanceType> E;
		// There are O(p2q + pq2) such critical values 
		
		//namely the distance between starting points and endpoints of P and Q 
		//(case a), , 
		E.push_back(lPQ[0][0]);
		E.push_back(lPQ[P.size()-1][Q.size()-1]);
		// the distances between vertices of one curve and edges of the other (case b)
		default_segment_distance<typename TrajectoryType::value_type> dseg;
		for (size_t i=0; i < P.size(); i++)
			for (size_t j=0; j < Q.size() -1; j++)
			    E.push_back(dseg(P[i],Q[j],Q[j+1]));
		for (size_t j=0; j < Q.size(); j++)
			for (size_t i=0; i < P.size() -1; i++)
			    E.push_back(dseg(Q[j],P[i],P[i+1]));
		//the common distance of two vertices of one curve to the intersection point 
		// of their bisector with some edge of the other (case c). 
		
		double xp,yp,xq,yq,xp1,yp1,xq1,yq1,yu,xu;
		
		for (size_t i=0; i< P.size() -1;i++)
		{
			xp  = P[i][0]; yp = P[i][1];
			xp1 = P[i+1][0]; yp1= P[i+1][1];
			for (size_t j=0; j<Q.size() -1; j++)
			{
				xq = Q[j][0]; yq=Q[j][1];
				xq1 = Q[j+1][0]; yq1=Q[j+1][1];
				for (size_t k=j+1; k < Q.size(); k++)
				{
					double den = bP[i][j]-bP[i][k];
					if (fabs(den) < 0.00001)
					{
						double u1 = (lPQ[i][k]-lPQ[i][j]) / den;
						if (u1 >= 0 && u1 <= 1)
						{
							double xu = xp + u1*(xp1-xp); double yu = yp + u1*(yp1-yp);
							double tst = sqrt((xu-xq)*(xu-xq)+ (yu-yq)*(yu-yq));
							E.push_back(tst);
						}
					}
				}
			}
			
		}
		for (size_t j=0; j < Q.size() -1; j++)
		{
			xq = Q[j][0]; yq=Q[j][1];
			xq1 = Q[j+1][0]; yq1=Q[j+1][1];
			for (size_t i=0; i < P.size() -1; i++)
			{
				xp  = P[i][0]; yp = P[i][1];
				xp1 = P[i+1][0]; yp1= P[i+1][1];
				for (size_t k=i+1; k < P.size(); k++)
				{
					double den = bQ[i][j]-bQ[k][j];
					if (fabs(den) < 0.00001)
					{
						double u2 = (lPQ[k][j]-lPQ[i][j]) / den;
						if (u2 >= 0 && u2 <= 1)
						{
							double xu = xq + u2*(xq1-xq); double yq = yq + u2*(yq1-yq);
							double tst = sqrt((xu-xp)*(xu-xp) + (yu-yp)*(yu-yp));
							 E.push_back(tst);
						}
					}
					
				}
				
			}
		}
		
		
		
		
		
		
		
		sort(E.begin(),E.end());
		//std::cout << "Critical Set" << tools::make_string(E) << endl;
		
		double low = 0, high = std::numeric_limits<double>::infinity();
		for (size_t i=0; i <E.size(); i++)
		{
			if (i>0 && E[i] == E[i-1])
				continue;
			if (!decide(P,Q,E[i]))
			{
				low = E[i];
			}else{
				high = E[i];
				break;
			}
		}
		std::cout << low << "<= d <=" <<high << std::endl;
		/*for (double a = low; a < high+(high-low)/2; a+=(high-low)/10)
			cout << "Decide @" << a << ":\t" << decide(P,Q,a) << endl;*/
		return high;
		
		
	}
	DistanceType frechet_compute_octave(TrajectoryType &P, TrajectoryType &Q)
	{
		bootstrap(P,Q);
		// Type A critical Values
		// all length between two points
		std::vector<DistanceType> E;
		for (size_t i=0; i < P.size(); i++)
		for (size_t j=0; j < Q.size(); j++)
			E.push_back(sqrt(lPQ[i][j]));
		// Type B critical values
		for (size_t i=0; i < P.size()-1; i++)
		{
			double ap = lP[i];
			for (size_t j=0; j<Q.size() -1; j++)
			{
				double aq = lQ[j];
				double bp = bP[i][j];
				double tst = -bp / (2*ap);
				if (tst >= 0 && tst <= 1)
				  E.push_back(sqrt(-(((.5*bp)*(.5*bp))/ap)+lPQ[i][j]));
				double bq = bQ[i][j];
				tst = -bq / (2*aq);
				if (tst >= 0 && tst <= 1)
				  E.push_back(sqrt(-(((.5*bq)*(.5*bq))/aq)+lPQ[i][j]));
				
			}
			
		}
		
		
		//cout << "Before sort" << endl;
		std::sort(E.begin(),E.end());
		//cout << "After sort" << endl;
		//std::cout << "Critical Set" << tools::make_string(E) << endl;
		// now linear search, binary would be better
		double low = 0, high = std::numeric_limits<double>::infinity();
		for (size_t i=0; i <E.size(); i++)
		{
			if (!decide(P,Q,E[i]))
			{
				low = E[i];
			}else{
				high = E[i];
				break;
			}
		}
		//cout << low << "<= d <=" <<high << endl;
		
		// Second pass ==> monotone increasing path from start to end
		double xp,yp,xq,yq,xp1,yp1,xq1,yq1,yu,xu;
		E.clear();
		for (size_t i=0; i< P.size() -1;i++)
		{
			xp  = P[i][0]; yp = P[i][1];
			xp1 = P[i+1][0]; yp1= P[i+1][1];
			for (size_t j=0; j<Q.size() -1; j++)
			{
				xq = Q[j][0]; yq=Q[j][1];
				xq1 = Q[j+1][0]; yq1=Q[j+1][1];
				for (size_t k=j+1; k < Q.size(); k++)
				{
					double den = bP[i][j]-bP[i][k];
					if (fabs(den) < 0.00001)
					{
						double u1 = (lPQ[i][k]-lPQ[i][j]) / den;
						if (u1 >= 0 && u1 <= 1)
						{
							double xu = xp + u1*(xp1-xp); double yu = yp + u1*(yp1-yp);
							double tst = sqrt((xu-xq)*(xu-xq)+ (yu-yq)*(yu-yq));
							if (tst > low && tst < high)
							   E.push_back(tst);
						}
					}
				}
			}
			
		}
		for (size_t j=0; j < Q.size() -1; j++)
		{
			xq = Q[j][0]; yq=Q[j][1];
			xq1 = Q[j+1][0]; yq1=Q[j+1][1];
			for (size_t i=0; i < P.size() -1; i++)
			{
				xp  = P[i][0]; yp = P[i][1];
				xp1 = P[i+1][0]; yp1= P[i+1][1];
				for (size_t k=i+1; k < P.size(); k++)
				{
					double den = bQ[i][j]-bQ[k][j];
					if (fabs(den) < 0.00001)
					{
						double u2 = (lPQ[k][j]-lPQ[i][j]) / den;
						if (u2 >= 0 && u2 <= 1)
						{
							double xu = xq + u2*(xq1-xq); double yq = yq + u2*(yq1-yq);
							double tst = sqrt((xu-xp)*(xu-xp) + (yu-yp)*(yu-yp));
							if (tst > low && tst < high)
							   E.push_back(tst);
						}
					}
					
				}
				
			}
		}
	    std::sort(E.begin(),E.end());
	    //cout << "Pass 2 critical values: "<<tools::make_string(E) << endl;
	    
	    for (size_t i=0; i < E.size(); i++)
	    {
			// the shortest one is the distance
			if (decide(P,Q,E[i]))
			{
			//	cout << "Second Pass did it" << E[i] << endl;
				return E[i];
			}
		}
		return high;
	}
	
	
	
	
	
    DistanceType operator()(TrajectoryType &u, TrajectoryType &v, double eps,
    element_distance<typename TrajectoryType::value_type> *_d2)
  {	
	  // Which squared distance?
	  this->d2 = _d2;
	  // Precalculate length
	  bootstrap(u,v);	
	  // either decide
	  if (eps >= 0)
	    return decide(u,v,eps);
	  // or scan (incomplete) 
	 return frechet_compute(u,v);
	  
  }
   
};      

template<class TrajectoryType>
double frechet(TrajectoryType &u, TrajectoryType &v, double eps = -1)
{
	frechet_distance_impl<TrajectoryType, double> df;
	default_element_distance_squared<typename TrajectoryType::value_type> d2;
	return df(u,v,eps,&d2);
}	

template<class TrajectoryType>
double frechet(TrajectoryType &u, TrajectoryType &v, double eps, double &maxDist)
{
	frechet_distance_impl<TrajectoryType, double> df;
	default_element_distance_squared<typename TrajectoryType::value_type> d2;
	double ret = df(u,v,eps,&d2);
	maxDist = df.getMaxPointDistance();
	return ret;
}	


}
#endif
