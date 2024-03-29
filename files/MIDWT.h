/*
File Name: MIDWT.c
Last Modification Date:	06/14/95	13:01:15
Current Version: MIDWT.c	2.4
File Creation Date: Wed Oct 12 08:44:43 1994
Author: Markus Lang  <lang@jazz.rice.edu>

Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
Created by Markus Lang, Department of ECE, Rice University. 

This software is distributed and licensed to you on a non-exclusive 
basis, free-of-charge. Redistribution and use in source and binary forms, 
with or without modification, are permitted provided that the following 
conditions are met:

1. Redistribution of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.
2. Redistribution in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software 
   must display the following acknowledgment: This product includes 
   software developed by Rice University, Houston, Texas and its contributors.
4. Neither the name of the University nor the names of its contributors 
   may be used to endorse or promote products derived from this software 
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS, 
AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY 
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For information on commercial licenses, contact Rice University's Office of 
Technology Transfer at techtran@rice.edu or (713) 348-6173

Change History: Fixed the code such that 1D vectors passed to it can be in
                either passed as a row or column vector. Also took care of 
		the code such that it will compile with both under standard
		C compilers as well as for ANSI C compilers
		Jan Erik Odegard <odegard@ece.rice.edu> Wed Jun 14 1995

		Fix minor bug to allow maximum number of levels

decription of the matlab call:
%y = midwt(x,h,L);
% 
% function computes the inverse discrete wavelet transform y for a 1D or 2D
% input signal x.
%
%    Input:
%	x    : finite length 1D or 2D input signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: mdwt, mrdwt, mirdwt

*/

#include <math.h>
#include <stdio.h>

#define mymax(A,B) (A > B ? A : B)
#define mymat(a, i, j) (*(a + (m*(j)+i)))  /* macro for matrix indices */

//void bpsconv(float *x_out, int lx, float *g0, float *g1, int lhm1, 
//	int lhhm1, float *x_inl, float *x_inh);
//void MIDWT(float *x, int m, int n, float *h, int lh, int L, float *y);


class MIDWT_cls{
public:
	static void MIDWT(float *y, int m, int n, float *h, int lh, int L, float *x){
		float  *g0, *g1, *ydummyl, *ydummyh, *xdummy;
		long i;

		int actual_L, actual_m, actual_n, r_o_a, c_o_a, ir, ic, lhm1, lhhm1, sample_f;
	
	/*	xdummy = new float[max(m,n)];
		ydummyl = new float[max(m,n)+lh/2-1];
		ydummyh = new float[max(m,n)+lh/2-1];
		g0 = new float[lh];
		g1 = new float[lh];
		cout << "MAX " << max(m,n) << endl;
		//*/

		xdummy = (float *)malloc(mymax(m,n)*sizeof(float));
		ydummyl = (float *)malloc((mymax(m,n)+lh/2-1)*sizeof(float));
		ydummyh = (float *)malloc((mymax(m,n)+lh/2-1)*sizeof(float));
		g0 = (float *)malloc(lh*sizeof(float));
		g1 = (float *)malloc(lh*sizeof(float));
	
	
		if (n==1){
			n = m;
			m = 1;
		}
		/* synthesis lowpass and highpass */
		for (i=0; i<lh; i++){
			g0[i] = h[i];
			g1[i] = h[lh-i-1];
		}
		for (i=1; i<=lh; i+=2)
			g1[i] = -g1[i];
  
		lhm1 = lh - 1;
		lhhm1 = lh/2 - 1;
		/* 2^L */
		sample_f = 1;
		for (i=1; i<L; i++)
			sample_f = sample_f*2;
  
		if (m>1)
			actual_m = m/sample_f;
		else 
			actual_m = 1;
		actual_n = n/sample_f;

		for (i=0; i<(m*n); i++)
			x[i] = y[i];
  
		/* main loop */
		for (actual_L=L; actual_L >= 1; actual_L--){
			r_o_a = actual_m/2;
			c_o_a = actual_n/2;
    
			/* go by columns in case of a 2D signal*/
			if (m>1){
				for (ic=0; ic<actual_n; ic++){            /* loop over column */
					/* store in dummy variables */
					ir = r_o_a;
					for (i=0; i<r_o_a; i++){    
						ydummyl[i+lhhm1] = mymat(x, i, ic);
						ydummyh[i+lhhm1] = mymat(x, ir++, ic);
					}
					/* perform filtering lowpass and highpass*/
					bpsconv(xdummy, r_o_a, g0, g1, lhm1, lhhm1, ydummyl, ydummyh); 
					/* restore dummy variables in matrix */
					for (i=0; i<actual_m; i++)
						mymat(x, i, ic) = xdummy[i];
				}
			}
			/* go by rows */
			for (ir=0; ir<actual_m; ir++){            /* loop over rows */
				/* store in dummy variable */
				ic = c_o_a;
				for  (i=0; i<c_o_a; i++){    
					ydummyl[i+lhhm1] = mymat(x, ir, i);
					ydummyh[i+lhhm1] = mymat(x, ir, ic++);
				} 
				/* perform filtering lowpass and highpass*/
				bpsconv(xdummy, c_o_a, g0, g1, lhm1, lhhm1, ydummyl, ydummyh); 
				/* restore dummy variables in matrices */
				for (i=0; i<actual_n; i++)
					mymat(x, ir, i) = xdummy[i];
			}  
			if (m==1)
				actual_m = 1;
			else
				actual_m = actual_m*2;
			actual_n = actual_n*2;
		}

		free(xdummy);free(ydummyh);free(ydummyl);free(g0);free(g1);
	}

	static void bpsconv(float *x_out, int lx, float *g0, float *g1, int lhm1,
		int lhhm1, float *x_inl, float *x_inh){
		int i, j, ind, tj;
		float x0, x1;

		for (i=lhhm1-1; i > -1; i--){
			x_inl[i] = x_inl[lx+i];
			x_inh[i] = x_inh[lx+i];
		}
		ind = 0;

		for (i=0; i<(lx); i++){
			x0 = 0;
			x1 = 0;
			tj = -2;
			for (j=0; j<=lhhm1; j++){
				tj+=2;
				x0 = x0 + x_inl[i+j]*g0[lhm1-1-tj] + x_inh[i+j]*g1[lhm1-1-tj] ;
				x1 = x1 + x_inl[i+j]*g0[lhm1-tj] + x_inh[i+j]*g1[lhm1-tj] ;
			}
			x_out[ind++] = x0;
			x_out[ind++] = x1;
		}
	}

};