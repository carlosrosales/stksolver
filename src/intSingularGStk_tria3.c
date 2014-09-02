/******************************************************************************
* Name      : intSingularGStk_tria3.c                                         *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function G when a collocation point is *
* one of the nodes of the integration element and returns the values for all  *
* shape functions in Int[]. Uses a regularization transformation that changes *
* the triangle into a degenerated square and Gauss-Jacobi quadrature. When    *
* the evaluation node falls in a triangle edge the element is divided in two  *
* and the the transformation is applied to both subtriangles.                 *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
*                                                                             *
* Int[0][1][2]    : contribution to u_x term from Phi_2 and node 3.           *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of stkSolver.                                             *
*                                                                             *
* stkSolver is free software; you can redistribute it and/or modify           *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* stkSolver is distributed in the hope that it will be useful,                *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with stkSolver; if not, write to the Free Software                    *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include "integral_tria3.h"

int intSingularGStk_tria3(unsigned int SinNode, double X[3][3], double *Xeval,
                       double Int[3][3][3])
{
    const unsigned int NODES_IN_ELEM = 3;
	unsigned int    i, j, k, m;
	double   A, B, dr, J, U1, U2, W1, W2;
	double   C[3], dx[3], L[2], N[3], XL[3];
 
    /* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        for(j = 0; j < 3; j++){
            Int[0][j][i] = 0.0;
            Int[1][j][i] = 0.0;
            Int[2][j][i] = 0.0;
        }
    }

    /* Integrate using Gauss-Jacobi quadrature */
    for(i = 0; i < NGAUSS; i++){
        U1 = Gauss[i][0];
        W1 = Gauss[i][1];
    
        for(j = 0; j < NGAUSS; j++){
            U2 = GaussJacobi[j][0];
            W2 = GaussJacobi[j][1];

            /* Obtain L values for these (U1,U2) set */
            if(SinNode == 1){
                L[0] = 0.5*(1.0 - U2);
                L[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }
            else if(SinNode == 2){
                L[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L[1] = 0.5*(1.0 - U2);
            }
            else{
                L[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L[1] = 0.25*(1.0 - U1)*(1.0 + U2);
            }

        	/* Contribution from (i,j) */
			shape_tria3(L,N);
			J = X2L_tria3(X,XL,L,N);
            dx[0] = Xeval[0] - XL[0];
            dx[1] = Xeval[1] - XL[1];
            dx[2] = Xeval[2] - XL[2];
            dr = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
            A = W1*W2*J*0.125/dr;
			B = A/(dr*dr);
			
			/* Get contributions from each node in the element */
		    C[0] = B*dx[0];
		    C[1] = B*dx[1];
		    C[2] = B*dx[2];
            for(k = 0; k < NODES_IN_ELEM; k++){
                for(m = 0; m < 3; m++){
                    Int[0][m][k] += C[0]*dx[m]*N[k];
                    Int[1][m][k] += C[1]*dx[m]*N[k];
                    Int[2][m][k] += C[2]*dx[m]*N[k];
                    switch(m){
                        case 0: Int[0][0][k] += A*N[k]; break;
                        case 1: Int[1][1][k] += A*N[k]; break;
                        case 2: Int[2][2][k] += A*N[k];
                    }
                }
            }
    	}
	}

	return 0;
}
