/******************************************************************************
* Name      : intSingularGStk_tria6.c                                         *
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
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

#include "integral_tria6.h"

int intSingularGStk_tria6(unsigned int SinNode, double X[6][3], double *Xeval, 
                          double Int[3][3][6])
{
    const unsigned int  NODES_IN_ELEM = 6;
    unsigned int    i, j, k, m, T;
    double   A, B, C, dr, J, U1, U2, W1, W2;
    double   D[3], dx[3], L1[2], L2[2], N1[6], N2[6], XL[3];
 
    /* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        for(j = 0; j < 3; j++){
            Int[0][j][i] = 0.0;
            Int[1][j][i] = 0.0;
            Int[2][j][i] = 0.0;
        }
    }

	/* Choose number of subelements */
    if(SinNode%2 == 1){
        T = 1;
        A = 0.125;
    }
    else{
        T = 2;
        A = 0.0625;
    }

    /* Integrate using Gauss-Jacobi quadrature */
    for(i = 0; i < NGAUSS; i++){
        U1 = Gauss[i][0];
        W1 = Gauss[i][1];
    
        for(j = 0; j < NGAUSS; j++){
            U2 = GaussJacobi[j][0];
            W2 = GaussJacobi[j][1];

            /* Obtain L values for this (U1,U2) set */
            if(SinNode == 1){
                L1[0] = 0.5*(1.0 - U2);
                L1[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }
            else if(SinNode == 2){
                L1[0] = 0.25*(1.0 - U2);
                L1[1] = 0.5 + 0.25*U1*(1.0 + U2);
                L2[0] = 0.5 - 0.25*U1*(1.0 + U2);
                L2[1] = 0.25*(1.0 - U2);
            }
            else if(SinNode == 3){
                L1[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L1[1] = 0.5*(1.0 - U2);
            }
            else if(SinNode == 4){
                L1[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L1[1] = 0.5 - 0.25*U1*(1.0 + U2);
                L2[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L2[1] = 0.25*(1.0 - U2);
            }
            else if(SinNode == 5){
                L1[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L1[1] = 0.25*(1.0 - U1)*(1.0 + U2);
            }
            else{
                L1[0] = 0.5 + 0.25*U1*(1.0 + U2);
                L1[1] = 0.25*(1.0 - U1)*(1.0 + U2);
                L2[0] = 0.25*(1 - U2);
                L2[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }

            /* Contribution from triangle 1 */
            shape_tria6(L1,N1);
            J = X2L_tria6(X,XL,L1,N1);
            dx[0] = Xeval[0] - XL[0];
            dx[1] = Xeval[1] - XL[1];
            dx[2] = Xeval[2] - XL[2];
            dr = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
            B = W1*W2*J*A/dr;
			C = B/(dr*dr);
	
			/* Get contribution from each node in the element */
		    D[0] = C*dx[0];
		    D[1] = C*dx[1];
		    D[2] = C*dx[2];
            for(k = 0; k < NODES_IN_ELEM; k++){
                for(m = 0; m < 3; m++){
                    Int[0][m][k] += D[0]*dx[m]*N1[k];
                    Int[1][m][k] += D[1]*dx[m]*N1[k];
                    Int[2][m][k] += D[2]*dx[m]*N1[k];
                }
                Int[0][0][k] += B*N1[k];
                Int[1][1][k] += B*N1[k];
                Int[2][2][k] += B*N1[k];
            }

			/* Contribution from triangle 2 (if necessary) */
			if(T == 2){
                shape_tria6(L2,N2);
                J = X2L_tria6(X,XL,L2,N2);
                dx[0] = Xeval[0] - XL[0];
                dx[1] = Xeval[1] - XL[1];
                dx[2] = Xeval[2] - XL[2];
                dr = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
                B = W1*W2*J*A/dr;
			    C = B/(dr*dr);
	
			    /* Get contribution from each node in the element */
		        D[0] = C*dx[0];
		        D[1] = C*dx[1];
		        D[2] = C*dx[2];
                for(k = 0; k < NODES_IN_ELEM; k++){
                    for(m = 0; m < 3; m++){
                        Int[0][m][k] += D[0]*dx[m]*N2[k];
                        Int[1][m][k] += D[1]*dx[m]*N2[k];
                        Int[2][m][k] += D[2]*dx[m]*N2[k];
                    }
                    Int[0][0][k] += B*N2[k];
                    Int[1][1][k] += B*N2[k];
                    Int[2][2][k] += B*N2[k];
                }
			}
    	}
	}

	return 0;
}
