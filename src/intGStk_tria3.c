/******************************************************************************
* Name      : intGStk_tria3.c                                                 *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function of the Stokes equation for    *
* one element and returns the integrals in Int[][][] for all shape functions. *
* Uses standard gaussian quadrature in a triangle with TNGAUSS points.        *
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

int intGStk_tria3(double X[3][3], double *Xeval, double Int[3][3][3])
{
    const unsigned int NODES_IN_ELEM = 3;
    unsigned int i, j, k;
    double A, B, J, dr, W;
    double C[3], dx[3], L[2], N[3], XL[3];

    /* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        for(j = 0; j < 3; j++){
            Int[0][j][i] = 0.0;
            Int[1][j][i] = 0.0;
            Int[2][j][i] = 0.0;
        }
    }

    /* Integrate using gaussian quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];

        /* Get shape functions N for these (L1,L2) values */
        shape_tria3(L,N);		
        J = X2L_tria3(X,XL,L,N);

        /* Get auxiliar quantities */
        dx[0] = Xeval[0] - XL[0];
        dx[1] = Xeval[1] - XL[1];
        dx[2] = Xeval[2] - XL[2];
        dr = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
        A = J*W/dr;
        B = A/(dr*dr);

	    /* Get contribution from each node in the element */
		C[0] = B*dx[0];
		C[1] = B*dx[1];
		C[2] = B*dx[2];
        for(j = 0; j < NODES_IN_ELEM; j++){
            for(k = 0; k < 3; k++){
                Int[0][k][j] += C[0]*dx[k]*N[j];
                Int[1][k][j] += C[1]*dx[k]*N[j];
                Int[2][k][j] += C[2]*dx[k]*N[j];
            }
            Int[0][0][j] += A*N[j];
            Int[1][1][j] += A*N[j];
            Int[2][2][j] += A*N[j];
        }
	}

	return 0;
}
