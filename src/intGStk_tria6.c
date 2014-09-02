/******************************************************************************
* Name      : intGStk_tria6.c                                                 *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function of the Stokes equation for    *
* one element and returns the integrals in Int[][][] for all shape functions: *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
* Uses standard gaussian quadrature in a triangle with TNGAUSS points.        *
*                                                                             *
* Int[0][1][2]    : contribution to u_x term from Phi_y at node 3.            *
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

int intGStk_tria6(double X[6][3], double *Xeval, double Int[3][3][6])
{
    const unsigned int  NODES_IN_ELEM = 6;
    unsigned int    i, j, k;
    double  A, B, Cx, Cy, Cz, J, dr, W;
    double  dx[3], L[2], N[6], XL[3];

    /* Initialize */
    for(i = 0; i < 3; i++){
        for(j = 0; j < NODES_IN_ELEM; j++){
            Int[0][i][j] = 0.0;
            Int[1][i][j] = 0.0;
            Int[2][i][j] = 0.0;
        }
    }

    /* Integrate using gaussian quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];

        /* Get shape functions N for these (L1,L2) values */
        shape_tria6(L,N);		
        J = X2L_tria6(X,XL,L,N);

        /* Get auxiliar quantities */
        dx[0] = Xeval[0] - XL[0];
        dx[1] = Xeval[1] - XL[1];
        dx[2] = Xeval[2] - XL[2];
        dr = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
        A = J*W/dr;
        B = A/(dr*dr);

	    /* Get contribution from each node in the element */
		Cx = B*dx[0];
		Cy = B*dx[1];
		Cz = B*dx[2];
        for(k = 0; k < NODES_IN_ELEM; k++){
            for(j = 0; j < 3; j++){
                Int[0][j][k] += Cx*dx[j]*N[k];
                Int[1][j][k] += Cy*dx[j]*N[k];
                Int[2][j][k] += Cz*dx[j]*N[k];
            }
            Int[0][0][k] += A*N[k];
            Int[1][1][k] += A*N[k];
            Int[2][2][k] += A*N[k];
        }
	}

	return 0;
}
