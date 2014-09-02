/******************************************************************************
* File      : intSingularHStk_tria3.c                                         *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function for the pressure from Stokes  *
* equation for one element when the collocation point is one of the nodes of  *
* the integration element, and returns the values for all shape functions in  *
* matrix Int[][].                                                             *
* Uses element subdivision and standard Gauss-Legendre quadrature integration.*
* The nodes are rotated until the singular point is in node 3 to simplify the *
* subdivision algorithm (in this case only requires dividing by 2 in each     *
* subsequent subdivision).                                                    *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
*                                                                             *
* I[1][2] = contribution to P from Phi_2 and node 3                           *
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

int intSingularHStk_tria3(unsigned int SinNode, double X[3][3], double *Xeval, 
                          double Int[3][3])
{
    const unsigned int  NODES_IN_ELEM = 3;
    unsigned int    i, j, k, m;
    double   A, dx, dy, dz, J, JL, dr, W, M1, M2, M3;
    double   B[3], L[2], XL[3], N[3], xs[3][3][2];
 
	/* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        Int[0][i] = 0.0;
        Int[1][i] = 0.0;
        Int[2][i] = 0.0;
    }
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            for(k = 0; k <2; k++) xs[i][j][k] = 0.0;
	
    /* Reorder so that singularity is at node 3 */
    if(SinNode == 1){
        for(i = 0; i < 3; i++) shift(X[2][i],X[0][i],X[1][i]);
    }
    else if(SinNode == 2){
        for(i = 0; i < 3; i++) shift(X[2][i],X[1][i],X[0][i]);
    }

    /* Setup the 3 initial triangles */
    xs[0][0][0] = 1.0;	xs[0][0][1] = 0.0;
    xs[0][1][0] = 0.5;	xs[0][1][1] = 0.5;
    xs[0][2][0] = 0.5;	xs[0][2][1] = 0.0;

    xs[1][0][0] = 0.5;	xs[1][0][1] = 0.0;
    xs[1][1][0] = 0.5;	xs[1][1][1] = 0.5;
    xs[1][2][0] = 0.0;	xs[1][2][1] = 0.5;
	
    xs[2][0][0] = 0.5;	xs[2][0][1] = 0.5;
    xs[2][1][0] = 0.0;	xs[2][1][1] = 1.0;
    xs[2][2][0] = 0.0;	xs[2][2][1] = 0.5;
	
    /* Subdivision loop */
    J = 1.0;
    for(m = 0; m < NSUBDIVISIONS; m++){
        J = J*0.25;
    
        /* Integrate over the 3 triangles furthest away from singular */
        /* point using standard Gauss-Legendre quadrature             */
        for(i = 0; i < 3; i++){
            for(j = 0; j < TSNGAUSS; j++){
                M1 = TSGauss[j][0];
                M2 = TSGauss[j][1];
                M3 = 1.0 - M1 - M2;
                W = TSGauss[j][2];
	    		
                /* Get L and N as functions of the variables */ 
                /* in the unit triangle (T1,T2)              */
                L[0] = xs[i][0][0]*M1 + xs[i][1][0]*M2 + xs[i][2][0]*M3;
                L[1] = xs[i][0][1]*M1 + xs[i][1][1]*M2 + xs[i][2][1]*M3;
                shape_tria3(L,N);
                JL = X2L_tria3(X,XL,L,N);

                /* Auxiliar quantities */
                dx = Xeval[0] - XL[0];
                dy = Xeval[1] - XL[1];
                dz = Xeval[2] - XL[2];
                dr = sqrt( dx*dx + dy*dy + dz*dz );
	    		A = J*JL*W/(dr*dr*dr);

                /* Get contribution from each node in the element */
                B[0] = A*dx;
                B[1] = A*dy;
                B[2] = A*dz;
                for(k = 0; k < NODES_IN_ELEM; k++){
                    Int[0][k] += B[0]*N[k];
                    Int[1][k] += B[1]*N[k];
                    Int[2][k] += B[2]*N[k];
                }
            }
        }
    
        /* Set triangles for new subdivision */
        for(i = 0; i < 3; i++) {
            for(j = 0; j < 3; j++){
			   	xs[i][j][0] = xs[i][j][0]*0.5;
			   	xs[i][j][1] = xs[i][j][1]*0.5;
            }
        }
    }

    /* Unscramble node contributions */
    if(SinNode == 1){
        for(i = 0; i < 3; i++) shift(Int[i][1],Int[i][0],Int[i][2]);
    }
    else if(SinNode == 2){
        for(i = 0; i < 3; i++) shift(Int[i][0],Int[i][1],Int[i][2]);
    }

    return 0;
}



















