/******************************************************************************
* File      : intSingularHStk_tria6.c                                         *
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
* If the singular node lies on an edge the original triangle is divided into  *
* two triangles both with the singular point at node 5. Otherwise the nodes   *
* are rotated until the singular point is in node 5 to simplify the           *
* subdivision algorithm (in this case only requires dividing by 2 in each     *
* subsequent subdivision).                                                    *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

#include "integral_tria6.h"

int intSingularHStk_tria6(unsigned int SinNode, double X[6][3], double *Xeval, 
                          double Int[3][6])
{
    const unsigned int  NODES_IN_ELEM = 6;
    unsigned int    i, j, k, m;
    double   A, A1, A2, dx, dx1, dx2, dy, dy1, dy2, dz, dz1, dz2, dr, dr1, dr2;
    double   J, JL, JL1, JL2, M1, M2, M3, U1, U2, U3, W;
    double   B[3], B1[3], B2[3], L[2], L1[2], L2[2], N[6], N1[6], N2[6];
    double   XL[3], XL1[3], XL2[3],  xs[3][3][2];
 
	/* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        Int[0][i] = 0.0;
        Int[1][i] = 0.0;
        Int[2][i] = 0.0;
    }
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            for(k = 0; k <2; k++) xs[i][j][k] = 0.0;
	
    /* Reorder so that singularity is at node 5 */
    if(SinNode == 1){
        for(i = 0; i < 3; i++){
            shift(X[4][i],X[0][i],X[2][i]);
            shift(X[5][i],X[1][i],X[3][i]);	    	
        }
    }
    else if(SinNode == 3){
        for(i = 0; i < 3; i++){
            shift(X[4][i],X[2][i],X[0][i]);
            shift(X[3][i],X[1][i],X[5][i]);
        }
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
                if(SinNode%2 == 1){        /* Node at corner */
                    L[0] = xs[i][0][0]*M1 + xs[i][1][0]*M2 + xs[i][2][0]*M3;
                    L[1] = xs[i][0][1]*M1 + xs[i][1][1]*M2 + xs[i][2][1]*M3;
                    shape_tria6(L,N);
                    JL = X2L_tria6(X,XL,L,N);
	    			
                    dx = Xeval[0] - XL[0];
                    dy = Xeval[1] - XL[1];
                    dz = Xeval[2] - XL[2];
                    dr = sqrt( dx*dx + dy*dy + dz*dz );
                    A = J*JL*W/(dr*dr*dr);
					B[0] = A*dx;
                    B[1] = A*dy;
                    B[2] = A*dz;
                    for(k = 0; k < NODES_IN_ELEM; k++){
                        Int[0][k] += B[0]*N[k];
                        Int[1][k] += B[1]*N[k];
                        Int[2][k] += B[2]*N[k];
                    }
                }
                else{
                    if(SinNode == 2){
                        U1 = 0.5*M3; U2 = M1 + 0.5*M3; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1 + 0.5*M3; U2 = M3; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
                    else if(SinNode == 4){
                        U1 = M1; U2 = 0.5*M3; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1; U2 = M2 + 0.5*M3; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
	    			else{
                        U1 = 0.5*M3; U2 = M1; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1 + 0.5*M3; U2 = M2; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
                    shape_tria6(L1,N1);
                    shape_tria6(L2,N2);
                    JL1 = X2L_tria6(X,XL1,L1,N1);
                    JL2 = X2L_tria6(X,XL2,L2,N2);	

                    dx1 = Xeval[0] - XL1[0];
                    dy1 = Xeval[1] - XL1[1];
                    dz1 = Xeval[2] - XL1[2];
                    dx2 = Xeval[0] - XL2[0];
                    dy2 = Xeval[1] - XL2[1];
                    dz2 = Xeval[2] - XL2[2];
                    dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
                    dr2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                    A1 = 0.5*J*JL1*W/(dr1*dr1*dr1);
                    A2 = 0.5*J*JL2*W/(dr2*dr2*dr2);
	    			
                    B1[0] = A1*dx1;
                    B1[1] = A1*dy1;
                    B1[2] = A1*dz1;
                    B2[0] = A2*dx2;
                    B2[1] = A2*dy2;
                    B2[2] = A2*dz2;
                    for(k = 0; k < NODES_IN_ELEM; k++){
                        Int[0][k] += B1[0]*N1[k] + B2[0]*N2[k];
                        Int[1][k] += B1[1]*N1[k] + B2[1]*N2[k];
                        Int[2][k] += B1[2]*N1[k] + B2[2]*N2[k];
                    }
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
        for(i = 0; i < 3; i++){
            shift(Int[i][2],Int[i][0],Int[i][4]);
            shift(Int[i][3],Int[i][1],Int[i][5]);
        }
    }
    else if(SinNode == 3){
        for(i = 0; i < 3; i++){
            shift(Int[i][0],Int[i][2],Int[i][4]);
            shift(Int[i][5],Int[i][1],Int[i][3]);
        }
    }

    return 0;
}



















