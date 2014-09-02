/******************************************************************************
* File      : stokesFormA_tria3.c                                             *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Forms the coefficient matrix mA by assembling submatrices of the element,   *
* SubA, for the IBEM Stokes problem with one rigid body inclusion.            *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
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

#include "constants.h"
#include "stokesFormA_tria3.h"

int  stokesFormA_tria3(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, double *vProbParam, double **mA, 
                       double **mBC, double *vB)
{
    const unsigned int NODES_IN_ELEM = 3, SIZE = 3*nNodes;
    int currentNode, i, j, k, SinNode, xCol, xRow, yCol, yRow, zCol, zRow;  
    double preFactor, viscosity;
    double Xeval[3];
    double X[3][3];
    double INT[3][3][3];

    /* Initialization */
    viscosity = vProbParam[0];
    preFactor = 8.0*pi*viscosity;
    for(i = 0; i < SIZE; i++)
        for(j = 0; j < SIZE; j++) mA[i][j] = 0.0;
    
    for(i = 0; i < nNodes; i++){
        Xeval[0] = mNodes[i][0];
        Xeval[1] = mNodes[i][1];
        Xeval[2] = mNodes[i][2];
        xRow = i;
        yRow = xRow + nNodes;
        zRow = yRow + nNodes;

        /* Check for singular case */
        for(j = 0; j < nElems; j++){
            SinNode = 0;
            for(k = 0; k < NODES_IN_ELEM; k++){
                currentNode = mElems[j][k] - 1;
                X[k][0] = mNodes[currentNode][0];
                X[k][1] = mNodes[currentNode][1];
                X[k][2] = mNodes[currentNode][2];
                if(currentNode == i) SinNode = k + 1;
            }

            if(SinNode == 0) intGStk_tria3(X,Xeval,INT);
            else intSingularGStk_tria3(SinNode,X,Xeval,INT);
    		
            /* Given velocity boundary condition */
            for(k = 0; k < NODES_IN_ELEM; k++){		
                xCol = mElems[j][k] - 1;
                yCol = xCol + nNodes;
                zCol = yCol + nNodes;
					
                mA[xRow][xCol] -= INT[0][0][k];
                mA[xRow][yCol] -= INT[0][1][k];
                mA[xRow][zCol] -= INT[0][2][k];
					
                mA[yRow][xCol] -= INT[1][0][k];
                mA[yRow][yCol] -= INT[1][1][k];
                mA[yRow][zCol] -= INT[1][2][k];
					
                mA[zRow][xCol] -= INT[2][0][k];
                mA[zRow][yCol] -= INT[2][1][k];
                mA[zRow][zCol] -= INT[2][2][k];
            }
    	}

    	/* Assemble right hand side vector */
    	vB[xRow] = mBC[i][0]*preFactor;
    	vB[yRow] = mBC[i][1]*preFactor;
    	vB[zRow] = mBC[i][2]*preFactor;
	}

	return 0;
}
