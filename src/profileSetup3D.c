/******************************************************************************
* File      : profileSetup3D.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Provides a fully developed flow profile for nodes marked as inlet and       *
* outlet (vBCType == 1). Active only if teh given velocity is negative.       *
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

#include <stdio.h>
#include <stdlib.h>

extern unsigned int nNodes;

int profileSetup3D(double Ly, double Lz, double maxU, double **mNodes,
                   double **mBC)
{
    int i;
    double ny, nz, y, z;
    
    ny = 2.0/Ly;
    nz = 2.0/Lz;
    for(i = 0; i < nNodes; i++){
        if(mBC[i][0] < 0){
            y = mNodes[i][1]*ny;
            z = mNodes[i][2]*nz;
            mBC[i][0] = maxU*( 1.0 - y*y )*( 1.0 - z*z );
            mBC[i][1] = 0.0;
            mBC[i][2] = 0.0;
        }
    }
    
    return 0;
}

