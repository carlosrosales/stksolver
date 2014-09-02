/******************************************************************************
* File      : stokesPostProcess_tria6.c                                       *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates velocity and pressure fields at the required points, and stores  *
* them in files 'velocity.dat' and 'presure.dat'.                             *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

#include "stokesPostProcess_tria6.h"

int stokesPostProcess_tria6(unsigned int ANALYSIS, double **Xinner, 
                            double **mNodes, unsigned int **mElems, 
                            double *vProbParam, unsigned int *vBCType, 
                            double *vB)
{
	FILE    *fv,*fp;
	unsigned int    i;
	double  P;
    double  U[3], Xin[3];

	switch(ANALYSIS){
		case 0:
			fprintf(file_log,"\n\tstokesPostProcess_tria6():");
			fprintf(file_log," doing velocity calculation ...");
			fv = fopen("velocity.dat","w");
			fprintf(fv,"x	y	z	Ux	Uy	Uz\n");
			break;
		case 1:
			fprintf(file_log,"\n\tstokesPostProcess_tria6():");
			fprintf(file_log," doing pressure calculation ...");
			fp = fopen("pressure.dat","w");
			fprintf(fp,"x	y	z	P\n");
			break;
		default:
			fprintf(file_log,"\n\tstokesPostProcess_tria6():");
			fprintf(file_log," doing velocity and pressure calculations ...");
			fv = fopen("velocity.dat","w");
			fp = fopen("pressure.dat","w");
			fprintf(fv,"x	y	z	Ux	Uy	Uz\n");
			fprintf(fp,"x	y	z	P\n");
			break;
	}

	/* Calculate velocity and pressure at the required points */
	for(i = 0; i < nInternalPoints; i++){
	   	Xin[0] = Xinner[i][0];
   		Xin[1] = Xinner[i][1];
   		Xin[2] = Xinner[i][2];
   		
        switch(ANALYSIS){
            case 0:	
                velocity_tria6(Xin,mNodes,mElems,vProbParam,vB,U);
                fprintf(fv,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fv,"%le\t%le\t%le\n",U[0],U[1],U[2]);
                break;
            case 1:
                P = pressure_tria6(Xin,mNodes,mElems,vB);
                fprintf(fp,"%le\t%le\t%le\t%le\n",Xin[0],Xin[1],Xin[2],P);
                break;
            default:
                velocity_tria6(Xin,mNodes,mElems,vProbParam,vB,U);
                P = pressure_tria6(Xin,mNodes,mElems,vB);
                fprintf(fv,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fv,"%le\t%le\t%le\n",U[0],U[1],U[2]);
                fprintf(fp,"%le\t%le\t%le\t%le\n",Xin[0],Xin[1],Xin[2],P);
        }
    }

    switch(ANALYSIS){
        case 0:	fclose(fv); break;
        case 1:	fclose(fp); break;
        default: fclose(fv); fclose(fp);
	}

	return 0;
}
