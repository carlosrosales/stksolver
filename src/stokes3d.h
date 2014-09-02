/******************************************************************************
* File      : stokes3d.h                                                      *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

/* Common function prototypes */
void errorHandler(char errorText[]);

void comFilter(char *FileName);

void freeDoubleMatrix(double **M, unsigned int ROWS);

void freeDoubleTensor(double ***T, unsigned int ROWS, unsigned int COLS);

void freeUintMatrix(unsigned int **M, unsigned int ROWS);

void solverGMRES(unsigned int preCond, unsigned int nInit, unsigned int size,
                 double **Amatrix, double *Rhs);

unsigned int *uintVector(unsigned int LENGTH, unsigned int INIT);

unsigned int **uintMatrix(unsigned int ROWS, unsigned int COLS, 
                            unsigned int INIT);

int elemType();
                       
int gaussBksb(unsigned int N, double **A, double *B);

int profileSetup3D(double Ly, double Lz, double maxU, double **mNodes,
                   double **mBC);
                   
int update(double **mNodes, double **mBC, double *Xp, double *Vp, double *Wp,
           double *Lp, double mass, double Rp, double *Fp, double *Tp);

int  stokesFormA_tria3(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, double *vProbParam, double **mA, 
                       double **mBC, double *vB);
                       
int  stokesFormA_tria6(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, double *vProbParam, double **mA, 
                       double **mBC, double *vB);
                       
int stokesPostProcess_tria3(unsigned int ANALYSIS, double **Xinner, 
                            double **mNodes, unsigned int **mElems, 
                            double *vProbParam, unsigned int *vBCType, 
                            double *vB);
                            
int stokesPostProcess_tria6(unsigned int ANALYSIS, double **Xinner, 
                            double **mNodes, unsigned int **mElems, 
                            double *vProbParam, unsigned int *vBCType, 
                            double *vB);

double *doubleVector(unsigned int LENGTH, unsigned int INIT);

double **doubleMatrix(unsigned int ROWS, unsigned int COLS, 
                       unsigned int INIT);

double ***doubleTensor(unsigned int ROWS, unsigned int COLS, 
                        unsigned int DEPTH, unsigned int INIT);

/* Define global variables */
FILE *file_log;
char cElemType[6];
double dt;
unsigned int nDim, nElems, nElemsWall, nElemType, nInternalPoints, nInterfaces,
             nMats, nNodes, nNodesWall, nNodesInElem, steps;
