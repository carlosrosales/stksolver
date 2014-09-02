/******************************************************************************
* File      : stokes3d.c                                                      *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Main function for Stokes 3D solver with triangular elements. Uses indirect  *
* boundary element method to solve for the case of given velocities only.     *
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
#include "gaussData.h"
#include "stokes3d.h"

int main()
{
    FILE *fIn, *fAux;
    time_t currentTime;
    struct tm *localTime;
    char cFilename[]="input.bem", cBuffer[32], cBCFile[32], cElemFile[32], 
         cNodeFile[32], cInternalPointsFile[32], cSection[20], cTmp[32], 
         solver[10];
    double exchanges, t1, t2, t3, t4, t_tot, t1s, t2s, t3s, t4s, t_tots, Rp, 
           mass, maxU, runtime, Ly, Lz;
    double Fp[3], Lp[3], Tp[3], Vp[3], Wp[3], Xp[3];
    double *vProbParam, *vB;
    double **mA, **mBC, **mNodes, **Xinner;
    int t1m, t2m, t3m, t4m, t_totm;
    unsigned int ANALYSIS, count, i, j, nBodies, nBuffer, nInit, preCond, size;
    unsigned int *vBCType, *indx;
    unsigned int **mElems;

    /* Open and initialize log file */
    file_log = fopen("bem.log","a");
    fprintf(file_log,"stokes3d: starting @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);

    /* Clean and open main input file */
    comFilter(cFilename);
    fIn=fopen(cFilename,"r");

    /* Read node section */
    fscanf(fIn,"%s %d %d %s",cSection,&nNodesWall,&nNodes,cNodeFile);

    /* Read elements section */
    fscanf(fIn,"%s %d %d ",cSection,&nElemsWall,&nElems);
    fscanf(fIn,"%s %s",cElemType,cElemFile);
    for(i = 0; i < 6; i++) cElemType[i] = tolower(cElemType[i]);

    /* Read problem section */
    fscanf(fIn,"%s",cSection);
    size = 3*nNodes;                 /* Characteristic size of the problem   */
    vProbParam = doubleVector(1,1);  /* Could be a double scalar - viscosity */
    vBCType = uintVector(nNodes,1);
    mBC = doubleMatrix(nNodes,3,1);
    fscanf(fIn,"%le %le %le %le",&vProbParam[0],&maxU,&Ly,&Lz);
    fscanf(fIn,"%s",cBCFile);

    /* Read Particle data section */
    fscanf(fIn,"%s",cSection);
    fscanf(fIn,"%le %le %le",&Xp[0],&Xp[1],&Xp[2]); // Position of the center of mass
    fscanf(fIn,"%le %le %le",&Vp[0],&Vp[1],&Vp[2]); // Velocity of the center of mass
    fscanf(fIn,"%le %le %le",&Wp[0],&Wp[1],&Wp[2]); // Angular velocity about body axes
    fscanf(fIn,"%le %le",&mass,&Rp);                // mass and radius
    fscanf(fIn,"%le %d",&dt,&steps);                // Time step and number of steps

    /* Clean and read BC file */
    comFilter(cBCFile);
    fAux = fopen(cBCFile,"r");
    for(j = 0; j < nNodes; j++){
        fscanf(fAux,"%d",&nBuffer);     /* Discard index */
        fscanf(fAux,"%d",&vBCType[j]);      
        fscanf(fAux,"%le %le %le",&mBC[j][0],&mBC[j][1],&mBC[j][2]);
        if(j > nNodesWall){
            mBC[j][0] = Vp[0];  /* To avoid rounding error introduced */
            mBC[j][1] = Vp[1];  /* by patran when saving files        */
            mBC[j][2] = Vp[2];
        }
    }
    fclose(fAux);

    /* Set element properties (nDim, nNodesInElem, nElemType) */
    nElemType = elemType();
    if((nElemType != 5) && (nElemType != 6))
        errorHandler("ERROR: incorrect element type");

    /* Read analysis section - default solver is Gauss-Jordan */
    fscanf(fIn,"%s %s ",cBuffer,solver);
    for(i = 0; i < strlen(solver); i++) solver[i] = toupper(solver[i]);
    if(strcmp(solver,"GMRES") == 0) /* Get initialization method for gmres */
        fscanf(fIn,"%d %d",&preCond,&nInit);
    fscanf(fIn," %d ",&ANALYSIS);
    
    /* Read internal points section (points where U and P are evaluated) */
    fscanf(fIn,"%s %d",cBuffer,&nInternalPoints);
    if(nInternalPoints != 0){
        fscanf(fIn,"%s",cInternalPointsFile);
        Xinner = doubleMatrix(nInternalPoints,nDim,0);
    } 
    fclose(fIn);

    /* Memory allocation for element connectivity and nodes cordinates */
    mNodes = doubleMatrix(nNodes,nDim,0);
    mElems = uintMatrix(nElems,nNodesInElem,0);

    /* Clean and read node file */
    comFilter(cNodeFile);
    fAux = fopen(cNodeFile,"r");
    for(i = 0; i < nNodes; i++){
        fscanf(fAux,"%d",&nBuffer);         /* Discard index */
        for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&mNodes[i][j]);
    }
    fclose(fAux);

    /* Clean and read element file */
    comFilter(cElemFile);
    fAux = fopen(cElemFile,"r");
    for(i = 0; i < nElems; i++){
    fscanf(fAux,"%d",&nBuffer);         /* Discard index */
        for(j = 0; j < nNodesInElem; j++) fscanf(fAux,"%d",&mElems[i][j]);
    }
    fclose(fAux);

    /* Clean and read internal points file */
    if(nInternalPoints != 0){
        comFilter(cInternalPointsFile);
        fAux = fopen(cInternalPointsFile,"r");
        for(i = 0; i < nInternalPoints; i++){
            fscanf(fAux,"%d",&nBuffer);
            for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&Xinner[i][j]);
        }
        fclose(fAux);
    }
    
    /* Use fully developed flow profile at inlet and outlet */
    /* Only active if the velocities are negative in mBC    */
    profileSetup3D(Ly,Lz,maxU,mNodes,mBC);

    /* Record program status at this stage in file log */
    fprintf(file_log,"\nUsing %d Nodes and %d %s ",nNodes,nElems,cElemType);
    fprintf(file_log,"Elems to solve 3D Stokes Flow Problem");
    fprintf(file_log,"\nUsing %s solver.",solver);
    t1 = (double)clock()/CLOCKS_PER_SEC;

    /****************************************************************/
    /***     INPUT DATA OBTAINED, SOLUTION OF PROBLEM FOLLOWS     ***/
    /****************************************************************/
 
    /* Memory allocation for coefficient matrix and right hand side vector */
    mA = doubleMatrix(size,size,0);
    vB = doubleVector(size,0);
    
    /* Save initial configuration and geometry */
    runtime = 0.0;
    for(i = 0; i < 3; i++){
        Fp[i] = 0.0;
        Tp[i] = 0.0;
        Lp[i] = 0.0;
    }
    fAux = fopen("update.dat","w");
    fprintf(fAux,"#it\ttime\tXp[0]\tXp[1]\tXp[2]\tVp[0]\tVp[1]\tVp[2]\t");
    fprintf(fAux,"Wp[0]\tWp[1]\tWp[2]\tFx\tFy\tFz\tTx\tTy\tTz\n");
    fprintf(fAux,"%d:\t%le\t",(int)runtime,runtime);
    fprintf(fAux,"%le\t%le\t%le\t",Xp[0],Xp[1],Xp[2]);
    fprintf(fAux,"%le\t%le\t%le\t",Vp[0],Vp[1],Vp[2]);
    fprintf(fAux,"%le\t%le\t%le\t",Wp[0],Wp[1],Wp[2]);
    fprintf(fAux,"%le\t%le\t%le\t",Fp[0],Fp[1],Fp[2]);
    fprintf(fAux,"%le\t%le\t%le\n",Tp[0],Tp[1],Tp[2]);
    fclose(fAux);
    
    fAux = fopen("particle.dat","w");
    fprintf(fAux,"#X\tY\n");
    for(i = 0; i < nNodes; i++)
        fprintf(fAux,"%le\t%le\t%le\n",mNodes[i][0],mNodes[i][1],mNodes[i][2]);
    fprintf(fAux,"\n");
    fclose(fAux);

fprintf(file_log,"\nInitiating main loop.");
fclose(file_log); file_log = fopen("bem.log","a");

    /* Iterate for the indicated number of steps */
    for(i = 1; i <= steps; i++){
    
        runtime = runtime + dt;
    
        /* Form coefficient matrix */
        if(nElemType == 5){
            fprintf(file_log,"\ndepStokes(): calling formMatrix_tria3() ...");
            stokesFormA_tria3(mNodes,mElems,vBCType,vProbParam,mA,mBC,vB);
        }
        else{
            fprintf(file_log,"\ndepStokes(): calling formMatrix_tria6() ...");
fclose(file_log); file_log = fopen("bem.log","a");
            stokesFormA_tria6(mNodes,mElems,vBCType,vProbParam,mA,mBC,vB);
        }

        /* Solve linear system */
        if(strcmp(solver,"GMRES") == 0){
            fprintf(file_log,"\ndepStokes(): calling solverGMRES() ...");
fclose(file_log); file_log = fopen("bem.log","a");
            solverGMRES(preCond,nInit,size,mA,vB);
        }
        else{
            fprintf(file_log,"\ndepStokes(): calling gaussBksb() ...");
            gaussBksb(size,mA,vB);
        }
        
        /* Get force and torque for this step */
        if(nElemType == 5){
            force_tria3(mNodes,mElems,vB,Fp);
            torque_tria3(mNodes,mElems,vB,Xp,Tp);
        }
        else{
            fprintf(file_log,"\ndepStokes(): calculating force and torque ...\n");
            force_tria6(mNodes,mElems,vB,Fp);
            torque_tria6(mNodes,mElems,vB,Xp,Tp);
        }
    
        /* Update all quantitites */
        update(mNodes,mBC,Xp,Vp,Wp,Lp,mass,Rp,Fp,Tp);

        fAux = fopen("update.dat","a");
        fprintf(fAux,"%d:\t%le\t",i,runtime);
        fprintf(fAux,"%le\t%le\t%le\t",Xp[0],Xp[1],Xp[2]);
        fprintf(fAux,"%le\t%le\t%le\t",Vp[0],Vp[1],Vp[2]);
        fprintf(fAux,"%le\t%le\t%le\t",Wp[0],Wp[1],Wp[2]);
        fprintf(fAux,"%le\t%le\t%le\t",Fp[0],Fp[1],Fp[2]);
        fprintf(fAux,"%le\t%le\t%le\n",Tp[0],Tp[1],Tp[2]);
        fclose(fAux);

        fAux = fopen("particle.dat","a");
        for(j = nNodesWall; j < nNodes; j++)
            fprintf(fAux,"%le\t%le\t%le\n",mNodes[j][0],mNodes[j][1],mNodes[j][2]);
        fprintf(fAux,"\n");
        fclose(fAux);
    }
    
    /**************************************************************/
    /***      SOLUTION OBTAINED, POST-PROCESSING FOLLOWS        ***/
    /**************************************************************/

    /* Calculate velocity and pressure at required points */
    if(nInternalPoints != 0){
        if(nElemType == 5){
            fprintf(file_log,"\ndepStokes(): calling postProcess_tria3() ...");
            stokesPostProcess_tria3(ANALYSIS,Xinner,mNodes,mElems,vProbParam,vBCType,vB);
        }
        else{
            fprintf(file_log,"\ndepStokes(): calling postProcess_tria6() ...");
            stokesPostProcess_tria6(ANALYSIS,Xinner,mNodes,mElems,vProbParam,vBCType,vB);
        }
    }
    t2 = (double)clock()/CLOCKS_PER_SEC-t1;

    /*******************************************************************/
    /***  POST-PROCESSING FINISHED, ANALYSIS OF PERFORMANCE FOLLOWS  ***/
    /*******************************************************************/

    fprintf(file_log,"\n\n*** depStokes(): calling performance analysis ***");
    if(t1 > 60){
        t1m = (int)(t1/60.0);
        t1s = t1-t1m*60.0;
        fprintf(file_log,"\nTime reading input: \t\t%d m %2.1f s",t1m,t1s);
    }
    else fprintf(file_log,"\nTime reading input: \t\t%2.1f seconds",t1);
    if(t2 > 60*steps){
        t2m = (int)(t2/(60.0*steps));
        t2s = (t2-t2m*60*steps)/steps;
        fprintf(file_log,"\nTime per iteration: \t\t%d m %2.1f s",t2m,t2s);
    }
    else fprintf(file_log,"\nTime per iteration: \t\t%2.1f seconds",t2/steps);
    t_tot = t1+t2;
    if(t_tot > 60){
        t_totm = (int)(t_tot/60.0);
        t_tots = t_tot-t_totm*60.0;
        fprintf(file_log,"\nTotal execution time: \t%d m %2.1f s",t_totm,t_tots);
    }
    else fprintf(file_log,"\nTotal execution time: \t\t%2.1f seconds",t_tot);

    /************************************************************************/
    /***  PERFORMANCE ANALYSIS FINISHED, CLOSE ALL FILES AND FREE MEMORY  ***/
    /************************************************************************/

    free(vB);
    free(vBCType);
    free(vProbParam);
    freeUintMatrix(mElems,nElems);
    freeDoubleMatrix(mBC,nNodes);
    freeDoubleMatrix(mA,size);
    freeDoubleMatrix(mNodes,nNodes);
    if(nInternalPoints != 0) freeDoubleMatrix(Xinner,nInternalPoints);
    
    fprintf(file_log,"\n\ndepStokes(): finished @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);
    fprintf(file_log,"\n\n**********************************");
    fprintf(file_log,"************************\n\n");
    fclose(file_log);

    return 0;
}
