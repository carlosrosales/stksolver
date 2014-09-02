/******************************************************************************
* File      : p2b-stk.c                                                 *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Creates the *.bem fields from the *.out inputs exported from Patran using   *
* bem_save_stk(). If a node has two different boundary conditions (with one   *
* bc for interface and one for source) only the bc for the source is kept.    *
* Number of nodes, elements and boundary conditions is printed to the file    *
* "mesh-info.dat". This is the version for Stokes.                            *
*                                                                             *
* NODES       : Number of nodes in file "nodes.out"                           *
* ELEMS       : Number of elements in fiel "elems.out"                        *
* ELEM_TYPE   : Type of element (tria3, tria6, quad4 or quad8)                *
* SCALING     : Scaling factor by wich the coordinate values are multiplied   *
* BCS_NUMBER  : Number of lines in file "bcs.out"                             *
* REORDER     : Need to reorder nodes in the elements? [0 == no | 1 == yes]   *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of p2b-stk.                                               *
*                                                                             *
* p2b-stk is free software; you can redistribute it and/or modify             *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* p2b-stk is distributed in the hope that it will be useful,                  *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with p2b-stk; if not, write to the Free Software                      *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "errorHandler.c"
#include "doubleMatrix.c"
#include "freeDoubleMatrix.c"
#include "freeUintMatrix.c"
#include "uintMatrix.c"
#include "uintVector.c"

int main(int argc, char *argv[])
{
    FILE    *bcs_in, *bcs_out, *elems_in, *elems_out, *nodes_in, *nodes_out;
    unsigned int    BCS_NUMBER, count, COLS, ROWS, i, id, j, k, m, NODES, ELEMS, 
                    FLAG, NODES_IN_ELEM, REORDER;
    unsigned int    **elemMat, *indx, *newIndx;
    double          SCALING_FACTOR, x, y, z;
    double          a[4][5], bc[4], **M1, **M2, **nodeMat;

    /* Check number and validity of arguments */
    if(argc == 2 && strcmp(argv[1],"-h") == 0){
        printf("\nCall as:\n\tp2b NODES ELEMS ELEM_TYPE SCALING BCS_NUMBER REORDER\n\n");
        printf("Creates the *.bem fields from the *.out inputs exported from Patran using\n");
        printf("bem_save_u(). This is the version for the stokes case.\n\n");
        printf("NODES\t: Number of nodes in file 'nodes.out'\n");
        printf("ELEMS\t: Number of elements in fiel 'elems.out'\n");
        printf("ELEM_TYPE\t: Type of element (tria3, tria6, quad4 or quad8)\n");
        printf("SCALING\t: Scaling factor by wich the coordinate values will be multiplied\n");
        printf("BCS_NUMBER\t: Number of lines in file 'bcs.out'\n");
        printf("REORDER\t: Need to reorder nodes in the elements? [0 == no | 1 == yes]\n\n");
        exit(0);
    }
    else if(argc != 7){
        printf("Error: Incorrect number of arguments!\nCorrect syntaxis is:\n\n");
        printf("p2b NODES ELEMS ELEM_TYPE SCALING_FACTOR BCS_NUMBER REORDER\n\n");
        errorHandler("Type './p2b -h' for help.\n\n");
    }
    else{
        NODES = atoi(argv[1]);
        ELEMS = atoi(argv[2]);
        SCALING_FACTOR = atof(argv[4]);
        BCS_NUMBER = atoi(argv[5]);
        REORDER = atoi(argv[6]);
    }
    if( (NODES < 0) || (ELEMS < 0) )
        errorHandler("Error: Number of nodes and elements must be positive!!");
        
    /* Check that files can be open */
    if((bcs_in = fopen("bcs.out","r")) == NULL)
        errorHandler("Error: Can't open input file bcs.out");
    if((elems_in = fopen("elems.out","r")) == NULL)
        errorHandler("Error: Can't open input file elems.out");
    if((nodes_in = fopen("nodes.out","r")) == NULL)
        errorHandler("Error: Can't open input file nodes.out");
    if((bcs_out = fopen("bcs.bem","w")) == NULL)
        errorHandler("Error: Can't open output file bcs.bem");
    if((elems_out = fopen("elems.bem","w")) == NULL)
        errorHandler("Error: Can't open output file elems.bem");
    if((nodes_out = fopen("nodes.bem","w")) == NULL)
        errorHandler("Error: Can't open output file nodes.bem");

    /* Get number of nodes in element */ 
    if(strcmp(argv[3],"tria3") == 0) NODES_IN_ELEM = 3;
    else if(strcmp(argv[3],"tria6") == 0) NODES_IN_ELEM = 6;
    else if(strcmp(argv[3],"quad4") == 0) NODES_IN_ELEM = 4;
    else if(strcmp(argv[3],"quad8") == 0) NODES_IN_ELEM = 8;
    else {
            fclose(nodes_in);
            fclose(elems_in);
            fclose(bcs_in);
            fclose(nodes_out);
            fclose(elems_out);
            fclose(bcs_out);
            errorHandler("Error: Incorrect element type");
    }
        
    /* Memory for temporary storage */
    COLS = 5;
    M1 = doubleMatrix(BCS_NUMBER,COLS,1);
    M2 = doubleMatrix(NODES,COLS,1);
    nodeMat = doubleMatrix(NODES,3,1);
    elemMat = uintMatrix(ELEMS,NODES_IN_ELEM,1);
    
    /* Read node and element data into temporary storage */
    for(i = 0; i < NODES; i++){
        fscanf(nodes_in,"%d %le ",&id,&nodeMat[i][0]);
        fscanf(nodes_in,"%le %le",&nodeMat[i][1],&nodeMat[i][2]);
    }
    for(i = 0; i < ELEMS; i++){
        fscanf(elems_in,"%d",&id);
        for(j = 0; j < NODES_IN_ELEM; j++) fscanf(elems_in,"%d",&elemMat[i][j]);
    }
    fclose(nodes_in);
    fclose(elems_in);
    
        
    /************************************************************/
    /***    WRITE SCALED NODE FILE WITHOUT REPEATED NODES     ***/
    /************************************************************/
    
    count = 0;
    for(i = 0; i < NODES; i++){
        x = nodeMat[i][0]*SCALING_FACTOR;
        y = nodeMat[i][1]*SCALING_FACTOR;
        z = nodeMat[i][2]*SCALING_FACTOR;
        fprintf(nodes_out,"%d\t%le\t%le\t%le\n",++count,x,y,z);
    }
    fclose(nodes_out);
    
    /************************************************************/
    /*** CONVERT CONNECTIVITY FILE TO ACCEPTABLE INPUT FORMAT ***/
    /************************************************************/
    
    if(REORDER == 1){
        for(i = 0; i < ELEMS; i++){
            if(NODES_IN_ELEM == 3){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\n",elemMat[i][2],elemMat[i][0]);
            }
            else if(NODES_IN_ELEM = 6){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][4],elemMat[i][2]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][5],elemMat[i][0]);
                fprintf(elems_out,"%d\n",elemMat[i][3]);
            }
            else if(NODES_IN_ELEM == 4){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][2],elemMat[i][3]);
                fprintf(elems_out,"%d\n",elemMat[i][0]);
            }
            else if(NODES_IN_ELEM == 8){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][5],elemMat[i][2]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][6],elemMat[i][3]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][7],elemMat[i][0]);
                fprintf(elems_out,"%d\n",elemMat[i][4]);
            }
        }
    }
    else{
        for(i = 0; i < ELEMS; i++){
            fprintf(elems_out,"%d",i+1);
            for(j = 0; j < NODES_IN_ELEM; j++){
                fprintf(elems_out,"\t%d",elemMat[i][j]);
            }
            fprintf(elems_out,"\n");
        }
    }
    fclose(elems_out);

    /************************************************************/
    /***      CONVERT BCS FILE TO ACCEPTABLE INPUT FORMAT     ***/
    /************************************************************/
    
    for(i = 0; i < BCS_NUMBER; i++){
        fscanf(bcs_in,"%le %le ",&M1[i][0],&M1[i][1]);
        fscanf(bcs_in,"%le %le %le",&M1[i][2],&M1[i][3],&M1[i][4]);
    }
    fclose(bcs_in);
    
    for(i = 0; i < NODES; i++){    
        FLAG = 0;
        for(j = 0; j < BCS_NUMBER; j++){
            if(M1[j][0] == i+1){
                if(M1[j][1] == 0.0) bc[FLAG] = 0.0;         
                else bc[FLAG] = 1.0;
                a[FLAG][0] = i+1;
                a[FLAG][1] = M1[j][1];
                a[FLAG][2] = M1[j][2];
                a[FLAG][3] = M1[j][3];
                a[FLAG][4] = M1[j][4];
                FLAG = FLAG + 1;
            }
        }
        if(FLAG == 2){
            M2[i][0] = a[0][0];
            M2[i][1] = a[0][1];
            M2[i][2] = a[0][2];
            M2[i][3] = a[0][3];
            M2[i][4] = a[0][4];     
        }else if(FLAG == 4){
            if(bc[0] == 1.0){
                M2[i][0] = a[0][0];
                M2[i][1] = a[0][1];
                M2[i][2] = a[0][2];
                M2[i][3] = a[0][3]; 
                M2[i][4] = a[0][4];  
            }else{
                M2[i][0] = a[2][0];
                M2[i][1] = a[2][1];
                M2[i][2] = a[2][2];
                M2[i][3] = a[2][3];
                M2[i][4] = a[2][4];  
            }
        }
        else errorHandler("Incorrect FLAG value");
    }

    /* Write bcs output file (only for non-repeated nodes) */
    for(i = 0; i < NODES; i++){
        fprintf(bcs_out,"%d\t%d\t",(int)M2[i][0],(int)M2[i][1]);
        fprintf(bcs_out,"%le\t%le\t%le\n",M2[i][2],M2[i][3],M2[i][4]);
    }
    fclose(bcs_out);
    
    freeDoubleMatrix(M1,BCS_NUMBER);
    freeDoubleMatrix(M2,NODES);
    freeDoubleMatrix(nodeMat,NODES);
    freeUintMatrix(elemMat,ELEMS);

    return 0;
}

 
