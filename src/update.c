/******************************************************************************
* Name      : update.c                                                        *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Updates particle surface nodes coordinates and velocity. Fp and Tp are      *
* inputs and remain unchanged on exit. First order scheme following Russ.     *
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
#include <math.h>

extern double dt;
extern unsigned int nNodesWall, nNodes;

int update(double **mNodes, double **mBC, double *Xp, double *Vp, double *Wp,
           double *Lp, double mass, double Rp, double *Fp, double *Tp)
{
    unsigned int i;
    double dvx, dvy, dvz, dx, dy, dz, wx, wy, wz, Iinv;
    double cw[3], sw[3], wrkspc[2];
    
    /* Initialization */
    Iinv = 2.5/(mass*Rp*Rp);
    dvx = Fp[0]/mass;
    dvy = Fp[1]/mass;
    dvz = Fp[2]/mass;
    if(fabs(dvx) < 1.0E-15) dvx = 0.0;
    if(fabs(dvy) < 1.0E-15) dvy = 0.0;
    if(fabs(dvz) < 1.0E-15) dvz = 0.0;
    
    /* Update angular velocity first */
    for(i = 0; i < 3; i++){
        Wp[i] = Lp[i]*Iinv;
        /*   For second order update scheme use this     */
        /* Wp[i] = Lp[i]*Iinv*(1.0 + 0.5*Tp[i]*Iinv*dt); */
        cw[i] = cos(Wp[i]*dt);
        sw[i] = sin(Wp[i]*dt);
    }
    
    /* Update velocity and position of nodes in the rigid body surface */
    for(i = nNodesWall; i < nNodes; i++){
    
        /* Update the speed of all nodes */
        dx = mNodes[i][0] - Xp[0];
        dy = mNodes[i][1] - Xp[1];
        dz = mNodes[i][2] - Xp[2];    
        wx = dz*Wp[1] - dy*Wp[2];    // Must be done with old positions for
        wy = dx*Wp[2] - dz*Wp[0];    // both center of masss and surface point
        wz = dy*Wp[0] - dx*Wp[1];        
        mBC[i][0] = Vp[0] + dvx*dt + wx;
        mBC[i][1] = Vp[1] + dvy*dt + wy;
        mBC[i][2] = Vp[2] + dvz*dt + wz;
        
        /* Apply rotation to the nodes */
        mNodes[i][0] = mNodes[i][0]*cw[1]*cw[2] + mNodes[i][1]*cw[1]*sw[2] - mNodes[i][2]*sw[1];
        wrkspc[0] = sw[0]*sw[1]*cw[2] - cw[0]*sw[2];
        wrkspc[1] = cw[0]*cw[2] + sw[0]*sw[1]*sw[2];
        mNodes[i][1] = mNodes[i][0]*wrkspc[0] + mNodes[i][1]*wrkspc[1] + mNodes[i][2]*sw[0]*cw[1];
        wrkspc[0] = sw[0]*sw[2] + cw[0]*sw[1]*cw[2];
        wrkspc[1] = cw[0]*sw[1]*sw[2] - sw[0]*cw[2];
        mNodes[i][2] = mNodes[i][0]*wrkspc[0] + mNodes[i][1]*wrkspc[1] + mNodes[i][2]*cw[0]*cw[1];
        
        /* And increase their velocity to match the center of mass */
        mNodes[i][0] = mNodes[i][0] + Vp[0]*dt;
        mNodes[i][1] = mNodes[i][1] + Vp[1]*dt;
        mNodes[i][2] = mNodes[i][2] + Vp[2]*dt;
    } 
 
    /* Update position and velocity of center of mass */
    Xp[0] = Xp[0] + Vp[0]*dt;
    Xp[1] = Xp[1] + Vp[1]*dt;
    Xp[2] = Xp[2] + Vp[2]*dt;
    Vp[0] = Vp[0] + dvx*dt;
    Vp[1] = Vp[1] + dvy*dt;
    Vp[2] = Vp[2] + dvz*dt;
    
    /* Update angular momentum */
    if(fabs(Tp[0]*Iinv) < 1.0E-15) Tp[0] = 0.0;
    else Lp[0] = Lp[0] + Tp[0]*dt;
    if(fabs(Tp[1]*Iinv) < 1.0E-15) Tp[1] = 0.0;
    else Lp[1] = Lp[1] + Tp[1]*dt;
    if(fabs(Tp[2]*Iinv) < 1.0E-15) Tp[2] = 0.0;
    else Lp[2] = Lp[2] + Tp[2]*dt;

    return 0;
}
