/******************************************************************************
* File      : constants.h                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Contains all global constant definitions (except gauss quadrature data) and * 
* pre-procesor macros for the Stokes problem.                                 *
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

#define eps0	8.854187817E-12    /* Electric permittivity of vacuum (F/m)   */
#define mu0	    1.256637061E-06    /* Magnetic permittivity of vacuum (N/A^2) */
#define qe		1.602176531E-19         /* Elementary (electron) charge (C) */
#define pi		3.141592653589793E+00   /* Pi */
#define pi2		6.283185307179586E+00   /* 2*Pi */
#define pi4		1.256637061435917E+01   /* 4*Pi */
#define pio2	1.570796326794896E+00   /* Pi/2 */
#define pio4	7.853981633974483E-01   /* Pi/4 */
#define ip		3.183098861837907E-01   /* 1/Pi */
#define ip2		1.591549430918953E-01   /* 1/(2Pi) */
#define ip4		7.957747154594767E-02   /* 1/(4Pi) */

#define shift(a,b,c){double d = (a); (a) = (b); (b) = (c); (c) = d;}
#define swap(a,b){double c = (a); (a) = (b); (b) = c;}

typedef struct{
    /* Constant quantities */
    double   radius,         /* radius */
            mass,           /* mass */
            area,           /* surface area */
            lift,           /* lift force for non-neutrally buoyant bodies */
            IbodyInv[3][3]; /* local moment of inertia tensor inverse */
				
    /* State variables */
    double   r[3],           /* position vector */
            R[3][3],        /* orientation matrix giving the rotation */
            L[3];           /* angular momentum */
				
    /* Derived quantities */
    double   v[3],           /* velocity vector */
            omega[3],       /* angular velocity vector */
            Iinv[3][3];     /* moment of inertia inverse */
	
    /* Computed quantities */
    double   F[3],           /* force vector */
            T[3];           /* torque vector */
} RigidBody;

typedef struct{
    RigidBody *rb;          /* array of pointers to rigid bodies */
    int n;                  /* number of rigid bodies */
} RigidBodySystem;

