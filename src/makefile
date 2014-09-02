# Makefile for stkSolver
#
# COPYRIGHT & LICENSE INFORMATION
# 
# Copyright 2006 Carlos Rosales Fernandez and The Institute of High
# Performance Computing (A*STAR)
#
# This file is part of stkSolver.
#
# stkSolver is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# stkSolver is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with stkSolver; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

SHELL = /bin/sh

# Compiler
CC = gcc

# Flags
LFLAGS = -lm -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3
CFLAGS = -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3

# Linking Step
stkSolver: comFilter.o dotProd.o doubleMatrix.o doublePointer.o doubleVector.o\
     elemType.o errorHandler.o force_tria3.o force_tria6.o freeDoubleMatrix.o \
     freeDoublePointer.o freeUintMatrix.o gaussBksb.o getLocalNormal_tria3.o  \
     getNormal_tria3.o initRes.o intGStk_tria3.o intHStk_tria3.o              \
     intSingularGStk_tria3.o intSingularHStk_tria3.o iterGMRES.o              \
     L2Norm.o matVectProd.o stokesPostProcess_tria3.o pressure_tria3.o        \
     shape_tria3.o solverGMRES.o uintMatrix.o uintVector.o X2L_tria3.o        \
     X2L_tria6.o stokesFormA_tria6.o getLocalNormal_tria6.o getNormal_tria6.o \
     intGStk_tria6.o intHStk_tria6.o intSingularGStk_tria6.o                  \
     intSingularHStk_tria6.o stokesPostProcess_tria6.o pressure_tria6.o       \
     shape_tria6.o update.o velocity_tria3.o velocity_tria6.o profileSetup3D.o\
     torque_tria3.o torque_tria6.o stokes3d.o stokesFormA_tria3.o 
	$(CC) $(LFLAGS) comFilter.o dotProd.o doubleMatrix.o doublePointer.o      \
     elemType.o errorHandler.o force_tria3.o force_tria6.o freeDoubleMatrix.o \
     freeDoublePointer.o freeUintMatrix.o gaussBksb.o getLocalNormal_tria3.o  \
     getNormal_tria3.o initRes.o intGStk_tria3.o intHStk_tria3.o stokes3d.o   \
     intSingularGStk_tria3.o intSingularHStk_tria3.o iterGMRES.o              \
     L2Norm.o matVectProd.o stokesPostProcess_tria3.o pressure_tria3.o        \
     shape_tria3.o solverGMRES.o uintMatrix.o uintVector.o X2L_tria3.o        \
     X2L_tria6.o stokesFormA_tria6.o getLocalNormal_tria6.o getNormal_tria6.o \
     intGStk_tria6.o intHStk_tria6.o intSingularGStk_tria6.o                  \
     intSingularHStk_tria6.o stokesPostProcess_tria6.o pressure_tria6.o       \
     shape_tria6.o update.o velocity_tria3.o velocity_tria6.o profileSetup3D.o\
     torque_tria3.o torque_tria6.o doubleVector.o stokesFormA_tria3.o -o stkSolver
	mv stkSolver ../bin/stkSolver 

stkSolver-sft: comFilter.o dotProd.o doubleMatrix.o doublePointer.o doubleVector.o\
     elemType.o errorHandler.o force_tria3.o force_tria6.o freeDoubleMatrix.o \
     freeDoublePointer.o freeUintMatrix.o gaussBksb.o getLocalNormal_tria3.o  \
     getNormal_tria3.o initRes.o intGStk_tria3.o intHStk_tria3.o              \
     intSingularGStk_tria3.o intSingularHStk_tria3.o iterGMRES.o              \
     L2Norm.o matVectProd.o stokesPostProcess_tria3.o pressure_tria3.o        \
     shape_tria3.o solverGMRES.o uintMatrix.o uintVector.o X2L_tria3.o        \
     X2L_tria6.o stokesFormA_tria6.o getLocalNormal_tria6.o getNormal_tria6.o \
     intGStk_tria6.o intHStk_tria6.o intSingularGStk_tria6.o                  \
     intSingularHStk_tria6.o stokesPostProcess_tria6.o pressure_tria6.o       \
     shape_tria6.o update.o velocity_tria3.o velocity_tria6.o profileSetup3D.o\
     torque_tria3.o torque_tria6.o stokes3d-shift.o stokesFormA_tria3.o
	$(CC) $(LFLAGS) comFilter.o dotProd.o doubleMatrix.o doublePointer.o      \
     elemType.o errorHandler.o force_tria3.o force_tria6.o freeDoubleMatrix.o \
     freeDoublePointer.o freeUintMatrix.o gaussBksb.o getLocalNormal_tria3.o  \
     getNormal_tria3.o initRes.o intGStk_tria3.o intHStk_tria3.o stokes3d-shift.o\
     intSingularGStk_tria3.o intSingularHStk_tria3.o iterGMRES.o              \
     L2Norm.o matVectProd.o stokesPostProcess_tria3.o pressure_tria3.o        \
     shape_tria3.o solverGMRES.o uintMatrix.o uintVector.o X2L_tria3.o        \
     X2L_tria6.o stokesFormA_tria6.o getLocalNormal_tria6.o getNormal_tria6.o \
     intGStk_tria6.o intHStk_tria6.o intSingularGStk_tria6.o                  \
     intSingularHStk_tria6.o stokesPostProcess_tria6.o pressure_tria6.o       \
     shape_tria6.o update.o velocity_tria3.o velocity_tria6.o profileSetup3D.o\
     torque_tria3.o torque_tria6.o doubleVector.o stokesFormA_tria3.o -o stkSolver-sft
	mv stkSolver-sft ../bin/stkSolver-sft

# Compiling Step
comFilter.o: comFilter.c
	$(CC) -c $(CFLAGS) comFilter.c -o comFilter.o

dotProd.o: dotProd.c
	$(CC) -c $(CFLAGS) dotProd.c -o dotProd.o

doubleMatrix.o: doubleMatrix.c 
	$(CC) -c $(CFLAGS) doubleMatrix.c -o doubleMatrix.o

doublePointer.o: doublePointer.c
	$(CC) -c $(CFLAGS) doublePointer.c -o doublePointer.o

doubleVector.o: doubleVector.c
	$(CC) -c $(CFLAGS) doubleVector.c -o doubleVector.o

elemType.o: elemType.c
	$(CC) -c $(CFLAGS) elemType.c -o elemType.o

errorHandler.o: errorHandler.c
	$(CC) -c $(CFLAGS) errorHandler.c -o errorHandler.o

force_tria3.o: force_tria3.c force_tria3.h 
	$(CC) -c $(CFLAGS) force_tria3.c -o force_tria3.o

force_tria6.o: force_tria6.c force_tria6.h integral_tria6.h
	$(CC) -c $(CFLAGS) force_tria6.c -o force_tria6.o

freeDoubleMatrix.o: freeDoubleMatrix.c
	$(CC) -c $(CFLAGS) freeDoubleMatrix.c -o freeDoubleMatrix.o

freeDoublePointer.o: freeDoublePointer.c
	$(CC) -c $(CFLAGS) freeDoublePointer.c -o freeDoublePointer.o

freeUintMatrix.o: freeUintMatrix.c
	$(CC) -c $(CFLAGS) freeUintMatrix.c -o freeUintMatrix.o

gaussBksb.o: gaussBksb.c
	$(CC) -c $(CFLAGS) gaussBksb.c -o gaussBksb.o

getLocalNormal_tria3.o: getLocalNormal_tria3.c
	$(CC) -c $(CFLAGS) getLocalNormal_tria3.c -o getLocalNormal_tria3.o

getLocalNormal_tria6.o: getLocalNormal_tria6.c
	$(CC) -c $(CFLAGS) getLocalNormal_tria6.c -o getLocalNormal_tria6.o

getNormal_tria3.o: getNormal_tria3.c
	$(CC) -c $(CFLAGS) getNormal_tria3.c -o getNormal_tria3.o

getNormal_tria6.o: getNormal_tria6.c
	$(CC) -c $(CFLAGS) getNormal_tria6.c -o getNormal_tria6.o

initRes.o: initRes.c initRes.h
	$(CC) -c $(CFLAGS) initRes.c -o initRes.o

intGStk_tria3.o: intGStk_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intGStk_tria3.c -o intGStk_tria3.o

intGStk_tria6.o: intGStk_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intGStk_tria6.c -o intGStk_tria6.o

intHStk_tria3.o: intHStk_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intHStk_tria3.c -o intHStk_tria3.o

intHStk_tria6.o: intHStk_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intHStk_tria6.c -o intHStk_tria6.o

intSingularGStk_tria3.o: intSingularGStk_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularGStk_tria3.c -o intSingularGStk_tria3.o

intSingularGStk_tria6.o: intSingularGStk_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularGStk_tria6.c -o intSingularGStk_tria6.o

intSingularHStk_tria3.o: intSingularHStk_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularHStk_tria3.c -o intSingularHStk_tria3.o

intSingularHStk_tria6.o: intSingularHStk_tria6.c integral_tria6.h
	$(CC) -c $(CFLAGS) intSingularHStk_tria6.c -o intSingularHStk_tria6.o

iterGMRES.o: iterGMRES.c iterGMRES.h
	$(CC) -c $(CFLAGS) iterGMRES.c -o iterGMRES.o

L2Norm.o: L2Norm.c
	$(CC) -c $(CFLAGS) L2Norm.c -o L2Norm.o

matVectProd.o: matVectProd.c
	$(CC) -c $(CFLAGS) matVectProd.c -o matVectProd.o

stokes3d.o: stokes3d.c stokes3d.h constants.h gaussData.h
	$(CC) -c $(CFLAGS) stokes3d.c -o stokes3d.o

stokes3d-shift.o: stokes3d-shift.c stokes3d.h constants.h gaussData.h
	$(CC) -c $(CFLAGS) stokes3d-shift.c -o stokes3d-shift.o

stokesFormA_tria3.o: stokesFormA_tria3.c stokesFormA_tria3.h constants.h
	$(CC) -c $(CFLAGS) stokesFormA_tria3.c -o stokesFormA_tria3.o

stokesFormA_tria6.o: stokesFormA_tria6.c stokesFormA_tria6.h constants.h
	$(CC) -c $(CFLAGS) stokesFormA_tria6.c -o stokesFormA_tria6.o

stokesPostProcess_tria3.o: stokesPostProcess_tria3.c stokesPostProcess_tria3.h
	$(CC) -c $(CFLAGS) stokesPostProcess_tria3.c -o stokesPostProcess_tria3.o

stokesPostProcess_tria6.o: stokesPostProcess_tria6.c stokesPostProcess_tria6.h
	$(CC) -c $(CFLAGS) stokesPostProcess_tria6.c -o stokesPostProcess_tria6.o

torque_tria3.o: torque_tria3.c torque_tria3.h integral_tria3.h
	$(CC) -c $(CFLAGS) torque_tria3.c -o torque_tria3.o

torque_tria6.o: torque_tria6.c torque_tria6.h integral_tria6.h
	$(CC) -c $(CFLAGS) torque_tria6.c -o torque_tria6.o

pressure_tria3.o: pressure_tria3.c pressure_tria3.h constants.h
	$(CC) -c $(CFLAGS) pressure_tria3.c -o pressure_tria3.o

pressure_tria6.o: pressure_tria6.c pressure_tria6.h constants.h
	$(CC) -c $(CFLAGS) pressure_tria6.c -o pressure_tria6.o

profileSetup3D.o: profileSetup3D.c
	$(CC) -c $(CFLAGS) profileSetup3D.c -o profileSetup3D.o

shape_tria3.o: shape_tria3.c
	$(CC) -c $(CFLAGS) shape_tria3.c -o shape_tria3.o

shape_tria6.o: shape_tria6.c
	$(CC) -c $(CFLAGS) shape_tria6.c -o shape_tria6.o

solverGMRES.o: solverGMRES.c solverGMRES.h
	$(CC) -c $(CFLAGS) solverGMRES.c -o solverGMRES.o

uintMatrix.o: uintMatrix.c
	$(CC) -c $(CFLAGS) uintMatrix.c -o uintMatrix.o

uintVector.o: uintVector.c
	$(CC) -c $(CFLAGS) uintVector.c -o uintVector.o

update.o: update.c
	$(CC) -c $(CFLAGS) update.c -o update.o

velocity_tria3.o: velocity_tria3.c velocity_tria3.h constants.h
	$(CC) -c $(CFLAGS) velocity_tria3.c -o velocity_tria3.o

velocity_tria6.o: velocity_tria6.c velocity_tria6.h constants.h
	$(CC) -c $(CFLAGS) velocity_tria6.c -o velocity_tria6.o

X2L_tria3.o: X2L_tria3.c
	$(CC) -c $(CFLAGS) X2L_tria3.c -o X2L_tria3.o

X2L_tria6.o: X2L_tria6.c
	$(CC) -c $(CFLAGS) X2L_tria6.c -o X2L_tria6.o

# Eliminate all object files
clean:
	rm -f *.o

