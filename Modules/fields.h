#ifndef FIELDS_H_
#define FIELDS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "real.h"

REAL*	create1Dfield(int size); //creates 1D field
REAL**	create2Dfield(int sizeX, int sizeY); //creates 2D field
int **Create2DIntegerField(int imax, int jmax);
REAL*	createVector(int len); //creates vector of length len
REAL**	createMatrix(int rows, int cols); //creates rowsxcols matrix

void 	destroy1Dfield(REAL* field); //frees field
void 	destroy2Dfield(REAL** field, int sizeX); //frees field
void	destroyVector(REAL* vector); //frees vector
void	destroyMatrix(REAL** matrix, int rows); //frees matrix

void	fill1Dfield(REAL value, REAL* field, int size); //fills the entrys of field with value
void	fill2Dfield(REAL value, REAL** field, int sizeX, int sizeY); //filles the entrys of field with value

int     isEqualScalar(REAL x, REAL y, REAL eps); //tests if |x-y|<eps
int 	isEqual1Dfield(REAL* field1, REAL* field2, int size, REAL eps); //tests if all entrys have distance <eps
int 	isEqual2Dfield(REAL** field1, REAL** field2, int sizeX, int sizeY,  REAL eps); //tests if all entrys have distance <eps

void	applyFunctionTo1Dfield(REAL (*func)(REAL), REAL* field, int size); //applys func to all entrys of field
void	applyFunctionTo2Dfield(REAL (*func)(REAL), REAL** field, int sizeX, int sizeY); //applys func to all entrys of field

void	print1Dfield(REAL* field, int size); //prints field
void	print2Dfield(REAL** field, int sizeX, int sizeY); //prints field
void	printVector(REAL* vector, int len); //prints vector
void	printMatrix(REAL** matrix, int rows, int cols); //prints matrix
void	print2Dfield_int(int** field, int sizeX, int sizeY);

void    writePlotfile(REAL**U, int row, char*filename);
void	write1Dfield(const char* fileName, REAL* field, int size); //writes size in binary and field into fileName
void	write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY); //writes sizeX, sizeY in binary and field into fileName

REAL* 	read1Dfield(const char* fileName, int* size); //reads 1D field from fileName
REAL**	read2Dfield(const char* fileName, int* sizeX, int* sizeY); //reads 2D field from fileName

void	writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy);
void	writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy);
void	writeVTKfileFor2DscalarField2(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy);

#endif /* FIELDS_H_ */
