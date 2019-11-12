#include "fields.h"
REAL* create1Dfield(int size)
{
    REAL* new_real = (REAL*)malloc(size * sizeof(REAL));
    return new_real;
}

REAL** create2Dfield(int sizeX, int sizeY)
{
    REAL** new_real = (REAL**)malloc(sizeX*sizeof(REAL*));
    for(int i=0;i<sizeX;++i){
    new_real[i] = (REAL*)malloc(sizeY*sizeof(REAL));
    }
    return new_real;
}

REAL* createVector(int len)
{
    return create1Dfield(len);
}

REAL** createMatrix(int rows, int cols)
{
    return create2Dfield(rows, cols);
}

void destroy1Dfield(REAL* field)
{
    free (field);
}

void destroy2Dfield(REAL** field, int sizeX)
{
    for(int i=0; i<sizeX; ++i){;
        free(field[i]);
    }
    free(field);
}
int **Create2DIntegerField(int imax, int jmax)
{
    int** field;
    printf("1\n");
    field=malloc(( imax)*sizeof(int*));
    printf("2\n");
    for(int i=0; i< imax; i++){
        field[i]=malloc(( jmax)*sizeof(int));
    }
    return field;
}
void destroyVector(REAL* vector)
{
    destroy1Dfield(vector);
}

void destroyMatrix(REAL** matrix, int rows)
{
    destroy2Dfield(matrix, rows);
}

void fill1Dfield(REAL value, REAL* field, int size)
{
    for(int i=0; i<size; ++i) {
        field[i] = value;
    }
}

void fill2Dfield(REAL value, REAL** field, int sizeX, int sizeY)
{
    for(int i=0; i<sizeX; ++i) {
        fill1Dfield(value,field[i],sizeY);
    }
}

int isEqualScalar(REAL x, REAL y, REAL eps)
{
    if(fabs(x-y)>eps){
        return 0;
    }
    else return 1;
}

int isEqual1Dfield(REAL* field1, REAL* field2, int size, REAL eps)
{
    for(int i=0; i<size; ++i){
        if(fabs(field1[i]-field2[i])>eps){
            return 0;
        }
    }
    return 1;
}

int isEqual2Dfield(REAL** field1, REAL** field2, int sizeX, int sizeY, REAL eps)
{
    for(int i=0; i<sizeX; ++i){
        for(int j=0; j<sizeY; ++j){
            if(fabs(field1[i][j]-field2[i][j])>eps){
                return 0;
            }
        }
    }
    return 1;
}

void	applyFunctionTo1Dfield(REAL (*func)(REAL), REAL* field, int size)
{
    for(int i=0; i<size; ++i){
        field[i]=func(field[i]);
    }
}

void	applyFunctionTo2Dfield(REAL (*func)(REAL), REAL** field, int sizeX, int sizeY)
{
    for(int i=0; i<sizeX; ++i){
        applyFunctionTo1Dfield(func, field[i],sizeY);
        }
}

void	print1Dfield(REAL* field, int size)
{
    for(int i=0; i<size; ++i){
        printf("%lf",field[i]);
    }
}

void	print2Dfield(REAL** field, int sizeX, int sizeY)
{
    for(int j=0; j<sizeX; ++j){
        for(int i=0; i<sizeY; ++i){
            printf(" %lf ",field[j][i]);
        }
        printf(" \n");
    }
}
void	print1Dfield_int(int* field, int size)
{
    for(int i=0; i<size; ++i){
        printf("%d",field[i]);
    }
}

void	print2Dfield_int(int** field, int sizeX, int sizeY)
{
    for(int j=0; j<sizeX; ++j){
        for(int i=0; i<sizeY; ++i){
            printf(" %d ",field[j][i]);
        }
        printf(" \n");
    }
}

void	printVector(REAL* vector, int len)
{
    print1Dfield(vector,len);
}

void	printMatrix(REAL** matrix, int rows, int cols)
{
    print2Dfield(matrix,rows,cols);
}

long int toBinary(int n){
    if(n<2){
        return n;
    }
    else{
        return toBinary((n/2)*10)+(n%2);
    }
}
int toDecimal(long int binarynum)
{
    int decimalnum=0, temp=0, remainder;
    while(binarynum!=0)
    {
        remainder= binarynum%10;
        binarynum=binarynum/10;
        decimalnum=decimalnum + remainder*pow(2,temp);
        ++temp;
    }
    return decimalnum;
}

void write1Dfield(const char* fileName, REAL* field, int size)
{
    FILE *datei;
    datei=fopen(fileName,"w");
    fprintf(datei,"%ld",toBinary(size));
    for(int i=0; i<size; ++i){
        fprintf(datei,"%lf",field[i]);
    }
    fclose(datei);
}

void write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY)
{
    FILE *datei;
    datei=fopen(fileName,"w");
    if(datei!=NULL){
        fprintf(datei,"%ld",toBinary(sizeX));
        fprintf(datei,"%ld",toBinary(sizeY));
        for(int i=0; i<sizeX; ++i){
            for(int j=0; j<sizeY; ++j){
                fprintf(datei,"%lf",field[i][j]);
            }
        }
        fclose(datei);
    }
}

REAL* read1Dfield(const char* fileName, int* size)
{
    FILE *data=NULL;
    data=fopen(fileName,"rb");
    REAL*A;
    if(data==NULL){
        printf("Could not open file/n");
    }else{
        fscanf(data,"d%n",size);
        A=create1Dfield(*size);
        for(int i=0; i<*size; ++i){
            fscanf(data,"%lf",&A[i]);;
        }
    }
    fclose(data);
    return A;
}

REAL** read2Dfield(const char* fileName, int* sizeX, int* sizeY)
{
    FILE* data=NULL;
    data=fopen(fileName,"rb");
    if(data==NULL){
        printf("Could not open file/n");
        return NULL;
    }else{
        fscanf(data,"%d",sizeX);
        fscanf(data,"%d",sizeY);
        REAL**A=create2Dfield(*sizeX,*sizeY);
        for(int i=0;i<*sizeX;++i){
            for(int j=0;j<*sizeY;++j){
                fscanf(data,"%lf",&A[i][j]);
            }
        }
    fclose(data);
    return A;
    }
}

void writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "w");

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Scalar Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", sizeX, sizeY);
    fprintf(vtkFile, "X_COORDINATES %d double\n", sizeX);
    for (int i=0;i<sizeX;i++)
        fprintf(vtkFile, "%lf ", dx*(double)i);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Y_COORDINATES %d double\n", sizeY);
    for (int j=0;j<sizeY;j++)
        fprintf(vtkFile, "%lf ", dy*(double)j);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Z_COORDINATES 1 double\n");
    fprintf(vtkFile, "0.0\n");

    fprintf(vtkFile, "POINT_DATA %d\n", 1 * sizeX * sizeY);
    fprintf(vtkFile, "SCALARS %s double 1\n", description);
    fprintf(vtkFile, "LOOKUP_TABLE default \n");
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%e\n", field[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "w");

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Vector Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", sizeX, sizeY);
    fprintf(vtkFile, "X_COORDINATES %d double\n", sizeX);
    for (int i=0;i<sizeX;i++)
        fprintf(vtkFile, "%lf ", dx*(double)i);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Y_COORDINATES %d double\n", sizeY);
    for (int j=0;j<sizeY;j++)
        fprintf(vtkFile, "%lf ", dy*(double)j);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Z_COORDINATES 1 double\n");
    fprintf(vtkFile, "0.0\n");


    fprintf(vtkFile, "POINT_DATA %d\n", 1 * sizeX * sizeY);
    fprintf(vtkFile, "VECTORS %s double \n", description);
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%e %e 0.0\n", fieldU[i][j], fieldV[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DscalarField2(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "a");
    fprintf(vtkFile," \n");
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%e\n", field[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writePlotfile(REAL**U, int row, char*filename)
{
    FILE* file=fopen(filename,"a");
    for(int i=0;i<row;i++){
        fprintf(file,"%e\n",U[65][i]);
    }

}