#include "visual.h"

void OUTPUTVEC(REAL **U, REAL **V, REAL **P, int imax, int jmax, REAL delx, REAL dely, int n)
{
    /* S und T enthalten sind temporäre Arrays, die die Werte von P bzw. U und V enthalten, aber ohne die Randbedingungen */
    REAL **S = createMatrix(imax,jmax);
    REAL **T = createMatrix(imax,jmax);
    char fileName[64];
    /* Zuerst wird P in der Datei ./PresureField_n.vtk ausgegeben, wobei n eine beliebige ganze Zahl ist */
    if (P != NULL)
    {
        if (n != 0) sprintf(fileName,"data//PressureField_%i.vtk",n);
        else sprintf(fileName,"data//PressureField.vtk");
        for (int i = 0;  i < imax; i++)
            for (int j = 0; j < jmax; j++)
                T[i][j] = P[i+1][j+1];
        writeVTKfileFor2DscalarField(fileName,"pressurefield",T,imax,jmax,delx,dely);
    }
    if (U == NULL || V == NULL)
    {
        destroyMatrix(T,imax);
        destroyMatrix(S,imax);
        return;
    }
    /* Anschließend werden U und V in die Datei ./MomentumField_n.vtk ausgegeben mit dem n von oben */
    if (n != 0) sprintf(fileName,"data//MomentumField_%i.vtk",n);
    else sprintf(fileName,"data//MomentumField.vtk");
    for (int i = 1;  i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            T[i-1][j-1] = (U[i][j] + U[i-1][j])/2;
            S[i-1][j-1] = (V[i][j] + V[i][j-1])/2;
        }
    writeVTKfileFor2DvectorField(fileName,"momentumfield",T,S,imax,jmax,delx,dely);
    /* Abschließend werden T und S wieder freigegeben */
    destroyMatrix(T,imax);
    destroyMatrix(S,imax);
    return;
}
 

void PLOT(REAL **U,int imax, int row, char*fileName, REAL Re, char*problem)
{   
    printf("reday to plot\n");
    char* file;
    file=malloc(255*sizeof(char*));
    sprintf(file,"%s_%s_%f.ipynb",problem,fileName,Re);
    FILE *datei;
    datei = fopen (file, "a");
    if (datei == NULL)
    {
        printf("Fehler beim Oeffnen der Datei.");
    }
    fprintf(datei,"{\n");
    fprintf(datei," \"cells\": [\n");
    fprintf(datei,"  {\n");
    fprintf(datei,"   \"cell_type\": \"code\",\n");
    fprintf(datei,"   \"execution_count\": null,\n");
    fprintf(datei,"   \"metadata\": {},\n");
    fprintf(datei,"   \"outputs\": [],\n");
    fprintf(datei,"   \"source\": [\n");
    fprintf(datei,"    \"import matplotlib.pyplot as plt\\n\",\n");
    fprintf(datei,"    \"plt.plot([");
    for (int j=0;j< row;j++){
        fprintf(datei,"%lf",((double)j/(double)( row-1)));
        if(j< row-2+1){
            fprintf(datei,",");
        }
    }
    fprintf(datei,"], [");
    for (int j=0;j< row;j++){
        fprintf(datei,"%lf", U[imax/2][j]); 
        if(j< row-1){
            fprintf(datei,",");
        }
    }

    fprintf(datei,"])\\n\",\n"); //hier können wir graphen veraendern
    /* 
    if( Re==(double)100 && strcmp( problem,"dcavity")==0){
        fprintf(datei,"    \"plt.plot([1,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0], [1,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.2109,-0.15662,-0.1015,-0.06434,-0.04775,-0.04192,-0.03717,0], 'ro')\\n\",\n");

    }
    if( Re==(double)1000 && strcmp( problem,"dcavity")==0){
        fprintf(datei,"    \"plt.plot([1,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0], [1,0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,-0.0608,-0.10648,-0.27805,-0.38289,-0.2973,-0.2222,-0.20196,-0.18109,0], 'ro')\\n\",\n");

    }
    if( Re==(double)5000 && strcmp( problem,"dcavity")==0){
        fprintf(datei,"    \"plt.plot([1,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0], [1,0.48223,0.4612,0.45992,0.46036,0.33556,0.20087,0.08183,-0.03039,-0.07404,-0.22855,-0.3305,-0.40435,-0.43643,-0.42901,-0.41165,0], 'ro')\\n\",\n");

    }*/
    fprintf(datei,"    \"plt.ylabel('Geschwindigkeit u zu(reynoldzahl=%lf')\\n\",\n", Re);
    fprintf(datei,"    \"plt.xlabel('y-Wert')\\n\",\n");
    fprintf(datei,"    \"plt.xlim((0,1))\\n\",\n");
    fprintf(datei,"    \"plt.show()\"\n");
    fprintf(datei,"   ]\n  }\n ],\n");
    fprintf(datei," \"metadata\": {\n");
    fprintf(datei,"  \"kernelspec\": {\n");
    fprintf(datei,"   \"display_name\": \"Python 3\",\n");
    fprintf(datei,"   \"language\": \"python\",\n");
    fprintf(datei,"   \"name\": \"python3\"\n");
    fprintf(datei,"  },\n");
    fprintf(datei,"  \"language_info\": {\n");
    fprintf(datei,"   \"codemirror_mode\": {\n");
    fprintf(datei,"    \"name\": \"ipython\",\n");
    fprintf(datei,"    \"version\": 3\n");
    fprintf(datei,"   },\n");
    fprintf(datei,"   \"file_extension\": \".py\",\n");
    fprintf(datei,"   \"mimetype\": \"text/x-python\",\n");
    fprintf(datei,"   \"name\": \"python\",\n");
    fprintf(datei,"   \"nbconvert_exporter\": \"python\",\n");
    fprintf(datei,"   \"pygments_lexer\": \"ipython3\",\n");
    fprintf(datei,"   \"version\": \"3.7.3\"\n");
    fprintf(datei,"  }\n");
    fprintf(datei," },\n");
    fprintf(datei," \"nbformat\": 4,\n");
    fprintf(datei," \"nbformat_minor\": 2\n");
    fprintf(datei,"}");


    
    fclose(datei);
}
