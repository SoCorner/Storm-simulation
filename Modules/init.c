#include "init.h"


void READ_PARAMETER(char *Inputfile, REAL* xlength, REAL* ylength, int* imax, int* jmax, REAL* delx, REAL* dely,REAL* delt, REAL* t_end,
                   REAL* tau, int* itermax, REAL* eps, REAL* omg, REAL* alpha, REAL* Re, REAL* GX, REAL* GY, REAL* UI, REAL* VI, REAL* PI){
    FILE* data;
    data = fopen(Inputfile,"r");
    REAL A[18];
    unsigned int i = 0;
    char line[1000];
    if(data == NULL){
        printf("Could not open file");
        
    }
    while(fgets(line,sizeof line,data)!=NULL){
        char* start = line;
        REAL field;
        int n;
        while(sscanf(start,"%lf%n",&field,&n) == 1){
            A[i] = field;
            ++i;
            start += n;
        }
    }
    *xlength = A[0];
    *ylength = A[1];
    *imax = A[2];
    *jmax = A[3];
    *delt = A[4];
    *t_end = A[5];
    *tau = A[6];
    *itermax = A[7];
    *eps = A[8];
    *omg = A[9];
    *alpha = A[10];
    *Re = A[11];
    *GX = A[12];
    *GY = A[13];
    *UI = A[14];
    *VI = A[15];
    *PI = A[16];
    *delx = *xlength/(*imax+2);
    *dely = *ylength/(*jmax+2);
    fclose(data);
}


void INIT_UVP(REAL** U, REAL** V, REAL** P, int imax,int jmax, REAL UI, REAL VI, REAL PI,char*problem){
    fill2Dfield(UI,U,imax+2,jmax+2);
    fill2Dfield(VI,V,imax+2,jmax+2);
    fill2Dfield(PI,P,imax+2,jmax+2);
}

// C_F = 0 for flow
// C_B = 1 for boundary

//===============================================| initialize flag |=====================================
void INIT_FLAG(char*problem, int**FLAG, int imax, int jmax)
{   
    printf("initing flag\n");
    //============ general situation =================
    int C_F=0, C_B=1;
    for(int i=1;i< imax+1;i++){ 
        for(int j=1;j< jmax+1;j++){
             FLAG[i][j]= C_F;
        }
    }
    for(int i=0;i< imax+2;i++){ 
         FLAG[i][0]= C_B;
         FLAG[i][ jmax+1]= C_B;
    }
    for(int j=0;j< jmax+2;j++){ 
         FLAG[0][j]= C_B;
         FLAG[ imax+1][j]= C_B;
    }

    //============ special situation =================

    //============ Karman's situation =================
    if(strcmp(problem,"karman")==0){
        int temp= (jmax)/5;
        printf("temp: %d\n",temp);
        for(int k=0;k<temp;k++){
             FLAG[3*temp-k-1][2*temp+k]= C_B; 
             FLAG[3*temp-k-1][2*temp+k+1]= C_B; 
             FLAG[3*temp-k-1][2*temp+k+2]= C_B;
        }
    }

    //============ Stufe's situation =================
    if(strcmp(  problem,"stairs")==0){
        int temp= jmax/2;
        for(int i=0;i<= temp;i++){
            for (int j=0;j<= temp;j++){
                 FLAG[i][j]= C_B; 
            }
        }
    }
    
    //========================|            generall setting         |================================

	
    int ost = 0b01111;  //15
    int west = 0b10111; //23
    int sued = 0b11011; //27
    int nord = 0b11101; //29
    int sued_ost=0b01011; //11
    int sued_west=0b10011;//19
    int nord_ost=0b01101; //13
    int nord_west=0b10101; //21
   
  
    //========================|            fuer innere Zeile         |================================
    for (int i=1;i< imax+1;i++){
        for(int j=1;j< jmax+1;j++){
            //hindernisse
            if( FLAG[i][j]==1){ 
                // ===================| fluid sued |=======================
                if( FLAG[i+1][j]==0){ 
                    //sued_nord
                    if(FLAG[i-1][j]==0){
                        continue;
                    }
                    if(FLAG[i][j+1]==0){
                        if(FLAG[i][j-1]==0){
                            continue;
                        }
                        FLAG[i][j]=sued_ost;
                        continue;
                    }
                    if(FLAG[i][j-1]==0){
                        FLAG[i][j]=sued_west;
                    }
                    else{
                        FLAG[i][j]=sued;
                    }
                }

                // ===================| fluid nord |=======================
                if( FLAG[i-1][j]==0){ 
                    //sued_nord
                    if(FLAG[i+1][j]==0){
                        continue;
                    }
                    if(FLAG[i][j+1]==0){
                        printf("here u go!\n");
                        if(FLAG[i][j-1]==0){
                            continue;
                        }
                        FLAG[i][j]=nord_ost;
                        continue;
                        //printf("nord_ost: %d\n",nord_ost);
                    }
                    if(FLAG[i][j-1]==0){
                        FLAG[i][j]=nord_west;
                    }
                    else{
                        FLAG[i][j]=nord;
                    }
                }

                // ===================| fluid ost |=======================
                if( FLAG[i][j+1]==0){ 
                    //sued_nord
                    if(FLAG[i][j-1]==0){
                        continue;
                    }
                    if(FLAG[i-1][j]==0){
                        if(FLAG[i+1][j]==0){
                            continue;
                        }
                        FLAG[i][j]=nord_ost;
                        continue;
                        printf("I am here!1\n");
                    }
                    if(FLAG[i+1][j]==0){
                        FLAG[i][j]=sued_ost;
                    }
                    else{
                        FLAG[i][j]=ost;
                        //printf("I am here!\n");
                    }
                }

                // ===================| fluid west |=======================
                if( FLAG[i][j-1]==0){ 
                    //sued_nord
                    if(FLAG[i][j+1]==0){
                        continue;
                    }
                    if(FLAG[i-1][j]==0){
                        if(FLAG[i+1][j]==0){
                            continue;
                        }
                        FLAG[i][j]=nord_west;
                        continue;
                    }
                    if(FLAG[i+1][j]==0){
                        FLAG[i][j]=sued_west;
                    }
                    else{
                        FLAG[i][j]=west;
                    }
                }
            }   
        }
    
    }
    //========================|            fuer Randzeile         |================================
    for(int i=1;i< imax+2;i++){ 
        if( FLAG[i][1]==0){ 
             FLAG[i][0]=  ost;
        }
        if( FLAG[i][ jmax]==0){ 
             FLAG[i][ jmax+1]=  west;
        }
    }
    
    for (int j=1;j< jmax+2;j++){
        if( FLAG[1][j]==0){ 
             FLAG[0][j]= sued;
        }
        if( FLAG[ imax][j]==0){ 
             FLAG[ imax+1][j]=  nord;
        }
    }
}

void READ_PARAMETER_NEW(char *Inputfile, REAL* xlength, REAL* ylength, int* imax, 
                    int* jmax, REAL* delx, REAL* dely,REAL* delt, REAL* t_end,
                   REAL* tau, int* itermax, REAL* eps, REAL* omg, REAL* alpha,
                    REAL* Re, REAL* GX, REAL* GY, REAL* UI, REAL* VI, REAL* PI,int* wl, int*wr, int* wt, int*wb)
{
    FILE* data;
    data = fopen(Inputfile,"r");
    REAL A[21];
    unsigned int i = 0;
    char line[1000];
    if(data == NULL){
        printf("Could not open file");
        
    }
    while(fgets(line,sizeof line,data)!=NULL){
        char* start = line;
        REAL field;
        int n;
        while(sscanf(start,"%lf%n",&field,&n) == 1){
            A[i] = field;
            ++i;
            start += n;
        }
    }
    *xlength = A[0];
    *ylength = A[1];
    *imax = A[2];
    *jmax = A[3];
    *delt = A[4];
    *t_end = A[5];
    *tau = A[6];
    *itermax = A[7];
    *eps = A[8];
    *omg = A[9];
    *alpha = A[10];
    *Re = A[11];
    *GX = A[12];
    *GY = A[13];
    *UI = A[14];
    *VI = A[15];
    *PI = A[16];
    *wl=A[17];
    *wr=A[18];
    *wt=A[19];
    *wb=A[20];
    *delx = *xlength/(*imax+2);
    *dely = *ylength/(*jmax+2);
    fclose(data);
}