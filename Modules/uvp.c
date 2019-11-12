#include "uvp.h"

REAL min2(REAL x, REAL y){
    if(x<y){
        return x;
    }
    else return y;
}
REAL min3(REAL x, REAL y, REAL z){
    if(x<y){
        return min2(x,z);
    }
    else return min2(y,z);
}


void COMP_delt(REAL* delt,int imax, int jmax, REAL delx, REAL dely, REAL** U, REAL** V, REAL Re, REAL tau){
    REAL umax = 0;
    REAL vmax = 0;
    //double temp = Re/2 * 1.0/(1.0/(delx*delx) + 1.0/(dely*dely));
    for(int i=0;i<imax+2;++i){
        for(int j=0;j<jmax+2;++j){
            if(fabs(U[i][j])>umax){
                umax = fabs(U[i][j]);
            }
            if(fabs(V[i][i])>vmax){
                vmax = fabs(V[i][j]);
            }
        }
    }
    if(tau>=0){
        //printf("the value of umax is : %lf, and vmax is : %lf\n",umax, dely/umax);
        *delt = tau * min3(Re/2 * 1.0/(1.0/(delx*delx) + 1.0/(dely*dely)), delx/fabs(umax),dely/fabs(vmax));
        //*delt=Re/2 * 1.0/(1.0/(delx*delx)+ 1.0/(dely*dely));
    }
}
void COMP_delt_NEW(REAL* delt,int imax, int jmax, REAL delx, REAL dely, REAL** U, REAL** V, REAL Re, REAL tau){
{
    if( tau>=0){
        REAL umax=0;
        REAL vmax=0;
        REAL delx2inv=(double)1/pow( delx,2);
        REAL dely2inv=(double)1/pow( dely,2);
        for(int i=1;i< imax+1;i++){
            for(int j=1;j< jmax+1;j++){
                if(fabs( U[i][j])>umax){
                        umax=fabs( U[i][j]);
                }
                if(fabs( V[i][j])>vmax){
                        vmax=fabs( V[i][j]);
                }
            }
        }
        if(umax>0){
            if(vmax>0){
                if(( Re/(double)2)*(double)1/(delx2inv+dely2inv)< delx/umax && ( Re/(double)2)*(double)1/(delx2inv+dely2inv)< dely/vmax){
                     *delt= tau*( Re/(double)2)*(double)1/(delx2inv+dely2inv);  //re<dx,dy
                }else{
                    if( delx/umax< dely/vmax){
                         *delt= tau* delx/umax; //dx<dy,re
                    }else{
                         *delt= tau* dely/vmax; //dy<re,dx
                        }
                    }
            }else{                                                                      //vmax=0
                if(( Re/(double)2)*(double)1/(delx2inv+dely2inv)< delx/umax){
                     *delt= tau*( Re/(double)2)*(double)1/(delx2inv+dely2inv);  //re<dx
                }else{
                     *delt= tau* delx/umax;                                //dx<re
                }
            }
        }else{                                                                          //umax=0
            if(vmax>0){
                if(( Re/(double)2)*(double)1/(delx2inv+dely2inv)< dely/vmax){
                     *delt=  tau*( Re/(double)2)*(double)1/(delx2inv+dely2inv); //re<dy
                }else{
                     *delt= tau* dely/vmax; //dy<re
                    }
            }else{
                 *delt=  tau*( Re/(double)2)*(double)1/(delx2inv+dely2inv); //nur re
            }
        }
    }
}
}
void COMP_FG(REAL**U, REAL**V, REAL**F, REAL**G, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re){
    /*double f1,f2,f3,f4,f5;
    double g1,g2,g3,g4,g5;*/
    for(int i=1;i<imax;++i){
        for(int j=1;j<jmax+1;++j){
            /*f1 = U[i][j] + delt * 1/Re * (1/(delx*delx) * (U[i+1][j] - 2*U[i][j] + U[i-1][j]) + 1/(dely*dely) * (U[i][j+1]-2*U[i][j] + U[i][j-1]));
            f2 =  - delt * 1/delx * (pow(0.5 * (U[i][j] + U[i+1][j]),2.0) - pow(0.5 * (U[i-1][j] + U[i][j]),2.0));
            f3 =  - delt * alpha/delx *(0.25 * fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) - 0.25 * fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]));
            f4 =  - delt * 1/dely * ((V[i][j] + V[i+1][j])/2 * (U[i][j] + U[i][j+1])/2 - (V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] + U[i][j])/2);
            f5 =  - delt * alpha/dely *(fabs(V[i][j] + V[i+1][j])/2 * (U[i][j] - U[i][j+1])/2 - fabs(V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] - U[i][j])/2);
            F[i][j] = f1+f2+f3+f4+f5 + delt * GX;*/
            REAL dudx2=(U[i+1][j]-2*U[i][j]+U[i-1][j])/pow(delx,2);
            REAL dudy2=(U[i][j+1]-2*U[i][j]+U[i][j-1])/pow(dely,2);
            //du^2/dx
            REAL du2dx=1/delx*(pow((U[i][j]+U[i+1][j])/2,2)-pow((U[i-1][j]+U[i][j])/2,2))+
                alpha/delx*(fabs(U[i][j]+U[i+1][j])/2*(U[i][j]-U[i+1][j])/2-fabs(U[i-1][j]+U[i][j])/2*(U[i-1][j]-U[i][j])/2);
            REAL duvdy=1/dely*((V[i][j]+V[i+1][j])/2*(U[i][j]+U[i][j+1])/2-(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]+U[i][j])/2)+alpha/delx*(fabs(V[i][j]+V[i+1][j])/2*(U[i][j]-U[i][j+1])/2-fabs(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]-U[i][j])/2);
            F[i][j]=U[i][j]+delt*((double)1/Re*(dudx2+dudy2)-du2dx-duvdy+GX);
        }
    }
    for(int i=1;i<imax+1;++i){
        for(int j=1;j<jmax;++j){
            /*g1 = V[i][j] + delt * 1/Re * (1/(delx*delx) * (V[i+1][j] - 2*V[i][j] + V[i-1][j]) + 1/(dely*dely) * (V[i][j+1]-2*V[i][j] + V[i][j-1]));
            g2 =  - delt * 1/dely * (pow(0.5 * (V[i][j] + V[i][j+1]),2) - pow(1/2 * (V[i][j-1] + V[i][j]),2));
            g3 =  - delt * alpha/dely *(0.25 * fabs(V[i][j] + V[i][j+1]) * (V[i][j] - V[i][j+1]) - 1/4 * fabs(V[i][j-1] + V[i][j]) * (V[i][j-1] - V[i][j]));
            g4 =  - delt * 1/delx * ((U[i][j] + U[i][j+1])/2 * (V[i][j] + V[i+1][j])/2 - (U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] + V[i][j])/2);
            g5 =  - delt * alpha/delx *(fabs(U[i][j] + U[i][j+1])/2 * (V[i][j] - V[i+1][j])/2 - fabs(U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] - V[i][j])/2);
            G[i][j] = g1+g2+g3+g4+g5 + delt * GY;*/
            REAL dvdx2=(V[i+1][j]-2*V[i][j]+V[i-1][j])/pow(delx,2);
            REAL dvdy2=(V[i][j+1]-2*V[i][j]+V[i][j-1])/pow(dely,2);
            REAL duvdx=1/delx*((U[i][j]+U[i][j+1])/2*(V[i][j]+V[i+1][j])/2-(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]+V[i][j])/2)+alpha/delx*(fabs(U[i][j]+U[i][j+1])/2*(V[i][j]-V[i+1][j])/2-fabs(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]-V[i][j])/2);
            REAL dv2dy=1/dely*(pow((V[i][j]+V[i][j+1])/2,2)-pow((V[i][j-1]+V[i][j])/2,2))+alpha/dely*(fabs(V[i][j]+V[i][j+1])/2*(V[i][j]-V[i][j+1])/2-fabs(V[i][j-1]+V[i][j])/2*(V[i][j-1]-V[i][j])/2);
            G[i][j]=V[i][j]+delt*((double)1/Re*(dvdx2+dvdy2)-duvdx-dv2dy+GY);
        }
    }
    for(int j=1;j<jmax+1;++j){
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
    }
    for(int i=0;i<imax+1;++i){
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];
    }
}

void COMP_RHS(REAL**F, REAL**G, REAL**RHS, int imax, int jmax, REAL delt, REAL delx, REAL dely)
{   
    //printf("ready\n");
    for(int i=1;i<imax+1;++i){
        for(int j=1;j<jmax+1;++j){
            //printf("%d_%d\n",i,j);
            RHS[i][j] = 1/delt * ((F[i][j] - F[i-1][j])/delx + (G[i][j] - G[i][j-1])/dely);
        }
    }
    //printf("done\n");
}

int POISSON(REAL**P, REAL**RHS, REAL**U, REAL**V, int imax, int jmax, REAL delx, REAL dely, REAL eps, int itermax, REAL omg, REAL *res){

    for(int k=0;k<itermax;++k){
        for(int j=1;j<jmax+1;++j){
            P[0][j] = P[1][j];
            P[imax+1][j] = P[imax][j];
        }
        for(int i=1;i<imax+1;++i){
            P[i][0] = P[i][1];
            P[i][jmax+1] = P[i][jmax];
        }
        for(int i=1;i<imax+1;++i){
            for(int j=1;j<jmax+1;++j){
                P[i][j] = (1-omg) * P[i][j] + omg/(2*(1/(delx*delx) + 1/(dely*dely))) * (1/(delx*delx) * (P[i+1][j] + P[i-1][j]) + 1/(dely*dely) * (P[i][j+1] + P[i][j-1]) - RHS[i][j]);
            }
        }
        //*res=0.0;
        for(int i=1;i<imax+1;++i){
            for(int j=1;j<jmax+1;++j){
                *res = 1.0/(imax*jmax) * pow((P[i+1][j] - 2*P[i][j] + P[i-1][j])/(delx*delx) + (P[i][j+1] -2*P[i][j] + P[i][j-1])/(dely*dely) - RHS[i][j],2);
            }
        }
        *res=sqrt(*res);
        if(*res<eps){
            return k;
        }
    }
    return itermax;
}

void ADAP_UV(REAL**U, REAL**V, REAL**F, REAL**G, REAL**P, int imax, int jmax, REAL delt, REAL delx, REAL dely){
    for(int i=1;i<imax;++i){
        for(int j=1;j<jmax+1;++j){
            U[i][j] = F[i][j] - delt/delx * (P[i+1][j] - P[i][j]);
        }
    }
    for(int i=1;i<imax+1;++i){
        for(int j=1;j<jmax;++j){
            V[i][j] = G[i][j] - delt/dely * (P[i][j+1] - P[i][j]);
        }
    }
}
//=====================================| new programms |==================================



//=====================================| new poisson |==================================

int POISSON_NEW(REAL**P, REAL**RHS, REAL**U, REAL**V, int** FLAG, int imax, int jmax, int flowNr, REAL delx, REAL dely, REAL eps, int itermax, REAL omg, REAL *res){

    //=====================================| count the number of fluid |==================================

    //caculate the number of fluids
    flowNr=0;
    for(int i=0;i< imax+2;i++){
        for(int j=0;j< jmax+2;j++){
            if(FLAG[i][j]==0){
                flowNr++;
            }
        }
    }

    for(int k=0;k<itermax;++k){

        //=====================================| initalize the press based on FLAG |=============================
        /*for(int j=1;j<jmax+1;++j){
            P[0][j] = P[1][j];
            P[imax+1][j] = P[imax][j];
        }
        for(int i=1;i<imax+1;++i){
            P[i][0] = P[i][1];
            P[i][jmax+1] = P[i][jmax];
        }*/
        *res=0;
        applyHomogeneousNeumannBC(P,imax,jmax);
        COMP_PRESS(P, FLAG, imax, jmax, delx, dely);
        for(int i=1;i<imax+1;++i){
            for(int j=1;j<jmax+1;++j){
                // only for flow
                if(FLAG[i][j]==0){
                    *res += 1.0/(flowNr) * pow((P[i+1][j] - 2*P[i][j] + P[i-1][j])/(delx*delx) + (P[i][j+1] -2*P[i][j] + P[i][j-1])/(dely*dely) - RHS[i][j],2);
                }
            }
        }
        *res=sqrt(*res);
        if(*res<eps){
            return k;
        }
       

        //=====================================| apply SOR only for flow |=============================
        for(int i=1;i<imax+1;++i){
            for(int j=1;j<jmax+1;++j){
                if(FLAG[i][j]==0){
                    P[i][j] = (1-omg) * P[i][j] + omg/(2*(1/(delx*delx) + 1/(dely*dely))) * (1/(delx*delx) * (P[i+1][j] + P[i-1][j]) + 1/(dely*dely) * (P[i][j+1] + P[i][j-1]) - RHS[i][j]);
                }
            }
        }
        
        //=====================================| compute the residum |=============================
        //*res=0.0;
        

        //=====================================| test for the tolerance |=============================
        
    }
    return itermax;
}

//========================================| for F and G |=====================================================
void COMP_FG_NEW(REAL**U, REAL**V, REAL**F, REAL**G, int** FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re)
{
    for (int i=1;i< imax;i++){
        for (int j=1;j< jmax+1;j++){
            if( FLAG[i][j]==0 &&  FLAG[i+1][j]==0){ //fluidzellen
            REAL dudx2=(U[i+1][j]-2*U[i][j]+U[i-1][j])/pow(delx,2);
            REAL dudy2=(U[i][j+1]-2*U[i][j]+U[i][j-1])/pow(dely,2);
            REAL du2dx=1/delx*(pow((U[i][j]+U[i+1][j])/2,2)-pow((U[i-1][j]+U[i][j])/2,2))+
                alpha/delx*(fabs(U[i][j]+U[i+1][j])/2*(U[i][j]-U[i+1][j])/2-fabs(U[i-1][j]+U[i][j])/2*(U[i-1][j]-U[i][j])/2);
            REAL duvdy=1/dely*((V[i][j]+V[i+1][j])/2*(U[i][j]+U[i][j+1])/2-(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]+U[i][j])/2)+alpha/delx*(fabs(V[i][j]+V[i+1][j])/2*(U[i][j]-U[i][j+1])/2-fabs(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]-U[i][j])/2);
            F[i][j]=U[i][j]+delt*((double)1/Re*(dudx2+dudy2)-du2dx-duvdy+GX);
            }
        }
    }
    for(int i=1;i< imax+1;i++){
        for(int j=1;j< jmax;j++){
            if( FLAG[i][j]==0 &&  FLAG[i][j+1]==0){ //fluidzellen
                REAL dvdx2=(V[i+1][j]-2*V[i][j]+V[i-1][j])/pow(delx,2);
                REAL dvdy2=(V[i][j+1]-2*V[i][j]+V[i][j-1])/pow(dely,2);
                REAL duvdx=1/delx*((U[i][j]+U[i][j+1])/2*(V[i][j]+V[i+1][j])/2-(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]+V[i][j])/2)+alpha/delx*(fabs(U[i][j]+U[i][j+1])/2*(V[i][j]-V[i+1][j])/2-fabs(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]-V[i][j])/2);
                REAL dv2dy=1/dely*(pow((V[i][j]+V[i][j+1])/2,2)-pow((V[i][j-1]+V[i][j])/2,2))+alpha/dely*(fabs(V[i][j]+V[i][j+1])/2*(V[i][j]-V[i][j+1])/2-fabs(V[i][j-1]+V[i][j])/2*(V[i][j-1]-V[i][j])/2);
                G[i][j]=V[i][j]+delt*((double)1/Re*(dvdx2+dvdy2)-duvdx-dv2dy+GY);
            }
        }
    }
    for(int j=1;j< jmax+1;j++){ //randwerte fuer F
         F[0][j]= U[0][j];
         F[ imax][j]= U[ imax][j];
    }
    for(int i=1;i< imax+2;i++){//randwerte fuer G
         G[i][0]= V[i][0];
         G[i][ jmax]= V[i][ jmax];
    }
    for(int i=0;i< imax+2;i++){ //F,G Randwerte fuer Randzellen
        for(int j=0;j< jmax+2;j++){
            switch( FLAG[i][j])
            {
            case 0b01011://sued, ost
                 F[i][j]= U[i][j];
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b10011://sued, west
                 F[i-1][j]= U[i-1][j];
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b01101: //nord, ost
                 F[i][j]= U[i][j];
                 G[i][j]= V[i][j];
                break;

            case 0b10101: //nord, west
                 F[i-1][j]= U[i-1][j];
                 G[i][j]= V[i][j];
                break;

            case 0b01111: //osten
                 F[i][j]= U[i][j];
                break;

            case 0b10111: //west
                 F[i-1][j]= U[i-1][j];
                break;

            case 0b11011: //sued
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b11101: //nord
                 G[i][j]= V[i][j];
                break;

            }
        }
    }
}

//==========================| output of U and V |==============================================

void ADAP_UV_NEW(REAL**U, REAL**V, REAL**F, REAL**G, REAL**P, int**FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely)
{
    for(int i=1;i< imax;i++){
        for(int j=1;j< jmax+1;j++){
            if( FLAG[i][j]==0 &&  FLAG[i+1][j]==0){//fluidfelder
                 U[i][j]= F[i][j]- delt/ delx*( P[i+1][j]- P[i][j]);
            }
        }
    }
    for(int i=1;i< imax+1;i++){
        for(int j=1;j< jmax;j++){
            if( FLAG[i][j]==0 &&  FLAG[i][j+1]==0){//fluidfelder
                 V[i][j]= G[i][j]- delt/ dely*( P[i][j+1]- P[i][j]);
            }
        }
    }
    
}

//==========================| set the value of P based on FLAG |==============================================

void COMP_PRESS(REAL** P, int** FLAG, int imax, int jmax, REAL delx, REAL dely)
{
    REAL dx2=pow( delx,2);
    REAL dy2=pow( dely,2);
    for(int i=0;i< imax+2;i++){
        for(int j=0;j< jmax+2;j++){
            switch( FLAG[i][j])
            {
            case 0b01011://sued, ost
                 P[i][j]=1/(dx2+dy2)*(dx2* P[i][j+1]+dy2* P[i-1][j]);
                break;

            case 0b10011://sued, west
                 P[i][j]=1/(dx2+dy2)*(dx2* P[i][j-1]+dy2* P[i-1][j]);
                break;

            case 0b01101: //nord, ost
                 P[i][j]=1/(dx2+dy2)*(dx2* P[i][j+1]+dy2* P[i+1][j]);
                break;

            case 0b10101: //nord, west
                 P[i][j]=1/(dx2+dy2)*(dx2* P[i][j+1]+dy2* P[i-1][j]);
                break;

            case 0b01111: //osten
                 P[i][j]= P[i+1][j];
                break;

            case 0b10111: //west
                 P[i][j]= P[i-1][j];
                break;

            case 0b11011: //sued
                 P[i][j]= P[i][j-1];
                break;

            case 0b11101: //nord
                 P[i][j]= P[i][j+1];
                break;
            }
        }
    }
}

void COMP_FG_NEW1(REAL**U, REAL**V, REAL**F, REAL**G, int** FLAG, int imax, int jmax, REAL delt, REAL delx, REAL dely, REAL GX, REAL GY, REAL alpha, REAL Re)
{
     F=create2Dfield( imax+2, jmax+2);
     G=create2Dfield( imax+2, jmax+2);
    for (int i=1;i< imax;i++){
        for (int j=1;j< jmax+1;j++){
            if( FLAG[i][j]==0 &&  FLAG[i+1][j]==0){ //fluidzellen
                REAL dudx2=( U[i+1][j]-2* U[i][j]+ U[i-1][j])/pow( delx,2);
                REAL dudy2=( U[i][j+1]-2* U[i][j]+ U[i][j-1])/pow( dely,2);
                REAL du2dx=1/ delx*(pow(( U[i][j]+ U[i+1][j])/2,2)-pow(( U[i-1][j]+ U[i][j])/2,2))+ alpha/ delx*(fabs( U[i][j]+ U[i+1][j])/2*( U[i][j]- U[i+1][j])/2-fabs( U[i-1][j]+ U[i][j])/2*( U[i-1][j]- U[i][j])/2);
                REAL duvdy=1/ dely*(( V[i][j]+ V[i+1][j])/2*( U[i][j]+ U[i][j+1])/2-( V[i][j-1]+ V[i+1][j-1])/2*( U[i][j-1]+ U[i][j])/2)+ alpha/ delx*(fabs( V[i][j]+ V[i+1][j])/2*( U[i][j]- U[i][j+1])/2-fabs( V[i][j-1]+ V[i+1][j-1])/2*( U[i][j-1]- U[i][j])/2);
                 F[i][j]= U[i][j]+ delt*((double)1/ Re*(dudx2+dudy2)-du2dx-duvdy+ GX);
            }
        }
    }
    for(int i=1;i< imax+1;i++){
        for(int j=1;j< jmax;j++){
            if( FLAG[i][j]==0 &&  FLAG[i][j+1]==0){ //fluidzellen
                REAL dvdx2=( V[i+1][j]-2* V[i][j]+ V[i-1][j])/pow( delx,2);
                REAL dvdy2=( V[i][j+2]-2* V[i][j]+ V[i][j-1])/pow( dely,2);
                REAL duvdx=1/ delx*(( U[i][j]+ U[i][j+1])/2*( V[i][j]+ V[i+1][j])/2-( U[i-1][j]+ U[i-1][j+1])/2*( V[i-1][j]+ V[i][j])/2)+ alpha/ delx*(fabs( U[i][j]+ U[i][j+1])/2*( V[i][j]- V[i+1][j])/2-fabs( U[i-1][j]+ U[i-1][j+1])/2*( V[i-1][j]- V[i][j])/2);
                REAL dv2dy=1/ dely*(pow(( V[i][j]+ V[i][j+1])/2,2)-pow(( V[i][j-1]+ V[i][j])/2,2))+ alpha/ dely*(fabs( V[i][j]+ V[i][j+1])/2*( V[i][j]- V[i][j+1])/2-fabs( V[i][j-1]+ V[i][j])/2*( V[i][j-1]- V[i][j])/2);
                 G[i][j]= V[i][j]+ delt*((double)1/ Re*(dvdx2+dvdy2)-duvdx-dv2dy+ GY);
            }
        }
    }
    for(int j=1;j< jmax+1;j++){ //randwerte fuer F
         F[0][j]= U[0][j];
         F[ imax][j]= U[ imax][j];
    }
    for(int i=1;i< imax+2;i++){//randwerte fuer G
         G[i][0]= V[i][0];
         G[i][ jmax]= V[i][ jmax];
    }
    for(int i=0;i< imax+2;i++){ //F,G Randwerte fuer Randzellen
        for(int j=0;j< jmax+2;j++){
            switch( FLAG[i][j])
            {
            case 0b10101://sued, ost
                 F[i][j]= U[i][j];
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b01101://sued, west
                 F[i-1][j]= U[i-1][j];
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b10011: //nord, ost
                 F[i][j]= U[i][j];
                 G[i][j]= V[i][j];
                break;

            case 0b01011: //nord, west
                 F[i-1][j]= U[i-1][j];
                 G[i][j]= V[i][j];
                break;

            case 0b10001: //osten
                 F[i][j]= U[i][j];
                break;

            case 0b01001: //west
                 F[i-1][j]= U[i-1][j];
                break;

            case 0b00101: //sued
                 G[i][j-1]= V[i][j-1];
                break;

            case 0b00011: //nord
                 G[i][j]= V[i][j];
                break;

            }
        }
    }
}

