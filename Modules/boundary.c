#include "boundary.h"

void SETBCOND(REAL**U, REAL**V, int imax, int jmax){
    for(int j=1;j<jmax+1;++j){
        U[0][j]=0;
        V[0][j]=-V[1][j];
        U[imax][j]=0;
        V[imax+1][j]=-V[imax][j];
    }
    for(int i=1;i<imax+1;++i){
        V[i][0]=0;
        V[i][jmax]=0;
        U[i][0]=-U[i][1];
        U[i][jmax+1]=-U[i][jmax];
    }

}

void SETBCOND_NEW(REAL**U, REAL**V, int**FLAG, int imax, int jmax, int wl, int wr, int wt, int wb)
{   

    //=================================| general bounds based on problem |==========================================
    //left bound
    switch(wl) 
    {
    //Haftbedingunge
    case 1: 
        for(int i=1;i<=imax;i++){
            U[i][0]=-U[i][1];
            V[i][0]=0;
        }
        break;
    //Rutschbedingungen
    case 2: 
        for(int i=1;i<=imax;i++){
            U[i][0]=U[i][1];
             V[i][0]=0;
        }
        break;
    //Ausstroembedingungen
    case 3: 
        for(int i=1;i<= imax;i++){
             U[i][0]= U[i][1];
             V[i][0]= V[i][1];
        }
        break;
    default:
        printf("Please enter a number from {1,2,3}!\n");
        break;
    }
    //right bound
    switch( wr) 
    {
    //Haftbedingungen
    case 1:
        for(int i=1;i<= imax;i++){
             U[i][ jmax+1]=- U[i][ jmax];
             V[i][ jmax]=0;
        }
        break;
    //Rutschbedingungen
    case 2:
        for(int i=1;i<= imax;i++){
             U[i][ jmax+1]= U[i][ jmax];
             V[i][ jmax]=0;
        }
        break;
    //Ausstroembedingungen
    case 3:
        for(int i=1;i<= imax;i++){
             U[i][ jmax+1]= U[i][ jmax];
             V[i][ jmax]= V[i][ jmax-1];
        }
        break;
    default:
        printf("Please enter a number from {1,2,3}!\n");
        break;
    }
    //top bound
    switch( wt) 
    {
    //Haftbedingungen
    case 1:
        for(int j=1;j<= jmax;j++){
             U[0][j]=0;;
             V[0][j]=- V[1][j];
        }
        break;
    //Rutschbedingungen
    case 2:
        for(int j=1;j<= jmax;j++){
             U[0][j]=0;
             V[0][j]= V[1][j];
        }
        break;
    //Ausstroembedingungen
    case 3:
        for(int j=1;j<= jmax;j++){
             U[0][j]= U[1][j];
             V[0][j]= V[1][j];
        }
        break;
    default:
        printf("Please enter a number from {1,2,3}!\n");
        break;
    }
    //bottom border
    switch( wb)
    {
    //Haftbedingungen
    case 1:
        for(int j=1;j<= jmax;j++){
             U[ imax][j]=0;
             V[ imax+1][j]=- V[ imax][j];
        }
        break;
    //Rutschbedingungen
    case 2:
        for(int j=1;j<= jmax;j++){
             U[ imax][j]=0;
             V[ imax+1][j]= V[ imax][j];
        }
        break;
    //Ausstroembedingungen
    case 3:
        for(int j=1;j<= jmax;j++){
             U[ imax][j]= U[ imax-1][j];
             V[ imax+1][j]= V[ imax][j];
        }
        break;
    default:
        printf("Please enter a number from {1,2,3}!\n");
        break;
    }
    
    //=================================| set bounds based on FLag |==========================================
    
   for(int i=0;i< imax+2;i++){
        for(int j=0;j< jmax+2;j++){
            //printf("here is ok 4\n");
            switch( FLAG[i][j])
            {
            case 0b01011:
            //printf("sued_ost: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    if(j==0){
                         U[i][j]=0;
                         V[i][j]=- V[i+1][j];
                    }else{
                        V[i][j-1]=0;
                        U[i][j]=0;
                        V[i][j]=- V[i+1][j];
                    }
                }else{
                    if(j==0){
                        U[i][j]=0;
                        V[i][j]=- V[i+1][j];
                        U[i-1][j]=- U[i-1][j-1];
                    }else{
                        U[i-1][j]=- U[i-1][j-1];
                        V[i][j-1]=0;
                        U[i][j]=0;
                        V[i][j]=- V[i+1][j];
                    }
                     
                }
                break;

            case 0b10011:
            //printf("sued_west: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    if(j==0){
                        break;
                    }else{
                        U[i][j]=- U[i][j-1];
                         V[i][j-1]=0;
                    }
                }else{
                    if(j==0){
                        U[i-1][j]=0;
                        V[i][j]=- V[i-1][j];
                    }else{
                        U[i-1][j]=0;
                        U[i][j]=- U[i][j-1];
                        V[i][j-1]=0;
                        V[i][j]=- V[i-1][j];
                    }
                }
                 
                break;

            case 0b01101: 
            //printf("nord_ost: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    if(j==0){
                        U[i][j]=0;
                        V[i][j]=0;
                    }else{
                        U[i][j]=0;
                        V[i][j]=0;
                        V[i][j-1]=- V[i+1][j-1];
                    }
                }else{
                    if(j==0){
                        U[i][j]=0;
                        V[i][j]=0;
                        U[i-1][j]=- U[i-1][j+1];
                    }else{
                        U[i][j]=0;
                        U[i-1][j]=- U[i-1][j+1];
                        V[i][j]=0;
                        V[i][j-1]=- V[i+1][j-1];
                    }
                }
                 
                break;

            case 0b10101: 
            //printf("nord_west: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    if(j==0){
                        U[i][j]=- U[i][j+1];
                        V[i][j]=0;
                    }else{
                        U[i][j]=- U[i][j+1];
                        V[i][j]=0;
                    }
                }else{
                    if(j==0){
                        U[i-1][j]=0;
                        V[i][j]=0;
                        U[i][j]=- U[i][j+1];
                    }else{
                        U[i-1][j]=0;
                        U[i][j]=- U[i][j+1];
                        V[i][j]=0;
                        V[i][j-1]=- V[i-1][j-1];
                    }
                }
                break;
            case 0b01111: 
                //printf("osten: FLAG[%d][%d].\n",i,j);
                if(j==0){
                    U[i][j]=0;
                    V[i][j]=- V[i+1][j];
                }else{
                    U[i][j]=0;
                    V[i][j]=- V[i+1][j];
                    V[i][j-1]=- V[i+1][j-1];
                }
                 
                break;

            case 0b10111: 
               // printf("west: FLAG[%d][%d].\n",i,j);
               if(i==0){
                   break;
               }else{
                   if(j==0){
                       U[i-1][j]=0;
                       V[i][j]=- V[i-1][j];
                   }else{
                       U[i-1][j]=0;
                        V[i][j-1]=- V[i-1][j-1];
                        V[i][j]=- V[i-1][j];
                   }
               }
                 
                break;

            case 0b11011: 
                //printf("sued: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    if(j==0){
                        break;
                    }else{
                        U[i][j]=- U[i][j-1];
                        V[i][j-1]=0;
                    }
                }else{
                    if(j==0){
                        break;
                    }else{
                        U[i-1][j]=- U[i-1][j-1];
                        U[i][j]=- U[i][j-1];
                        V[i][j-1]=0;
                    }
                }
                 
                break;
            
            case 0b11101: 
                //printf("nord: FLAG[%d][%d].\n",i,j);
                if(i==0){
                    U[i][j]=- U[i][j+1];
                    V[i][j]=0;
                }else{
                    U[i-1][j]=- U[i-1][j+1];
                    U[i][j]=- U[i][j+1];
                    V[i][j]=0;
                }
                 
                break;
            
            }
        }
    }
    
}

//================================|set bound for special problem|===============================


void SETSPECBOND(REAL** U,  REAL**V, int imax, int jmax, char* problem)
{   
    //printf("now we are dcavity\n");
    if(strcmp( problem,"dcavity")==0){ 
    //printf("now we are dcavity.");
        for(int i=1;i<= imax;i++){
             U[i][jmax+1]=(REAL)2- U[i][ jmax]; 
        }
    }
    if(strcmp( problem,"karman")==0){
        for(int j=1;j<= jmax;j++){
             U[0][j]=1; 
        }  
		for(int j = 1 ; j < jmax + 1 ; j++ ){
        U[imax][j] = U[imax - 1][j] ;
		//U[imax][j]=U[imax-1][j];
		}    
		//bot no slip
		for(int i = 1 ; i < imax +1 ; i++ ){
			U[i][0] =  - U[i][1];
		}
		// top no slip
		for(int i = 1 ; i < imax +1 ; i++ ){
			U[i][jmax+1] =  - U[i][jmax];
		}
		 for(int j = 1 ; j < jmax + 1 ; j++ ){
            V[0][j] =  V[1][j];
		}
		// right bound: outflow
		for(int j = 1 ; j < jmax + 1 ; j++ ){
				V[imax + 1][j] =  V[imax][j];
		}
		// top bound: no slip: 
		for(int i = 1 ; i < imax +1 ; i++ ){
				V[i][jmax] = 0; 
		}
		// bot bound: no slip: 
		for(int i = 1 ; i < imax +1 ; i++ ){
				V[i][0] = 0;
		}
    }
    if(strcmp( problem,"stairs")==0){
		int temp= jmax/2;
        for(int i=0;i<= temp;i++){
            for (int j=0;j<= temp;j++){
                 U[i][j]= 0; 
				 V[i][j]= 0; 
            }
        }
        for(int j=1;j<= jmax;j++){
             U[0][j]=1; //geschwindigkeit auf 1 gesetzt
        }  
		for(int j = 1 ; j < jmax + 1 ; j++ ){
        U[imax][j] = U[imax - 1][j] ;
		//U[imax][j]=U[imax-1][j];
		}    
		//bot no slip
		for(int i = 1 ; i < imax +1 ; i++ ){
			U[i][0] =  - U[i][1];
		}
		// top no slip
		for(int i = 1 ; i < imax +1 ; i++ ){
			U[i][jmax+1] =  - U[i][jmax];
		}
		 for(int j = 1 ; j < jmax + 1 ; j++ ){
            V[0][j] =  V[1][j];
		}
		// right bound: outflow
		for(int j = 1 ; j < jmax + 1 ; j++ ){
				V[imax + 1][j] =  V[imax][j];
		}
		// top bound: no slip: 
		for(int i = 1 ; i < imax +1 ; i++ ){
				V[i][jmax] = 0; 
		}
		// bot bound: no slip: 
		for(int i = 1 ; i < imax +1 ; i++ ){
				V[i][0] = 0;
		}

	}
}

