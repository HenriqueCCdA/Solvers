#include<Solv.h>
/**********************************************************************
 * PCG  : metodo do gradiente conjugado com precondiconador diagonal  *
 * (M-1Ax=M-1b) (matriz simentrica)                                   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   z  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 *maxIt -> numero maximo de iteracoes                                 *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
*********************************************************************/
void pcg(INT const neq      ,INT const nad  
        ,INT *restrict ia   ,INT *restrict ja
        ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
        ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
        ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
        ,unsigned int maxIt ,bool newX          
        ,FILE* fLog         ,bool log
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  int k=0;
  DOUBLE alpha,beta,d,conv,energy;
  DOUBLE timei,timef;
  timei = getTimeC();

/* chute inicial*/
  if(newX)  
    for(i = 0; i < neq; i++)  
      x[i] = 0.e0;
      
  matvec(neq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < neq; i++)   {
    r[i] = b[i] - z[i];
    z[i] = r[i] / m[i];
    b[i] = z[i];
  }
  d    = dot(r,z,neq);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(neq,ia,ja,al,ad,b,z);
    alpha = d / dot(b,z,neq);

    for(i = 0; i < neq; i++)   {
      x[i] +=  alpha * b[i];
      r[i] -=  alpha * z[i];
      z[i] = r[i] / m[i];
    }

    beta = dot(r,z,neq)/d;

    for(i = 0; i < neq; i++)   {
      b[i] = z[i] + beta * b[i];
    }
    d = beta * d;
    if (sqrt(fabs(d)) < conv) break;
    if( k == 1000){ 
      printf("it: %d %e %e\n",j,sqrt(fabs(d)),conv);
      k = 0; 
    }
    k++;
  }
/* -------------------------------------------------------*/
  matvec(neq,ia,ja,al,ad,x,z);
/*norma de energia = xTAx* */
  energy = dot(x,z,neq);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;   
  printf("\tnad         :      %20d\n"  ,nad);
  printf("\tSolver conv :      %20.2e\n",conv);
  printf("\tSolver tol  :      %20.2e\n",tol);
  printf(" (PCG) solver:\n"
         "\tEquations   =      %20d\n"
         "\tIterarions  =      %20d\n"
	 "\tEnergy norm =      %20.12e\n"
	 "\tCPU time(s) =      %20.5lf\n" 
	 ,neq,j+1,energy,timef);
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/



/**********************************************************************/
void ic0cg(INT const neq      ,INT const nad  
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict mIcd 
          ,DOUBLE *restrict mIcCsr
          ,DOUBLE *restrict mIcCsc
          ,INT *restrict iaMCsr,INT *restrict jaMCsr
          ,INT *restrict iaMCsc,INT *restrict jaMCsc
          ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
          ,unsigned int maxIt ,bool newX          
          ,FILE* fLog         ,bool log
          ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  int k=0;
  DOUBLE alpha,beta,d,conv,energy;
  DOUBLE timei,timef;
  timei = getTimeC();

/* chute inicial*/
  if(newX)  
    for(i = 0; i < neq; i++)  
      x[i] = 0.e0;
      
  matvec(neq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < neq; i++)
    r[i] = b[i] - z[i];
  
/*Mz=r*/
  ic0Solv(mIcd  
         ,mIcCsr
         ,mIcCsc  
         ,z
         ,r
         ,iaMCsr  
         ,jaMCsr  
         ,iaMCsc  
         ,jaMCsc  
         ,neq
         ,nad);
  
  for(i = 0; i < neq; i++)   
    b[i] = z[i];


  d    = dot(r,z,neq);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(neq,ia,ja,al,ad,b,z);
    alpha = d / dot(b,z,neq);

    for(i = 0; i < neq; i++)   {
      x[i] +=  alpha * b[i];
      r[i] -=  alpha * z[i];
    }

/*... Mz=r*/
    ic0Solv(mIcd  
           ,mIcCsr
           ,mIcCsc  
           ,z
           ,r
           ,iaMCsr  
           ,jaMCsr  
           ,iaMCsc  
           ,jaMCsc  
           ,neq
           ,nad);
/*.........................................................*/
    beta = dot(r,z,neq)/d;

    for(i = 0; i < neq; i++)   {
      b[i] = z[i] + beta * b[i];
    }
    d = beta * d;
    if (sqrt(fabs(d)) < conv) break;
    if( k == 1000){ 
      printf("it: %d %e %e\n",j,sqrt(fabs(d)),conv);
      k = 0; 
    }
    k++;
  }
/* -------------------------------------------------------*/
  matvec(neq,ia,ja,al,ad,x,z);
/*norma de energia = xTAx* */
  energy = dot(x,z,neq);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;   
  printf("\tnad         :      %20d\n"  ,nad);
  printf("\tSolver conv :      %20.2e\n",conv);
  printf("\tSolver tol  :      %20.2e\n",tol);
  printf(" (ICCG) solver:\n"
         "\tEquations   =      %20d\n"
         "\tIterarions  =      %20d\n"
	 "\tEnergy norm =      %20.12e\n"
	 "\tCPU time(s) =      %20.5lf\n" 
	 ,neq,j+1,energy,timef);
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/


/**********************************************************************/
  
void ic0Solv(DOUBLE *restrict md 
            ,DOUBLE *restrict mCsrUp
            ,DOUBLE *restrict mCscUp
            ,DOUBLE *restrict x
            ,DOUBLE *restrict y
            ,INT *restrict iaCsrUp
            ,INT *restrict jaCsrUp
            ,INT *restrict iaCscUp
            ,INT *restrict jaCscUp
            ,INT const neq
            ,INT const nad){

  int ia1,ia2;
  int jak1;
  int i,j;  

  
  for(i=0;i<neq;i++)
    x[i] = y[i];

/*... triangular inferior*/
  x[0] /= md[0];
  for(i=1;i<neq;i++){
    ia1 = iaCscUp[i];
    ia2 = iaCscUp[i+1];
    for(j=ia1;j<ia2;j++){
      jak1  = jaCscUp[j];
      x[i] -= mCscUp[j]*x[jak1]; 
    }
    x[i] /= md[i];
  }      

/*... triangular superior*/
  x[neq-1] /= md[neq-1];
  for(i=neq-2;i>-1;i--){
     ia1 = iaCsrUp[i];
     ia2 = iaCsrUp[i+1];
     for(j=ia1;j<ia2;j++){
       jak1  = jaCsrUp[j];
       x[i] -= mCsrUp[j]*x[jak1]; 
     }
    x[i] /= md[i];
  }      


}  
/**********************************************************************/



/**********************************************************************/
void icZeroCsr(DOUBLE *restrict m, DOUBLE *restrict ad
              ,DOUBLE *restrict aUp
              ,DOUBLE *restrict mLo
              ,DOUBLE *restrict w   
              ,INT *restrict iaUp,INT *restrict jaUp
              ,INT *restrict iaLo,INT *restrict jaLo
              ,INT    *restrict ws1
              ,INT    *restrict ws2
              ,INT    *restrict ws3
              ,INT const neq,INT const nad){
  
  int iaUp1,iaUp2;
  int iakk;
  int iaUp3,iaUp4;
  int jak1,jak2;
  int i,j,k,ii;
  double tmp = 0.0;
  double timei=0.0; 
//for(k = 0; k < neq+1;k++){
//  printf("%d %d %d\n",k+1,iaUp[k],iaLo[k]);
//}
  
//for(k = 0; k < nad;k++){
//  printf("%d %d %d %lf\n",k+1,jaLo[k],jaUp[k],aUp[k]);
//}
  
  for(k = 0; k < neq;k++){
     m[nad+k]  = ad[k];   
     ws1[k  ]  = -1;
     ws2[k  ]  = -1;
     ws3[k  ]  =  0;
     w  [k  ]  =  0.0;
  }
  
  for(k = 0; k < nad;k++){
     m[k]   =  aUp[k];   
  }

  for(j = 0; j< neq  ;j++){
    iaUp1 = iaUp[j  ];
    iaUp2 = iaUp[j+1];
    
    m[nad+j] = sqrt(m[nad+j]);
      
    for(i=iaUp1,ii=0;i<iaUp2;i++,ii++){
       jak1        = jaUp[i];
         w[jak1]   = m[i];
    }
//    printf("%lf %lf %lf %lf %lf\n",w[0],w[1],w[2],w[3],w[4]);
    
    for(i=iaLo[j];i<iaLo[j+1];i++){
       jak1        = jaLo[i];
       ws1[jak1]   = i;
//       printf("%d %d %d\n",j,jak1,i);
    }
//  printf("%d %d %d %d %d\n",ws1[0],ws1[1],ws1[2],ws1[3],ws1[4]);
    
    timei = getTimeC() - timei;   
    for(k=0;k<j;k++){
      if(ws1[k] != -1){ 
        iaUp3        = iaUp[k  ];
        iaUp4        = iaUp[k+1];
        for(i=iaUp3;i<iaUp4;i++){
          jak1         = jaUp[i];
          ws2[jak1]    = i;
        }
//      printf("%d %d %d %d %d\n",ws2[0],ws2[1],ws2[2],ws2[3],ws2[4]);
        tmp          =   mLo[iaLo[j]+k];
//    printf("%d %d %d %lf\n",j+1,iaUp3,iaUp4,tmp);
        for(i=j+1;i<neq;i++){
          iakk = ws2[i];
          if( iakk != -1) {
//          printf("%d %d %d %d %lf %lf %lf\n",j+1,k+1,i,iakk,w[i],m[iakk],tmp);
            w[i] -= m[iakk]*tmp; 
          }
        }
        for(i=iaUp3;i<iaUp4;i++){
          jak1      = jaUp[i];
          ws2[jak1] = -1;
        }
      }
    }
    timei = getTimeC() - timei;   
    
//  printf("%lf %lf %lf %lf %lf\n",w[0],w[1],w[2],w[3],w[4]);
    for(i=iaUp1,ii=0;i<iaUp2;i++,ii++){
      jak1        = jaUp[i];
      m[i]        = w[jak1];
    }
    
    for(i=iaLo[j];i<iaLo[j+1];i++){
      jak1        = jaLo[i];
      ws1[jak1]   = -1;
//       printf("%d %d %d\n",j,jak1,i);
    }

    for(i=iaUp1;i<iaUp2;i++){
//    printf("%lf\n",m[i]);
      jak1                      = jaUp[i];
      m[i]                     /= m[j+nad];
      m[nad+jak1]              -= m[i]*m[i];
      mLo[iaLo[jak1]+ws3[jak1]] = m[i];
      ws3[jak1]++;     
    }
/*  
    printf("LCsc\n");
    for(i=0;i<nad;i++){
      printf("%d %lf\n",i,m[i]);
    }
  
    printf("LCsr\n");
    for(i=0;i<nad;i++){
      printf("%d %lf\n",i,mLo[i]);
    }  
    printf("================================================\n");
*/   
  } 
  printf("%lf\n",timei);
//exit(0);
}

/**********************************************************************/
void ic0Csr(DOUBLE *restrict md    
           ,DOUBLE *restrict mCsrUp
           ,DOUBLE *restrict mCscUp
           ,DOUBLE *restrict ad
           ,DOUBLE *restrict aCsrUp
           ,DOUBLE *restrict w   
           ,INT *restrict iaCsrUp,INT *restrict jaCsrUp
           ,INT *restrict iaCscUp,INT *restrict jaCscUp
           ,INT    *restrict ws1
           ,INT    *restrict ws2
           ,INT    *restrict ws3
           ,INT const neq,INT const nad){
  
  int iaCsrUp1,iaCsrUp2;
  int iakk;
  int iaCsrUp3,iaCsrUp4;
  int jak1;
  int i,j,k,ii;
  double tmp = 0.0;
  double timei=0.0; 
  
  for(k = 0; k < neq;k++){
      md[k]  = ad[k];   
     ws1[k]  = -1;
     ws2[k]  = -1;
     ws3[k]  =  0;
     w  [k]  =  0.0;
  }
  
  for(k = 0; k < nad;k++){
     mCsrUp[k]   =  aCsrUp[k];   
  }

  for(j = 0; j< neq    ;j++){
    iaCsrUp1 = iaCsrUp[j  ];
    iaCsrUp2 = iaCsrUp[j+1];
    
    md[j] = sqrt(md[j]);
      
    for(i=iaCsrUp1,ii=0;i<iaCsrUp2;i++,ii++){
       jak1        = jaCsrUp[i];
         w[jak1]   = mCsrUp[i];
    }
    
    for(i=iaCscUp[j];i<iaCscUp[j+1];i++){
       jak1        = jaCscUp[i];
       ws1[jak1]   = i;
    }
    
    for(k=0;k<j;k++){
      if(ws1[k] != -1){ 
        iaCsrUp3        = iaCsrUp[k  ];
        iaCsrUp4        = iaCsrUp[k+1];
        for(i=iaCsrUp3;i<iaCsrUp4;i++){
          jak1         = jaCsrUp[i];
          ws2[jak1]    = i;
        }
        tmp          =   mCscUp[iaCscUp[j]+k];
        for(i=j+1;i<neq;i++){
          iakk = ws2[i];
          if( iakk != -1) {
            w[i] -= mCsrUp[iakk]*tmp; 
          }
        }
        for(i=iaCsrUp3;i<iaCsrUp4;i++){
          jak1      = jaCsrUp[i];
          ws2[jak1] = -1;
        }
      }
    }
    
    for(i=iaCsrUp1,ii=0;i<iaCsrUp2;i++,ii++){
      jak1        = jaCsrUp[i];
      mCsrUp[i]   = w[jak1];
    }
    
    for(i=iaCscUp[j];i<iaCscUp[j+1];i++){
      jak1        = jaCscUp[i];
      ws1[jak1]   = -1;
    }

    for(i=iaCsrUp1;i<iaCsrUp2;i++){
      jak1                              = jaCsrUp[i];
      mCsrUp[i]                        /= md[j];
      md[jak1]                         -= mCsrUp[i]*mCsrUp[i];
      mCscUp[iaCscUp[jak1]+ws3[jak1]]   = mCsrUp[i];
      ws3[jak1]++;     
    }
  } 

}

/**********************************************************************/
void ic0Csr2(DOUBLE *restrict md    
            ,DOUBLE *restrict mCsrUp
            ,DOUBLE *restrict mCscUp
            ,DOUBLE *restrict ad
            ,DOUBLE *restrict aCsrUp
            ,DOUBLE *restrict w   
            ,INT *restrict iaCsrUp,INT *restrict jaCsrUp
            ,INT *restrict iaCscUp,INT *restrict jaCscUp
            ,INT    *restrict ws1
            ,INT    *restrict ws2
            ,INT    *restrict ws3
            ,INT const neq,INT const nad){
  
  int iaCsrUp1,iaCsrUp2;
  int iakk;
  int iaCsrUp3,iaCsrUp4;
  int jak1;
  int i,j,k,ii;
  double tmp = 0.0;
  double timei=0.0; 
  
  for(k = 0; k < neq;k++){
      md[k]  = ad[k];   
     ws1[k]  = -1;
     ws2[k]  = -1;
     ws3[k]  =  0;
     w  [k]  =  0.0;
  }
  
  for(k = 0; k < nad;k++){
     mCsrUp[k]   =  aCsrUp[k];   
  }

  for(j = 0; j< neq    ;j++){
    iaCsrUp1 = iaCsrUp[j  ];
    iaCsrUp2 = iaCsrUp[j+1];
    
    md[j] = sqrt(md[j]);
      
    for(i=iaCsrUp1,ii=0;i<iaCsrUp2;i++,ii++){
       jak1        = jaCsrUp[i];
       w[jak1]     = mCsrUp[i];
    }
    
    for(i=iaCscUp[j];i<iaCscUp[j+1];i++){
       jak1        = jaCscUp[i];
       ws1[jak1]   = i;
    }
    
    for(k=0;k<j;k++){
      if(ws1[k] != -1){ 
        iaCsrUp3        = iaCsrUp[k  ];
        iaCsrUp4        = iaCsrUp[k+1];
        for(i=iaCsrUp3;i<iaCsrUp4;i++){
          jak1         = jaCsrUp[i];
          ws2[jak1]    = i;
        }
        tmp          =   mCscUp[iaCscUp[j]+k];
        for(i=j+1;i<neq;i++){
          iakk = ws2[i];
          if( iakk != -1) {
            w[i] -= mCsrUp[iakk]*tmp; 
          }
        }
        for(i=iaCsrUp3;i<iaCsrUp4;i++){
          jak1      = jaCsrUp[i];
          ws2[jak1] = -1;
        }
      }
    }
    
    for(i=iaCsrUp1,ii=0;i<iaCsrUp2;i++,ii++){
      jak1        = jaCsrUp[i];
      mCsrUp[i]   = w[jak1];
    }
    
    for(i=iaCscUp[j];i<iaCscUp[j+1];i++){
      jak1        = jaCscUp[i];
      ws1[jak1]   = -1;
    }

    for(i=iaCsrUp1;i<iaCsrUp2;i++){
      jak1                              = jaCsrUp[i];
      mCsrUp[i]                        /= md[j];
      md[jak1]                         -= mCsrUp[i]*mCsrUp[i];
      mCscUp[iaCscUp[jak1]+ws3[jak1]]   = mCsrUp[i];
      ws3[jak1]++;     
    }
  } 

}
    
