#include<Solv.h>
//#include<mkl.h>
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
        ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict b0 
        ,DOUBLE *restrict x ,DOUBLE *restrict v
        ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
        ,unsigned int maxIt ,bool newX          
        ,FILE* fLog         ,bool log
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  int k=0;
  DOUBLE alpha,beta,d,conv,bNormA,rNormA,energy;
  DOUBLE timei,timef;
  timei = getTimeC();

/* chute inicial*/
  if(newX)  
    for(i = 0; i < neq; i++)  
      x[i] = 0.e0;
  
/*... ||b||a*/
  matvec(neq,ia,ja,al,ad,b,z);
  bNormA = sqrt(dot(b,z,neq));
/*........................................................*/
  
  matvec(neq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < neq; i++)   {
    r[i] = b[i] - z[i];
    z[i] = r[i] / m[i];
    b[i] = z[i];
  }
  d     = dot(r,z,neq);
  conv  = tol * sqrt(fabs(d));
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
    if(log){
/*.. ||b-Axit||a*/
      matvec(neq,ia,ja,al,ad,x,z);
      for(i = 0; i < neq; i++)   
        z[i] = b0[i] - z[i];
      matvec(neq,ia,ja,al,ad,z,v);
      rNormA = sqrt(dot(z,v,neq));
/*........................................................*/
      fprintf(fLog,"%9d %e %e %e\n"          
             ,j+1,sqrt(fabs(d)),alpha,rNormA/bNormA);
    }
  }
/* -------------------------------------------------------*/
  matvec(neq,ia,ja,al,ad,x,z);
/*norma de energia = xTAx* */
  energy = dot(x,z,neq);
/* -------------------------------------------------------*/
  
/*...*/
  for(i = 0; i < neq; i++)
    r[i]  = b0[i] - z[i];
  
  matvec(neq,ia,ja,al,ad,r,z);
  printf("\t||b - Ax|| = %20.2e\n"  , sqrt(dot(r,z,neq)));
  printf("\t||b||      = %20.2e\n"  , bNormA);
/*.........................................................*/

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
/*if(log)
    fprintf(fLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);*/
}
/**********************************************************************/

/**********************************************************************
 * DPCG  : metodo do gradiente conjugado com precondiconador diagonal  *
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
void dpcg(INT const neq      ,INT const nad  
         ,INT *restrict ia   ,INT *restrict ja
         ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
         ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict b0  
         ,DOUBLE *restrict x
         ,DOUBLE *restrict s ,DOUBLE *restrict r ,DOUBLE *restrict v
         ,DOUBLE const tol
         ,DOUBLE *restrict zt        
         ,DOUBLE *restrict e         ,DOUBLE *restrict ztz
         ,DOUBLE *restrict az        ,DOUBLE *restrict x0
         ,DOUBLE *restrict q1        ,DOUBLE *restrict q2                  
         ,INT const nDVector                              
         ,unsigned int maxIt ,bool newX          
         ,FILE* fLog         ,bool log
         ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  INT    *iaAzCsr=NULL,*jaAzCsr=NULL,nnzAz;
  DOUBLE *AzCsr=NULL;
  int k=0;
  bool flag = false;
  DOUBLE alpha,beta,d,conv,energy;
  DOUBLE timei,timef,timeE,timeCh,timeChS,timeBlasZt,timeBlasAz,timeBlasT;
  timei = getTimeC();

  
/*... E = zT*A*z */
  timeE = getTimeC();   
  mkZtAz(e     ,az 
        ,zt    ,s 
        ,ad    ,al                   
        ,ia    ,ja
        ,neq   ,nDVector
        ,matvec,false);
  timeE = getTimeC() - timeE;   
/*........................................................*/

/**/
  nnzAz = countAzNz(az,neq,nDVector,1.0e-16);
  iaAzCsr = (INT*)    malloc(sizeof(INT)*(neq+1));
  jaAzCsr = (INT*)    malloc(sizeof(INT)*nnzAz);
  AzCsr   = (DOUBLE*) malloc(sizeof(DOUBLE)*nnzAz);
  mktAzCsr(az,neq,nDVector,iaAzCsr,jaAzCsr,AzCsr,1.0e-16);
/*........................................................*/

/*... fatoracao GGt*/
  timeCh= getTimeC();   
  fatGGt(e,nDVector);
//fatLU(e,nDVector,LUKIJ);
  timeCh= getTimeC() - timeCh;   
/*........................................................*/
 
/*...*/
  mkWtw(zt,ztz,s,v,nDVector,neq,dot);
/*... fatoracao wTw*/
  timeCh= getTimeC();   
  fatGGt(ztz,nDVector);
//fatLU(e,nDVector,LUKIJ);
  timeCh= getTimeC() - timeCh;   
/*........................................................*/

/* 
  for(i=0;i<nDVector;i++){
    for(j=0;j<nDVector;j++)
      printf("%lf ",MAT2D(i,j,e,nDVector));
    printf("\n");
  }  
*/

/* chute inicial*/
  if(newX)  
    for(i = 0; i < neq; i++)  
      x[i] = 0.e0;
      
  matvec(neq,ia,ja,al,ad,x,s);
  
  for(i = 0; i < neq; i++)   {
    r[i]  = b[i] - s[i];
    x0[i] = x[i];
    x[i]  = 0.e0;
  }

/*...zT*r(0)*/
  timeBlasZt= getTimeC();   
//matVecFull(NOTRANSP,dfT,r,q1,nDVector,neq);
  matVecFullO2(zt,r,q1,nDVector,neq);
  timeBlasZt= getTimeC() - timeBlasZt;   
/*........................................................*/
  
/*... solver: E q1=zt*r(0)*/
  timeChS =  getTimeC();
  solvCholesky(e,q1,nDVector);  
//solvLU(e,q1,nDVector);
  timeChS = getTimeC() -timeChS ;
/*........................................................*/

/*...r = r - (az)q1*/
  timeBlasAz= getTimeC();   
//matVecFull(NOTRANSP,aDf,q1,z,neq,nDVector);
//matVecFullO2(az,q1,s,neq,nDVector);
  matVecCsr(neq           
           ,iaAzCsr  ,jaAzCsr
           ,AzCsr    ,q1
           ,s);
  timeBlasAz= getTimeC() - timeBlasAz;   
  
  for(i = 0; i < neq; i++)   {
    r[i]  -= s[i];
/*... predcondicionador diagonal*/
    s[i]  = r[i]/m[i];
    b[i]  = s[i];
  }
/*........................................................*/
  d    = dot(r,s,neq);
  conv = tol * sqrt(fabs(d));
/*...*/   
  for(j = 0; j < maxIt; j++)   {

    matvec(neq,ia,ja,al,ad,b,s);

/*...zt*s*/
    timeBlasZt= getTimeC() - timeBlasZt;   
//  matVecFull(NOTRANSP,dfT,z,q2,nDVector,neq);
    matVecFullO2I2(zt,s,q2,nDVector,neq);
//    cblas_dgemv(CblasRowMajor
//               ,CblasNoTrans,nDVector,neq,1.0e0,zt
//              ,neq,s,1,0.e0,q2,1);
    timeBlasZt= getTimeC() - timeBlasZt;   
/*........................................................*/
  
/*... solver: E q3=zT*s*/
    timeChS =  getTimeC() - timeChS;
    solvCholesky(e,q2,nDVector);  
//  solvLU(e,q2,nDVector);
    timeChS =  getTimeC() - timeChS;
/*........................................................*/
    
/*...v = v - (az)q3*/
    timeBlasAz= getTimeC() - timeBlasAz;   
//  matVecFull(NOTRANSP,aDf,q2,v,neq,nDVector);
//  matVecFullO2I2(az,q2,v,neq,nDVector);
    matVecCsr(neq           
             ,iaAzCsr  ,jaAzCsr
             ,AzCsr    ,q2
             ,v);
    for(i = 0; i < neq; i++)   
      s[i]  -= v[i];
//  cblas_dgemv(CblasRowMajor
//             ,CblasNoTrans,neq,nDVector,-1.0e0,az  
//            ,nDVector,q2,1,1.e0,s,1);
    timeBlasAz= getTimeC() - timeBlasAz;   
/*........................................................*/
    
    alpha = d / dot(b,s,neq);
    
/*... reortogonalizado de rj em relacao w*/
    if(flag) {
      for(i = 0; i < neq; i++)   {
        x[i] +=  alpha * b[i];
        r[i] -=  alpha * s[i];
      }
/*... wtrj*/
      matVecFull(NOTRANSP,zt,r,q1,nDVector,neq);
/*... (wtw)*q = wt*rj*/
      solvCholesky(ztz,q1,nDVector);
/*... wq*/
      matVecFull(TRANSP,zt,q1,s,neq,nDVector);
      for(i = 0; i < neq; i++)   {
        r[i]  -= s[i];
/*... precondicionador diagonal*/
        s[i]  =   r[i] / m[i];
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sem reortogonalizado de rj em relacao w*/
    else{
      for(i = 0; i < neq; i++){
        x[i] +=  alpha * b[i];
        r[i] -=  alpha * s[i];
/*... precondicionador diagonal*/
        s[i]  =   r[i] / m[i];
      }
/*...................................................................*/
    } 

    beta = dot(r,s,neq)/d;

    for(i = 0; i < neq; i++)   {
      b[i] = s[i] + beta * b[i];
    }

    d = beta * d;

    if (sqrt(fabs(d)) < conv) break;
    if( k == 1000){ 
      printf("it: %d %e %e\n",j,sqrt(fabs(d)),conv);
      k = 0; 
    }
    k++;
    if(log)
      fprintf(fLog,"%9d %e %e %e\n"          
             ,j+1,sqrt(fabs(d)),alpha,sqrt(dot(b,b,neq)));
  }
/*........................................................*/

/*... Ax*/
  matvec(neq,ia,ja,al,ad,x,s);

/*...zT*s*/
  timeBlasZt= getTimeC() - timeBlasZt;   
  matVecFull(NOTRANSP,zt,s,q2,nDVector,neq);
  timeBlasZt= getTimeC() - timeBlasZt;   
/*........................................................*/
  
/*... solver: Eq2=zT*s*/
  timeChS =  getTimeC() - timeChS;
  solvCholesky(e,q2,nDVector); 
//solvLU(e,q2,nDVector);
  timeChS =  getTimeC() - timeChS;
/*........................................................*/
    
/*... */
  for(i = 0; i < nDVector; i++)   
    q1[i]  -= q2[i];
/*........................................................*/
    
  timeBlasT= getTimeC();  
  matVecFull(TRANSP,zt,q1,v,neq,nDVector);
  timeBlasT= getTimeC() - timeBlasT;   

/*... */
  for(i = 0; i < neq; i++)   
    x[i]  += x0[i] + v[i];
/*........................................................ */

/* ------------------------------------------------------- */
  matvec(neq,ia,ja,al,ad,x,s);
/*norma de energia = xTAx* */
  energy = dot(x,s,neq);
/* ....................................................... */

/*...*/
  for(i = 0; i < neq; i++)
    r[i]  = b0[i] - s[i];
  matvec(neq,ia,ja,al,ad,r,s);
  printf("\t||b - Ax|| = %20.2e\n"  , sqrt(dot(r,s,neq)));
/* .......................................................*/
  timef = getTimeC() - timei;   
  printf("\tnad         :      %20d\n"  ,nad);
  printf("\tSolver conv :      %20.2e\n",conv);
  printf("\tSolver tol  :      %20.2e\n",tol);
  printf(" (DPCG1) solver:\n"
         "\tEquations   =      %20d\n"
         "\tIterarions  =      %20d\n"
	 "\tEnergy norm  =      %20.12e\n"
	 "\tCPU    time(s) =      %20.5lf\n" 
	 "\tBlasZt time(s) =      %20.5lf\n" 
	 "\tBlasAz time(s) =      %20.5lf\n" 
	 "\tBlasT  time(s) =      %20.5lf\n" 
	 "\tChS    time(s) =      %20.5lf\n" 
	 "\tCh     time(s) =      %20.5lf\n" 
	 "\tE      time(s) =      %20.5lf\n" 
	 ,neq,j+1,energy,timef,timeBlasZt,timeBlasAz,timeBlasT,timeChS,timeCh,timeE);
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
//    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"DPCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/

/**********************************************************************
 * DPCG2 : metodo do gradiente conjugado com precondiconador diagonal  *
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
void dpcg2(INT const neq      ,INT const nad  
         ,INT *restrict ia   ,INT *restrict ja
         ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
         ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict b0
         ,DOUBLE *restrict x
         ,DOUBLE *restrict s ,DOUBLE *restrict r ,DOUBLE *restrict v
         ,DOUBLE const tol
         ,DOUBLE *restrict wt        
         ,DOUBLE *restrict e         ,DOUBLE *restrict wtw
         ,DOUBLE *restrict aw        ,DOUBLE *restrict x0
         ,DOUBLE *restrict q1        ,DOUBLE *restrict q2                  
         ,INT const nDVector                              
         ,unsigned int maxIt ,bool newX          
         ,FILE* fLog         ,bool log
         ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  INT    *iaWtACsr=NULL,*jaWtACsr=NULL,nnzWtA;
  DOUBLE *wTaCsr=NULL;
  int k=0;
  bool flag=true;
  DOUBLE alpha,beta,d,bNormA,rNormA,conv,energy;
  DOUBLE timei,timef,timeE,timeCh,timeChS,timeBlasZt,timeBlasAz,timeBlasT;
  timei = getTimeC();

/*... E = wT*A*w */
  timeE = getTimeC();   
  mkZtAz(e     ,aw 
        ,wt    ,s 
        ,ad    ,al                   
        ,ia    ,ja
        ,neq   ,nDVector
        ,matvec,true);
  timeE = getTimeC() - timeE;   
/*........................................................*/

/**/
//  nnzWtA = countAzNz(wTa,nDVector,neq,1.0e-16);
//  iaWtACsr = (INT*)   malloc(sizeof(INT)*(neq+1));
//  jaWtACsr = (INT*)    malloc(sizeof(INT)*nnzAz);
//  wTaCsr   = (DOUBLE*) malloc(sizeof(DOUBLE)*nnzAz);
//  mktAzCsr(wTa,nDVector,neq, iaWtACsr,jaWtACsr,wTaCsr,1.0e-16);
/*........................................................*/

/*... fatoracao GGt*/
  timeCh= getTimeC();   
  fatGGt(e,nDVector);
//fatLU(e,nDVector,LUKIJ);
  timeCh= getTimeC() - timeCh;   
/*........................................................*/

/*...*/
  mkWtw(wt,wtw,s,v,nDVector,neq,dot);
/*... fatoracao wTw*/
  timeCh= getTimeC();   
  fatGGt(wtw,nDVector);
//fatLU(e,nDVector,LUKIJ);
  timeCh= getTimeC() - timeCh;   
/*........................................................*/

/* chute inicial*/
  if(newX)  
    for(i = 0; i < neq; i++)  
      x[i] = 0.e0;

/*... ||b||a*/
  matvec(neq,ia,ja,al,ad,b,v);
  bNormA = sqrt(dot(b,v,neq));
/*........................................................*/
    
  matvec(neq,ia,ja,al,ad,x,v);
  
  for(i = 0; i < neq; i++)   {
    r[i]  = b[i] - v[i];
  }
  matvec(neq,ia,ja,al,ad,x,v);
  
 
/* 
  for(i=0;i<nDVector;i++){
    for(j=0;j<nDVector;j++)
      printf("%lf ",MAT2D(i,j,e,nDVector));
    printf("\n");
  }  
*/

/*...wT*r(0)*/
  timeBlasZt= getTimeC();   
  matVecFull(NOTRANSP,wt,r,q1,nDVector,neq);
  timeBlasZt= getTimeC() - timeBlasZt;   
/*........................................................*/
  
/*... solver: E q1=wt*r(0)*/
  timeChS =  getTimeC();
  solvCholesky(e,q1,nDVector);  
//solvLU(e,q1,nDVector);
  timeChS = getTimeC() -timeChS ;
/*........................................................*/

/*... x0 = x-1 + W*q1*/  
  matVecFull(TRANSP,wt,q1,v,neq,nDVector);
  
  for(i = 0; i < neq; i++)   
   x[i] += v[i]; 
/*........................................................*/

/*... A*x0*/
  matvec(neq,ia,ja,al,ad,x,s);
  for(i = 0; i < neq; i++)   {
/*... r0 = b - A*x0*/
    r[i]  = b[i] - s[i];
/*... predcondicionador diagonal*/
    s[i]  = r[i]/m[i];
    b[i]  = s[i];
  }
/*........................................................*/

/*... Ad=wtA*z*/
    matVecFull(TRANSP,aw,s,q1,nDVector,neq);
//  matVecCsr(nDVector           
//           ,iaWtACsr ,jaWtACsr
//           ,wTaCsr   ,q2
//           ,v);
  solvCholesky(e,q1,nDVector);
/*...p0 = -zd+s*/  
  matVecFull(TRANSP,wt,q1,v,neq,nDVector);
  for(i = 0; i < neq; i++)
    b[i] = -v[i] + s[i];  
/*........................................................*/
  d    = dot(r,s,neq);
  conv = tol * sqrt(fabs(d));

/*...*/   
  for(j = 0; j < maxIt; j++)   {

    matvec(neq,ia,ja,al,ad,b,s);
    
    alpha = d / dot(b,s,neq);

    
/*... reortogonalizado de rj em relacao w*/
    if(flag) {
      for(i = 0; i < neq; i++)   {
        x[i] +=  alpha * b[i];
        r[i] -=  alpha * s[i];
      }
/*... wtrj*/
      matVecFull(NOTRANSP,wt,r,q1,nDVector,neq);
/*... (wtw)*q = wt*rj*/
      solvCholesky(wtw,q1,nDVector);
/*... wq*/
      matVecFull(TRANSP,wt,q1,s,neq,nDVector);
      for(i = 0; i < neq; i++)   {
        r[i]  -=  s[i];
/*... precondicionador diagonal*/
        s[i]  =   r[i] / m[i];
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sem reortogonalizado de rj em relacao w*/
    else{
      for(i = 0; i < neq; i++){
        x[i] +=  alpha * b[i];
        r[i] -=  alpha * s[i];
/*... precondicionador diagonal*/
        s[i]  =   r[i] / m[i];
      }
/*...................................................................*/
    } 

/*... Ad=ztA*s*/
    matVecFull(TRANSP,aw,s,q1,nDVector,neq);
    solvCholesky(e,q1,nDVector);
/*... W*d*/
    matVecFull(TRANSP,wt,q1,v,neq,nDVector);
    
    beta = dot(r,s,neq)/d;

    for(i = 0; i < neq; i++)   {
      b[i] = s[i] + beta * b[i] - v[i];
    }

    d = beta * d;

    if (sqrt(fabs(d)) < conv) break;
    if( k == 1000){ 
      printf("it: %d %e %e\n",j,sqrt(fabs(d)),conv);
      k = 0; 
    }
    k++;
    if(log){
/*.. ||b-Axit||a*/
      matvec(neq,ia,ja,al,ad,x,s);
      for(i = 0; i < neq; i++)
        s[i] = b0[i] - s[i];
      matvec(neq,ia,ja,al,ad,r,v);
      rNormA = sqrt(dot(s,v,neq));
/*........................................................*/
      fprintf(fLog,"%9d %e %e %e\n"          
             ,j+1,sqrt(fabs(d)),alpha
             ,rNormA/bNormA);
    }
  }
/*........................................................*/

/*... Ax*/
  matvec(neq,ia,ja,al,ad,x,s);
/*norma de energia = xTAx* */
  energy = dot(x,s,neq);
/* .......................................................*/

/*...*/
  for(i = 0; i < neq; i++)
    r[i]  = b0[i] - s[i];
  matvec(neq,ia,ja,al,ad,r,s);
  printf("\t||b - Ax|| = %20.2e\n"  , sqrt(dot(r,s,neq)));
  printf("\t||b||      = %20.2e\n"  , bNormA);
/* .......................................................*/


  timef = getTimeC() - timei;   
  printf("\tnad         :      %20d\n"  ,nad);
  printf("\tSolver conv :      %20.2e\n",conv);
  printf("\tSolver tol  :      %20.2e\n",tol);
  printf(" (DPCG2) solver:\n"
         "\tEquations   =      %20d\n"
         "\tIterarions  =      %20d\n"
	 "\tEnergy norm  =      %20.12e\n"
	 "\tCPU    time(s) =      %20.5lf\n" 
	 "\tBlasZt time(s) =      %20.5lf\n" 
	 "\tBlasAz time(s) =      %20.5lf\n" 
	 "\tBlasT  time(s) =      %20.5lf\n" 
	 "\tChS    time(s) =      %20.5lf\n" 
	 "\tCh     time(s) =      %20.5lf\n" 
	 "\tE      time(s) =      %20.5lf\n" 
	 ,neq,j+1,energy,timef,timeBlasZt,timeBlasAz,timeBlasT,timeChS,timeCh,timeE);
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
//    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"DPCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/


/**********************************************************************
 * IC0CG : metodo do gradiente conjugado com precondiconador IC0      *
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
 *mIcd  -> martiz do precondicionador                                 *
 *mIcCsr-> matriz do precondicionador                                 *
 *mIcCsc-> matriz do precondicionador                                 *
 *iaMCsr-> ponteiro da matriz do precondicionador                     *
 *jaMCsr-> ponteiro da matriz do precondicionador                     *
 *iaMCsc-> ponteiro da matriz do precondicionador                     *
 *jaMCsc-> ponteiro da matriz do precondicionador                     *
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

void mkZtAz(DOUBLE *restrict e,DOUBLE *restrict az
           ,DOUBLE *restrict zt,DOUBLE *restrict y 
           ,DOUBLE *restrict ad,DOUBLE *restrict al                   
           ,INT *restrict ia   ,INT *restrict ja
           ,INT const neq      ,INT const nDVector
           ,void(*matvec)()    ,bool const cod)
{
  int i,j;
  int dl;
  double tmp;
/*... calculo da matriz az*/
  for(i=0;i<nDVector;i++){
    dl  = i*neq;
    matvec(neq,ia,ja,al,ad,&zt[dl],y);
    for(j=0;j<neq;j++){
      MAT2D(j,i,az,nDVector)  =  y[j];
    }
  }
/*... caclulo da matriz e = ZtAz*/
  hccaDgemm(nDVector,neq,nDVector,zt,az,e);   
/*.....................................................................*/ 
  
/*... calculo da matriz zta*/
//  if(cod){
//    for(i=0;i<neq;i++){
//      for(j=i;j<nDVector;j++){
//        tmp               = MAT2D(j,i,az,nDVector);
//        MAT2D(i,j,az,neq) = tmp;
//
//      }
//    }

//  }
/*.....................................................................*/
}
   
 
INT countAzNz(DOUBLE *restrict az,INT nl,INT nc,DOUBLE tol)
{

  INT i,j,count;
  DOUBLE tmp;

  count = 0;
  for(i=0;i<nl;i++){
    for(j=0;j<nc;j++){
      tmp = MAT2D(i,j,az,nc);
      if( fabs(tmp) > tol) 
        count++;
    }
  }
  return count;

}
void  mktAzCsr(DOUBLE *az ,INT nl      , INT nc         
              ,INT *iaAzCsr,INT *jaAzCsr,DOUBLE *azCsr
              ,DOUBLE tol)
{
  INT i,j,count;
  DOUBLE tmp;

  count = 0;
  for(i=0;i<nl;i++){
    iaAzCsr[i] = count;
    for(j=0;j<nc;j++){
      tmp = MAT2D(i,j,az,nc);
      if( fabs(tmp) > tol){
        azCsr[count]   = tmp;
        jaAzCsr[count] = j;
        count++;
      }
    }
  }
  iaAzCsr[i] = count;
/*
  for(i=0;i<nl+1;i++){
    printf("ia %d\n",iaAzCsr[i]);
  } 
  for(i=0;i<iaAzCsr[nl];i++){
    printf("ja %d\n",jaAzCsr[i]);
  } 
*/
}

void  mkWtw(DOUBLE *restrict wt , DOUBLE *restrict wtw
           ,DOUBLE *restrict c1 ,DOUBLE *restrict c2
           ,INT const nDVector ,INT const nEq
           ,DOUBLE(*dot)()){
  INT i,j,k;

/*...*/
  for(i=0;i<nDVector;i++){
/*... w(i,*)*/
    for(k=0;k<nEq;k++)
      c1[k] = MAT2D(i,k,wt,nEq);
    for(j=0;j<nDVector;j++){
/*... w(j,*)*/
      for(k=0;k<nEq;k++)
        c2[k] = MAT2D(j,k,wt,nEq);
      MAT2D(i,j,wtw,nDVector) = dot(c1,c2,nEq);
    }
  }
/*...................................................................*/

}

