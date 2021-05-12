#ifndef _SOLV_H
  #define _SOLV_H
/*... solver*/
  #define PCG        1
  #define PBICGSTAB  2
/*...................................................................*/

/*...................................................................*/

/*...*/
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<HccaBlas.h>
  #include<Define.h>
  #include<Memoria.h>
  #include<Gauss.h>
/*...................................................................*/
 
/*....*/
  typedef struct Solv{
    DOUBLE            tol;
    unsigned int    maxIt;
    short          solver;
    FILE        *fileSolv;
    bool              log;
  }Solv;
/*...................................................................*/

/*....*/
  void solverC(Memoria *m          ,INT const neq     ,INT const nad
              ,INT *ia             ,INT *ja   
              ,DOUBLE *al          ,DOUBLE *ad        ,DOUBLE *au
              ,DOUBLE *b           ,DOUBLE *x
              ,DOUBLE tol          ,unsigned int maxit
              ,short const storage ,short const solver
              ,FILE* fileSolvLog   ,bool const fLog
              ,bool const newX     ,bool const openmp  
              ,bool const unsym    ,bool const loopwise);
/*...................................................................*/

/*========================= Iterativos ==============================*/
/*... gradiente conjugado precondicionado*/
  void pcg(INT const neq      ,INT const nad  
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict b0 
          ,DOUBLE *restrict x ,DOUBLE *restrict v 
          ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
          ,unsigned int maxit ,bool newX          
          ,FILE* fileSolvLog  ,bool log
          ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado precondicionado*/
  void dpcg(INT const neq       ,INT const nad  
           ,INT *restrict ia    ,INT *restrict ja
           ,DOUBLE *restrict al ,DOUBLE *restrict ad,DOUBLE *restrict au
           ,DOUBLE *restrict m  ,DOUBLE *restrict b ,DOUBLE *restrict b0  
           ,DOUBLE *restrict x
           ,DOUBLE *restrict z  ,DOUBLE *restrict r ,DOUBLE *restrict v
           ,DOUBLE const tol
           ,DOUBLE *restrict df 
           ,DOUBLE *restrict e  ,DOUBLE *restrict wtw
           ,DOUBLE *restrict aDf,DOUBLE *restrict x0 
           ,DOUBLE *restrict q1 ,DOUBLE *restrict q2 
           ,INT const nDVector    
           ,unsigned int maxit ,bool newX          
           ,FILE* fileSolvLog  ,bool log
           ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado precondicionado*/
  void dpcg2(INT const neq       ,INT const nad  
            ,INT *restrict ia    ,INT *restrict ja
            ,DOUBLE *restrict al ,DOUBLE *restrict ad,DOUBLE *restrict au
            ,DOUBLE *restrict m  ,DOUBLE *restrict b ,DOUBLE *restrict b0 
            ,DOUBLE *restrict x
            ,DOUBLE *restrict z  ,DOUBLE *restrict r ,DOUBLE *restrict v
            ,DOUBLE const tol
            ,DOUBLE *restrict df  
            ,DOUBLE *restrict e  ,DOUBLE *restrict wtw
            ,DOUBLE *restrict aDf,DOUBLE *restrict x0 
            ,DOUBLE *restrict q1 ,DOUBLE *restrict q2 
            ,INT const nDVector    
            ,unsigned int maxit ,bool newX          
            ,FILE* fileSolvLog  ,bool log
            ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/


/*... gradiente conjugado precondicionado IC(0)*/
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
          ,void(*matvec)()    ,DOUBLE(*dot)());
/*===================================================================*/

/*========================= Diretos =================================*/
/*... solver direto para fatoracao incompleta de cholesky (IC(0) */
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
              ,INT const nad);
/*...................................................................*/
/*===================================================================*/

/*=================== Fatoracao incompleta ==========================*/
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
           ,INT const neq,INT const nad);
/*...................................................................*/
/*===================================================================*/

void mkZtAz(DOUBLE *restrict e,DOUBLE *restrict az
           ,DOUBLE *restrict zt,DOUBLE *restrict y 
           ,DOUBLE *restrict ad,DOUBLE *restrict al                   
           ,INT *restrict ia   ,INT *restrict ja
           ,INT const neq      ,INT const nDVector
           ,void(*matvec)()    ,bool const cod);

INT countAzNz(DOUBLE *restrict az,INT nl,INT nc,DOUBLE tol);

void  mktAzCsr(DOUBLE *az ,INT neq     , INT nDVector
              ,INT *iaAzCsr,INT *jaAzCsr,DOUBLE *azCsr
              ,DOUBLE tol);

void  mkWtw(DOUBLE *restrict wt , DOUBLE *restrict wtw
           ,DOUBLE *restrict c1 ,DOUBLE *restrict c2
           ,INT const nDVector ,INT const nEq
           ,DOUBLE(*dot)());



#endif/*_SOLV_H*/
