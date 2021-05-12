#ifndef _HCCABLAS_
  #define _HCCABLAS_
  #define NFUNC           159
  #define HCCABLASZERO 1.0e-10
/*...*/
  #include<Define.h>
/*...*/
  #include<stdio.h>
  
  #ifdef _OPENMP
    #include<omp.h>
    DOUBLE tmpDotOmp;
    DOUBLE tmpDotOmp1;
    DOUBLE tmpDotOmp2;
    DOUBLE tmpDotOmp3;
    DOUBLE tmpDotOmp4;
    DOUBLE tmpDotOmp5;
    DOUBLE tmpDotOmp6;
    DOUBLE tmpDotOmp7;
    DOUBLE tmpDotOmp8;
  #endif

/*...*/
  void getNameHccaBlas(void);
  long flopMatVecFull(INT nLin,INT nCol);
  long flopMatVecCsr(INT neq,INT nad,short ty);
  long flopDot(INT nDim);
/*produto vetorial*/  
  void prodVet(DOUBLE *restrict a,DOUBLE *restrict b
              ,DOUBLE *restrict c);
  INT xDiffY(DOUBLE *restrict x,DOUBLE *restrict y
            ,DOUBLE const tol  , INT n);
/*...................................................................*/

/*... normas*/
  LDOUBLE lnormInf(LDOUBLE *restrict x,INT const lin,INT const col);
  DOUBLE normOne(DOUBLE *restrict x,INT const lin,INT const col);
  DOUBLE normInf(DOUBLE *restrict x,INT const lin,INT const col);
/*...................................................................*/

/*level 1*/
  void alphaProdVector(DOUBLE const alpha,DOUBLE *restrict a
                      ,INT const nDim    ,DOUBLE *restrict c); 
  void addVector(DOUBLE const alpha,DOUBLE *restrict a
                ,DOUBLE const beta ,DOUBLE *restrict b
                ,INT const nDim    ,DOUBLE *restrict c); 

  DOUBLE dot(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE dotO2L2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotL2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotL4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotL6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotL8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
#if _OPENMP
  DOUBLE     dotOmp(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE dotOmpO2L2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE dotOmpO2L4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpL2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpL4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpL6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpL8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
  DOUBLE   dotOmpO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT n);
#endif
/*...................................................................*/

/*level 2*/
/* ... matriz cheia*/
  void matVecFull(short const code
                 ,DOUBLE *restrict a
                 ,DOUBLE *restrict x
                 ,DOUBLE *restrict y
                 ,INT nLin          ,INT nCol);
  void matVecFullO2(DOUBLE *restrict a
                 ,DOUBLE *restrict x
                 ,DOUBLE *restrict y
                 ,INT nLin          ,INT nCol);
  void matVecFullO4(DOUBLE *restrict a
                   ,DOUBLE *restrict x
                   ,DOUBLE *restrict y
                   ,INT nLin          ,INT nCol);
  void matVecFullO2I2(DOUBLE *restrict a
                     ,DOUBLE *restrict x
                     ,DOUBLE *restrict y
                     ,INT nLin          ,INT nCol);
  void matVecFullO4I2(DOUBLE *restrict a
                     ,DOUBLE *restrict x
                     ,DOUBLE *restrict y
                     ,INT nLin          ,INT nCol);
  void matVecFullO4I4(DOUBLE *restrict a
                     ,DOUBLE *restrict x
                     ,DOUBLE *restrict y
                     ,INT nLin          ,INT nCol);
  void lmatVecFull(LDOUBLE *restrict a
                  ,LDOUBLE *restrict x
                  ,LDOUBLE *restrict y
                  ,INT nLin       
                  ,INT nCol);
/*...................................................................*/

/*...Csr*/
//typedef enum {csr=1,csrD=2,csrC=3} typeCsr;

/*... CsrD*/ 
  void     matVecCsrD(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void matVecCsrDO2I2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
/*... matriz simetrica*/
  void matVecCsrDSym(INT const neq           
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict al,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void matVecCsrDSymUpLower(INT const neq           
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict au,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y);

#ifdef _OPENMP
/*... CsrDOmp*/ 
  void       matVecCsrDOmp(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpI2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpI4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpI6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpO2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpO4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpO6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void       matVecCsrDOmpO2I2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
/*balanciamento manual do trabalho entre as threads*/  
  void matVecCsrDOmpBal(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO4(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO6(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO2I2(INT neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y
                     ,INT  *thBegin     ,INT *thEnd  );


#endif         
/*... CsrC*/ 
  void     matVecCsrC(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCI2(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCI4(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCI6(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCO2(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCO4(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void   matVecCsrCO6(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void matVecCsrCO2I2(INT neq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict au,DOUBLE *restrict ad
                     ,DOUBLE *restrict al
                     ,DOUBLE *restrict x ,DOUBLE *restrict y);
#ifdef _OPENMP
  void     matVecCsrCOmp(INT neq           
                        ,INT *restrict ia    ,INT *restrict ja
                        ,DOUBLE *restrict au ,DOUBLE *restrict ad
                        ,DOUBLE *restrict al
                        ,DOUBLE *restrict x  ,DOUBLE *restrict y
                        ,INT  *thBegin       ,INT *thEnd  
                        ,INT  *thHeight    
                        ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpI2(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpI4(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpI6(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpO2(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpO4(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpO6(INT neq           
                          ,INT *restrict ia    ,INT *restrict ja
                          ,DOUBLE *restrict au ,DOUBLE *restrict ad
                          ,DOUBLE *restrict al
                          ,DOUBLE *restrict x  ,DOUBLE *restrict y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *restrict thY,int nThreads);          
  void     matVecCsrCOmpO2I2(INT neq           
                            ,INT *restrict ia    ,INT *restrict ja
                            ,DOUBLE *restrict au ,DOUBLE *restrict ad
                            ,DOUBLE *restrict al
                            ,DOUBLE *restrict x  ,DOUBLE *restrict y
                            ,INT  *thBegin       ,INT *thEnd  
                            ,INT  *thHeight    
                            ,DOUBLE *restrict thY,int nThreads);          
#endif         
/*... Csr*/ 
  void     matVecCsr(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrI2(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrI4(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrI6(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrO2(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrO4(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void   matVecCsrO6(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
  void matVecCsrO2I2(INT neq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict x
                    ,DOUBLE *restrict y);
#ifdef _OPENMP
/*... CsrOmp*/ 
  void       matVecCsrOmp(INT neq           
                         ,INT *restrict ia  ,INT *restrict ja
                         ,DOUBLE *restrict a,DOUBLE *restrict x
                         ,DOUBLE *restrict y);
  void       matVecCsrOmpI2(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void       matVecCsrOmpI4(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void       matVecCsrOmpI6(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void       matVecCsrOmpO2(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void       matVecCsrOmpO4(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void       matVecCsrOmpO6(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
  void     matVecCsrOmpO2I2(INT neq           
                           ,INT *restrict ia  ,INT *restrict ja
                           ,DOUBLE *restrict a,DOUBLE *restrict x
                           ,DOUBLE *restrict y);
#endif
/*...................................................................*/

/*level 3*/
  void hccaDgemm(INT const ni,INT const nj,INT const nk
            ,DOUBLE *restrict a,DOUBLE *restrict b,DOUBLE *restrict c);
/*...................................................................*/

#endif
