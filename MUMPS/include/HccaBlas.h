#ifndef _HCCABLAS_H
  #define _HCCABLAS_H

/*...*/  
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/*...*/
  long flopDot(INT nDim);
  void prodVet(DOUBLE *restrict a
              ,DOUBLE *restrict b
              ,DOUBLE *restrict c);
/*...................................................................*/


/*======================== level 1 =================================*/
/*... soma de vetores*/
  void addVector(DOUBLE const alpha,DOUBLE *restrict a
                ,DOUBLE const beta ,DOUBLE *restrict b
                ,INT const nDim    ,DOUBLE *restrict c);
/*... multiplicacao de vetor por um escalar*/  
  void  alphaProdVector(DOUBLE const alpha,DOUBLE *restrict a
                       ,INT const nDim    ,DOUBLE *restrict c);
/*...................................................................*/

/*... produto enterno*/
  DOUBLE dot(DOUBLE *restrict a,DOUBLE *restrict b,INT const nDim);
  DOUBLE dotO2L2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotL2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotL4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotL6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotL8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
/*...................................................................*/
/*==================================================================*/

/*======================== level 2 =================================*/
/*CSRD*/
  void matVecCsrD   (INT const neq           
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict a ,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void matVecCsrDSym(INT const neq           
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict al,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y);
  void matVecCsrDSymUpLower(INT const neq           
                           ,INT *restrict ia   ,INT *restrict ja
                           ,DOUBLE *restrict au,DOUBLE *restrict ad
                           ,DOUBLE *restrict x ,DOUBLE *restrict y);
/*==================================================================*/

/*======================== level 3 =================================*/
/*==================================================================*/

#endif
