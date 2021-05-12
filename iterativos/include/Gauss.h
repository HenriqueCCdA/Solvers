#ifndef _GAUSS_H
  #define _GAUSS_H
/*...*/
  #include<math.h>
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #define NOTRANSP 1
  #define TRANSP   0
/*...................................................................*/

/*...*/
  #define LDMT           1
  #define LDLT           2
  #define GGT            3
  #define LUKIJ          4
  #define LUIKJ          5
  #define LUKJI          6
  #define LUJKI          7
  #define DOOLITTLEL     8
  #define DOOLITTLELC    9
  #define DOOLITTLECL   10
  #define LUKIJPP       11
  #define DOOLITTLECLPP 12
/*...................................................................*/

/*... */
  #include<Define.h>
  #include<HccaBlas.h>
  #include<HccaStdBool.h>
/*...................................................................*/
 

/*...*/
  void printFat(DOUBLE *a,DOUBLE *d,INT const n,short code);
/*...................................................................*/

/*...*/
  void initIterImprovmentData(LDOUBLE *restrict ar,LDOUBLE *restrict br
                             ,DOUBLE  *restrict a ,DOUBLE *restrict b
                             ,INT const nEq);
  void iterImprovement(LDOUBLE *restrict ar ,LDOUBLE *restrict br
                      ,LDOUBLE *restrict itr,LDOUBLE *restrict xr
                      ,DOUBLE  *restrict x
                      ,DOUBLE  *restrict a  ,DOUBLE *restrict b
                      ,DOUBLE  *restrict w                     
                      ,INT *restrict p      ,INT const nEq);
/*...................................................................*/

/*... decomposicoes*/  
  void fatLDMt(DOUBLE *restrict a,DOUBLE *restrict d,DOUBLE *restrict r
              ,DOUBLE *restrict w,INT const nEq);
  void fatLDLt(DOUBLE *restrict a,DOUBLE *restrict d,DOUBLE *restrict r
              ,INT const nEq);
  void  fatGGt(DOUBLE *restrict a,INT const nEq);
  void   fatLU(DOUBLE *restrict a,INT const nEq,short const code);
  void fatLUpp(DOUBLE *restrict a,INT    *restrict p    
              ,INT const nEq     ,short const code);
/*...................................................................*/
  
/*... decomposicoes imcomplretas*/
  void  ic0   (DOUBLE *restrict a,INT const nEq);
  
/*... solver*/
  void solvTri(DOUBLE *restrict a,DOUBLE *restrict d
              ,DOUBLE *restrict b,INT const nEq
              ,short const code);
  void solvCholesky(DOUBLE *restrict a,DOUBLE *restrict b,INT const nEq);
  void solvLU(DOUBLE *restrict a      ,DOUBLE *restrict b,INT const nEq);
  void solvLUpp(DOUBLE *restrict a    ,DOUBLE *restrict b
               ,DOUBLE *restrict r    ,INT *restrict p 
               ,INT const nEq);
/*...................................................................*/

/*...*/
  void   diagScaling(DOUBLE *restrict a,DOUBLE *scaling
                    ,INT const nEq);
  void scalingSystem(DOUBLE *restrict a      ,DOUBLE *restrict b
                    ,DOUBLE *restrict scaling
                    ,INT const nEq);
/*...................................................................*/

/*...*/
  void solverD(DOUBLE *a,DOUBLE *b,INT const nEq
              ,short const code   ,bool fscaling);
/*...................................................................*/


#endif /*GAUSS_H*/
