#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
/*...*/
#include<mmio.h>
#include<Solv.h>
#include<Coo.h>
/*...................................................................*/
#include<MUMPS/dmumps_c.h>
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] 

void icZeroCsr(DOUBLE *restrict m, DOUBLE *restrict ad
              ,DOUBLE *restrict aLo
              ,DOUBLE *restrict aUp 
              ,DOUBLE *restrict w   
              ,INT *restrict iaUp,INT *restrict jaUp
              ,INT *restrict iaLo,INT *restrict jaLo
              ,INT    *restrict ws1
              ,INT    *restrict ws2
              ,INT    *restrict ws3
              ,INT const neq,INT const nad);

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


int readB(double *b,FILE *fileIn,int const n){
  int i;
  int erro;
  double value;
  for(i=0;i<n;i++){
    erro = fscanf(fileIn,"%lf",&value);
    if(erro != 1) return -1;
    b[i] = value;
  }
  return 0; 
} 

void matVecCooSym(int neq  ,int nnz ,int *ia     ,int *ja
                 ,double *a,double *dum, double *x
                 ,double *y){
  
  int i,nLin,nCol;  
  double value;
  

  for(i=0;i<neq;i++)
    y[i] = 0.0e0;

  
  for(i=0;i<nnz;i++){
    nLin     = ia[i] - 1;
    nCol     = ja[i] - 1;
    value    = a[i];
    if( nLin == nCol )
      y[nLin] += value*x[nCol];
    else{
      y[nLin] += value*x[nCol];
      y[nCol] += value*x[nLin];
    }
  }


}

int main(int argc, char **argv){
    
    DMUMPS_STRUC_C id;
    MM_typecode matCode;
    int    i;
    int    nEq=0,nnz=0,nLin,nCol;
    int    sym=0;
    int    *ia=NULL,*ja=NULL;
    double *a=NULL,*b=NULL;
    double timeMUMPS = 0.0e0;
    double timeIc0   = 0.0e0;
    double *e=NULL;
/*... PCG*/
    double timePCG   = 0.0e0;
    int    *iaCsrUp = NULL,*jaCsrUp = NULL;
    double *adCsrUp = NULL,*aCsrUp = NULL;
    int    *iaCsrLo = NULL,*jaCsrLo = NULL;
    double *adCsrLo = NULL  ,*aCsrLo = NULL;
    double *r=NULL,*z=NULL,*x=NULL,*m=NULL,*f=NULL;
    double *mIcCsr=NULL,*mIcCsc=NULL,*mIcd=NULL;
    double *w1=NULL;
    int    *ws1=NULL,*ws2=NULL,*ws3=NULL;
    bool   cg = false;
    int   nadCsr=0;
    int *aux=NULL;
    double tol=1.5e-16;
    int maxIt=1500000;
    FILE *fLog=NULL;
/*..................................................................*/
    FILE   *fileCoo=NULL,*fileB=NULL,*fileOut=NULL;
    int myId,ierr=0;
    short erro=0,TotalError=0;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myId);

/*... abertura dos arquivos*/
    if(!myId){
        fileCoo = fopen(argv[1],"r");
        if(fileCoo == NULL){ 
            printf("Erro na abertura da arquivo %s.\n",argv[1]);
            erro=1;
        }
        fileB   = fopen(argv[2],"r");
        if(fileB == NULL){ 
            printf("Erro na abertura da arquivo %s.\n",argv[2]);
            erro=1;
        }
    }
/*... verificacao de possivel erro na alocacao da memoria*/
    ierr = MPI_Allreduce(&erro    , &TotalError, 1,
                         MPI_SHORT, MPI_SUM    , MPI_COMM_WORLD);
    
    if(TotalError){
      ierr = MPI_Finalize();
      exit(EXIT_FAILURE);
    }
/*...................................................................*/
    
    if(!myId){
/*... leitura da matriz no formato COO*/
      mm_read_banner(fileCoo,&matCode);
      mm_read_mtx_crd_size(fileCoo,&nLin,&nCol,&nnz);
        
      nEq = nLin;
      if(mm_is_symmetric(matCode)) sym = 2;
/*...................................................................*/

/*... alocacao da memoria*/
/*... vetor de linhas*/
      ia = (int*) malloc(sizeof(int)*nnz);
      if( ia == NULL ){
        printf("Erro na alocaco do vetor ia.\n");
        erro=1;
      }
/*... vetor de colunas*/
      ja = (int*) malloc(sizeof(int)*nnz);
      if( ja == NULL ){
        printf("Erro na alocaco do vetor ja.\n");
        erro=1;
      }
/*... vetor dos valores dos coeficientes*/
      a  = (double*) malloc(sizeof(double)*nnz);
      if(  a == NULL ){
        printf("Erro na alocaco do vetor  a.\n");
        erro=1;
      }
/*... vetor dos valores dos coeficientes*/
      b  = (double*) malloc(sizeof(double)*nLin);
      if(  b == NULL ){
        printf("Erro na alocaco do vetor  b.\n");
        erro=1;
      }
    }
/*...................................................................*/

/*... verificacao de possivel erro na alocacao da memoria*/
    MPI_Allreduce(&erro    , &TotalError, 1,
                  MPI_SHORT, MPI_SUM    , MPI_COMM_WORLD);
    
    if(TotalError){
      ierr = MPI_Finalize();
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*...*/   
     MPI_Bcast(&sym, 1, MPI_INT, 0,MPI_COMM_WORLD );
/*...................................................................*/

    if(!myId){
/*... leitura da matriz*/
        mm_read_mtx_crd_data(fileCoo,nLin,nCol,nnz,ia,ja,a,matCode);
/*...................................................................*/
        
/*... leitura da matriz no formato COO*/
        mm_read_mtx_array_size(fileB,&nLin,&nCol);

/*... leitura do vetor de forcas */
        readB(b,fileB,nLin);
/*...................................................................*/
        
/*...*/
        fclose(fileCoo);
        fclose(fileB  );
/*...................................................................*/
    }
/*...................................................................*/


/*... PCG*/    
    if(!myId){
/*... alocando vetorores par o PCG*/
      m  = (double*) malloc(sizeof(double)*nEq);
      if(  m == NULL ){
          printf("Erro na alocaco do vetor  m.\n");
          erro=1;
      }
      z  = (double*) malloc(sizeof(double)*nEq);
      if(  z == NULL ){
          printf("Erro na alocaco do vetor  z.\n");
          erro=1;
      }
      r  = (double*) malloc(sizeof(double)*nEq);
      if(  r == NULL ){
          printf("Erro na alocaco do vetor  r.\n");
          erro=1;
      }
      x  = (double*) malloc(sizeof(double)*nEq);
      if(  x == NULL ){
          printf("Erro na alocaco do vetor  x.\n");
          erro=1;
      }
      f  = (double*) malloc(sizeof(double)*nEq);
      if(  f == NULL ){
          printf("Erro na alocaco do vetor  f.\n");
          erro=1;
      }
/*...................................................................*/

/*... CsrD*/
      nadCsr = nnz - nEq;
/* ... csr lower*/
      iaCsrLo  = (int*) malloc(sizeof(int)*(nEq+1));
      if(  iaCsrLo == NULL ){
        printf("Erro na alocaco do vetor  ia.\n");
        erro=1;
      }
      jaCsrLo  = (int*) malloc(sizeof(int)*nadCsr);
      if(  jaCsrLo == NULL ){
        printf("Erro na alocaco do vetor  ja.\n");
        erro=1;
      }
      adCsrLo  = (double*) malloc(sizeof(double)*nEq);
      if(  adCsrLo == NULL ){
        printf("Erro na alocaco do vetor  ad.\n");
        erro=1;
      }
      aCsrLo  = (double*) malloc(sizeof(double)*nadCsr);
      if(  aCsrLo == NULL ){
        printf("Erro na alocaco do vetor  al.\n");
        erro=1;
      }
/*...................................................................*/

/* ... csr upper*/
      iaCsrUp  = (int*) malloc(sizeof(int)*(nEq+1));
      if(  iaCsrUp == NULL ){
        printf("Erro na alocaco do vetor  ia.\n");
        erro=1;
      }
      jaCsrUp  = (int*) malloc(sizeof(int)*nadCsr);
      if(  jaCsrUp == NULL ){
        printf("Erro na alocaco do vetor  ja.\n");
        erro=1;
      }
      adCsrUp  = (double*) malloc(sizeof(double)*nEq);
      if(  adCsrUp == NULL ){
        printf("Erro na alocaco do vetor  ad.\n");
        erro=1;
      }
      aCsrUp  = (double*) malloc(sizeof(double)*nadCsr);
      if(  aCsrUp == NULL ){
        printf("Erro na alocaco do vetor  al.\n");
        erro=1;
      }
/*...................................................................*/
      
/*...*/
      mIcd    = (double*) malloc(sizeof(double)*(nEq));
      if( mIcd   == NULL ){
        printf("erro na alocaco do vetor  mic.\n");
        erro=1;
      }

      mIcCsr  = (double*) malloc(sizeof(double)*(nadCsr));
      if( mIcCsr == NULL ){
        printf("erro na alocaco do vetor  mic.\n");
        erro=1;
      }

      mIcCsc  = (double*) malloc(sizeof(double)*(nadCsr));
      if( mIcCsc == NULL ){
        printf("erro na alocaco do vetor  mic.\n");
        erro=1;
      }

      ws1     = (int   *) malloc(sizeof(int   )*nEq   );
      if( ws1 == NULL ){
        printf("Erro na alocaco do vetor  w.\n");
        erro=1;
      }
      ws2     = (int   *) malloc(sizeof(int   )*nEq   );
      if( ws2== NULL ){
        printf("Erro na alocaco do vetor  w.\n");
        erro=1;
      }
      
      ws3     = (int   *) malloc(sizeof(int   )*nEq   );
      if( ws3== NULL ){
        printf("Erro na alocaco do vetor  w.\n");
        erro=1;
      }
      
      w1      = (double*) malloc(sizeof(double)*nEq   );
      if( w1 == NULL ){
        printf("Erro na alocaco do vetor  w.\n");
        erro=1;
      }
      
/*...................................................................*/
/*...................................................................*/

/*...................................................................*/

/*...*/
      aux    = (int   *) malloc(sizeof(int)*nEq);
      if(  aux == NULL ){
        printf("Erro na alocaco do vetor  aux.\n");
        erro=1;
      }
/*...................................................................*/
    }

/*... verificacao de possivel erro na alocacao da memoria*/
    ierr = MPI_Allreduce(&erro    , &TotalError, 1,
                         MPI_SHORT, MPI_SUM    , MPI_COMM_WORLD);
    
    if(TotalError){
      ierr = MPI_Finalize();
      exit(EXIT_FAILURE);
    }
/*...................................................................*/
    
/*...*/
    if(!myId && cg){

/*... converter do coo -> CsrD - parte inferior*/
      cooToCsr(ia        ,ja    ,a
              ,iaCsrLo   ,jaCsrLo
              ,aCsrLo    ,adCsrLo,aCsrLo
              ,nEq       ,nnz   ,CSRD
              ,aux      
              ,false     ,false ,true);

/*... converter do coo -> CsrD - parte superior*/
      cooToCsr(ia        ,ja    ,a
              ,iaCsrUp   ,jaCsrUp
              ,aCsrUp    ,adCsrUp,aCsrUp
              ,nEq       ,nnz   ,CSRD
              ,aux      
              ,true      ,false ,false);
      free(aux);
/*...................................................................*/

/*... b -> f*/
      alphaProdVector(1.0e0,b ,nEq,f);
/*... ad-> m*/
      alphaProdVector(1.0e0,adCsrUp,nEq,m);
//    for(i=0;i<nEq;i++)
//      m[i] = 0.0;               
/*... pcg*/
      timePCG = MPI_Wtime();
      pcg(nEq    ,nadCsr
//       ,iaCsrUp,jaCsrUp
//       ,aCsrUp ,adCsrUp ,aCsrUp 
         ,iaCsrLo,jaCsrLo
         ,aCsrLo ,adCsrLo ,aCsrLo 
         ,m      ,f ,x
         ,z      ,r ,tol
         ,maxIt  ,1          
         ,fLog   ,0
         ,matVecCsrDSymUpLower,dot);  
      timePCG = MPI_Wtime() - timePCG;
/*
      for(i=0;i<nEq;i++)    
        printf("%9d   %.8f\n",i,x[i]); 
*/ 
/*...................................................................*/

/*... ic0cg*/
      printf("\n");
      alphaProdVector(1.0e0,b ,nEq,f);
//    for(i=0;i<nEq;i++)
//      mIcd[i] = 0.0;               
      timeIc0 = MPI_Wtime();
      ic0Csr(mIcd
            ,mIcCsr
            ,mIcCsc
            ,adCsrUp
            ,aCsrUp
            ,w1     
            ,iaCsrUp
            ,jaCsrUp
            ,iaCsrLo
            ,jaCsrLo
            ,ws1
            ,ws2
            ,ws3
            ,nEq,nadCsr);
  
/*    printf("Md\n");
      for(i=0;i<nEq;i++)    
        printf("%9d   %20.6f\n",i,mIcd[i]);    
      printf("Csr\n");
      for(i=0;i<nadCsr;i++)    
        printf("%9d   %20.6f\n",i,mIcCsr[i]);    
      printf("Csc\n");
      for(i=0;i<nadCsr;i++)    
        printf("%9d   %20.6f\n",i,mIcCsc[i]);
      exit(0);*/
      ic0cg(nEq      ,nadCsr  
          ,iaCsrLo  ,jaCsrLo
          ,aCsrLo  ,adCsrLo,aCsrLo
          ,mIcd 
          ,mIcCsr
          ,mIcCsc
          ,iaCsrUp,jaCsrUp
          ,iaCsrLo,jaCsrLo
          ,f      ,x
          ,z      ,r
          ,tol
          ,maxIt ,        1          
          ,fLog  ,        0
          ,matVecCsrDSymUpLower,dot);
      timeIc0 = MPI_Wtime() - timeIc0;
      
/*...*/
      printf("timePCG %lf\ntimeICCG %lf\n",timePCG,timeIc0);
/*
      printf("%lf %lf\n",timePCG,timeIc0);
      printf("Md\n");
      for(i=0;i<nEq;i++)    
        printf("%9d   %20.6f\n",i,mIcd[i]);    
      printf("Csr\n");
      for(i=0;i<nadCsr;i++)    
        printf("%9d   %20.6f\n",i,mIcCsr[i]);    
      printf("Csc\n");
      for(i=0;i<nadCsr;i++)    
        printf("%9d   %20.6f\n",i,mIcCsc[i]);
*/
/*
      for(i=0;i<nEq;i++)    
        printf("%9d   %.8f\n",i,x[i]);    
*/

    
//    fileOut = fopen(argv[4],"w");
//    for(i=0;i<nEq;i++)    
//      fprintf(fileOut,"%9d   %20.6e\n",i+1,x[i]);    
//    fclose(fileOut);
/*...*/
      
      free(iaCsrLo);
      free(jaCsrLo);
      free(adCsrLo);
      free(aCsrLo);
      free(m);
      free(r);
      free(f);
/*...................................................................*/
    }
/*...................................................................*/

/*...  */
    id.job = JOB_INIT;
    id.sym = sym;
    id.comm_fortran=USE_COMM_WORLD;
    dmumps_c(&id);
/*...................................................................*/

/*... definindo o problema no host*/
    if(!myId){
      id.n   = nEq;
      id.nz  = nnz;
      id.irn = ia;
      id.jcn = ja;
      id.a   = a;  
      id.rhs = b;  
    }
/*...................................................................*/

/*... output do MUMPS*/
    id.ICNTL(1) = 0;
    id.ICNTL(2) = 6;
    id.ICNTL(3) = 6;
    id.ICNTL(4) = 0;
/*...................................................................*/

/*... */
//  id.ICNTL(7) = 5;
//  id.ICNTL(8)  = 77;
//  id.ICNTL(10) = 0;
//  id.ICNTL(11) = 0;
/*...................................................................*/
/*... resolvendo Ax=b*/
    timeMUMPS = MPI_Wtime();
    id.job=6;
    dmumps_c(&id);
    timeMUMPS = MPI_Wtime() - timeMUMPS;
/*...................................................................*/

/*... finalizando MUMPS*/
    id.job=JOB_END; 
    dmumps_c(&id); 
/*...................................................................*/

/*...*/         
/*  if(!myId){
      fileOut = fopen(argv[3],"w");
      for(i=0;i<nEq;i++)    
        fprintf(fileOut,"%9d   %.14e\n",i+1,b[i]);    
      fclose(fileOut);
    }*/
/*...................................................................*/
    exit(0);

/*...*/
    if(!myId){
      e  = (double*) malloc(sizeof(double)*nEq);
      if(  b == NULL ){
        printf("Erro na alocaco do vetor  e.\n");
        erro=1;
      }
      for(i=0;i<nEq;i++){
         e[i] = x[i] - b[i];     
      }    
      printf("Time MUMPS: %16.6lf\n"
             "Time PCG  : %16.6lf\n"
             ,timeMUMPS,timePCG);    
      printf("Erro relativo: %.14e\n"
             ,dot(e,e,nEq)/dot(b,b,nEq));    
    }
/*...................................................................*/
    ierr = MPI_Finalize();


    return EXIT_SUCCESS; 
}

