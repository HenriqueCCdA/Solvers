#include<stdio.h>
#include<stdlib.h>
/*...*/
#include<mmio.h>
#include<Solv.h>
#include<Coo.h>
/*...................................................................*/

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

void deflatInit(DOUBLE *zt,INT const nEq,int const nDeflation)
{
  int nDiv  = nEq/nDeflation;
  int resto = nEq%nDeflation;
  int i,j,k;
  
  for(i=0;i<nDeflation-1;i++){
    for(j=0;j<nEq;j++){
      MAT2D(i,j,zt,nEq) = 0.e0;
    }
  }

//printf("%d %d\n",nDiv,resto);
  k = 0;
  for(j=0;j<resto;j++){
    MAT2D(0,j,zt,nEq) = 1.e0;
    k++;
  }
  

  for(i=0;i<nDeflation;i++){
    for(j=0;j<nDiv;j++){
      MAT2D(i,k,zt,nEq) = 1.e0;
      k++;
    }
  }

//  for(i=0;i<nDeflation;i++){
//    for(j=0;j<nEq;j++){
//      printf("%f ",MAT2D(i,j,zt,nEq));
//    }
//    printf("\n");
// }
//  exit(0);
}

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

int main(int argc, char **argv){
    
    MM_typecode matCode;
    int    i;
    int    nEq=0,nnz=0,nLin,nCol;
    int    *ia=NULL,*ja=NULL;
    double *a=NULL,*b=NULL,*b0=NULL;
    double timeIc0   = 0.0e0;
    double timeDpcg1 = 0.0e0;
    double timeDpcg2 = 0.0e0;
    int    sym = 0;
    double *xPcg=NULL,*xIcCg=NULL,*xDpcg1=NULL,*xDpcg2=NULL;
    double *e1=NULL,*e2=NULL,*e3=NULL;
    double normInfxPcg;
/*... PCG*/
    double timePCG   = 0.0e0;
    int    *iaCsrUp = NULL,*jaCsrUp = NULL;
    double *adCsrUp = NULL,*aCsrUp = NULL;
    int    *iaCsrLo = NULL,*jaCsrLo = NULL;
    double *adCsrLo = NULL  ,*aCsrLo = NULL;
    double *r=NULL,*z=NULL,*m=NULL,*f=NULL;
/*... DPCG*/
    double *deflation = NULL,*e=NULL,*wtw=NULL,*az=NULL;
    double *q1=NULL,*q2=NULL,*x0=NULL,*v=NULL;
    int nDeflation= 45;
    double *mIcCsr=NULL,*mIcCsc=NULL,*mIcd=NULL;
    double *w1=NULL;
    int    *ws1=NULL,*ws2=NULL,*ws3=NULL;
    bool   fcg    = true;
    bool   fdcg   = true;
    bool   fic0cg = false;
    int   nadCsr=0;
    int *aux=NULL;
    double tolPcg;
    double tolDpcg; 
    double tolDx=1.0e-06;
    int maxIt=50000;
    FILE *fLog=NULL;
/*..................................................................*/
    FILE   *fileCoo=NULL,*fileB=NULL,*fileOut=NULL;
    FILE   *filePcg=NULL,*fileDpcg1=NULL,*fileDpcg2=NULL;

    if( argc < 6){
      printf("%s: input output nDeflation tolPcg tolDpcg.\n",argv[0]);
      exit(EXIT_FAILURE);
    }

    tolPcg  = atof(argv[4]); 
    tolDpcg = atof(argv[5]); 

/*... abertura dos arquivos*/
    fileCoo = fopen(argv[1],"r");
    if(fileCoo == NULL){ 
      printf("Erro na abertura da arquivo %s.\n",argv[1]);
      exit(EXIT_FAILURE);
    }
    fileB   = fopen(argv[2],"r");
    if(fileB == NULL){ 
      printf("Erro na abertura da arquivo %s.\n",argv[2]);
      exit(EXIT_FAILURE);
    }
    
    nDeflation = atol(argv[3]);


    filePcg   = fopen("pcgLog.txt","w");
    fileDpcg1 = fopen("dpcgLog1.txt","w");
    fileDpcg2 = fopen("dpcgLog2.txt","w");
/*... leitura da matriz no formato COO*/
    mm_read_banner(fileCoo,&matCode);
    mm_read_mtx_crd_size(fileCoo,&nLin,&nCol,&nnz);
        
    nEq = nLin;
    if(mm_is_symmetric(matCode)) sym = 1;
/*...................................................................*/

/*... alocacao da memoria*/
/*... vetor de linhas*/
    ia = (int*) malloc(sizeof(int)*nnz);
    if( ia == NULL ){
      printf("Erro na alocaco do vetor ia.\n");
      exit(EXIT_FAILURE);
    }
/*... vetor de colunas*/
    ja = (int*) malloc(sizeof(int)*nnz);
    if( ja == NULL ){
      printf("Erro na alocaco do vetor ja.\n");
      exit(EXIT_FAILURE);
    }
/*... vetor dos valores dos coeficientes*/
    a  = (double*) malloc(sizeof(double)*nnz);
    if(  a == NULL ){
      printf("Erro na alocaco do vetor  a.\n");
      exit(EXIT_FAILURE);
    }
/*... vetor dos valores dos coeficientes*/
    b  = (double*) malloc(sizeof(double)*nLin);
    if(  b == NULL ){
      printf("Erro na alocaco do vetor  b.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*... vetor dos valores dos coeficientes*/
    b0 = (double*) malloc(sizeof(double)*nLin);
    if( b0 == NULL ){
      printf("Erro na alocaco do vetor  b0.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

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


/*... PCG*/    
/*... alocando vetorores par o PCG*/
    m  = (double*) malloc(sizeof(double)*nEq);
    if(  m == NULL ){
      printf("Erro na alocaco do vetor  m.\n");
      exit(EXIT_FAILURE);
    }
    z  = (double*) malloc(sizeof(double)*nEq);
    if(  z == NULL ){
      printf("Erro na alocaco do vetor  z.\n");
      exit(EXIT_FAILURE);
    }
    r  = (double*) malloc(sizeof(double)*nEq);
    if(  r == NULL ){
      printf("Erro na alocaco do vetor  r.\n");
      exit(EXIT_FAILURE);
    }
    xPcg  = (double*) malloc(sizeof(double)*nEq);
    if(  xPcg == NULL ){
      printf("Erro na alocaco do vetor  xPcg.\n");
      exit(EXIT_FAILURE);
    }
    f  = (double*) malloc(sizeof(double)*nEq);
    if(  f == NULL ){
      printf("Erro na alocaco do vetor  f.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*... CsrD*/
    nadCsr = nnz - nEq;
/* ... csr lower*/
    iaCsrLo  = (int*) malloc(sizeof(int)*(nEq+1));
    if(  iaCsrLo == NULL ){
      printf("Erro na alocaco do vetor  ia.\n");
      exit(EXIT_FAILURE);
    }
    jaCsrLo  = (int*) malloc(sizeof(int)*nadCsr);
    if(  jaCsrLo == NULL ){
      printf("Erro na alocaco do vetor  ja.\n");
      exit(EXIT_FAILURE);
    }
    adCsrLo  = (double*) malloc(sizeof(double)*nEq);
    if(  adCsrLo == NULL ){
      printf("Erro na alocaco do vetor  ad.\n");
      exit(EXIT_FAILURE);
    }
    aCsrLo  = (double*) malloc(sizeof(double)*nadCsr);
    if(  aCsrLo == NULL ){
      printf("Erro na alocaco do vetor  al.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/* ... csr upper*/
    iaCsrUp  = (int*) malloc(sizeof(int)*(nEq+1));
    if(  iaCsrUp == NULL ){
      printf("Erro na alocaco do vetor  ia.\n");
      exit(EXIT_FAILURE);
    }
    jaCsrUp  = (int*) malloc(sizeof(int)*nadCsr);
    if(  jaCsrUp == NULL ){
      printf("Erro na alocaco do vetor  ja.\n");
      exit(EXIT_FAILURE);
    }
    adCsrUp  = (double*) malloc(sizeof(double)*nEq);
    if(  adCsrUp == NULL ){
      printf("Erro na alocaco do vetor  ad.\n");
      exit(EXIT_FAILURE);
    }
    aCsrUp  = (double*) malloc(sizeof(double)*nadCsr);
    if(  aCsrUp == NULL ){
      printf("Erro na alocaco do vetor  al.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/
      
/*... ICCG*/
    mIcd    = (double*) malloc(sizeof(double)*(nEq));
    if( mIcd   == NULL ){
      printf("erro na alocaco do vetor  mic.\n");
      exit(EXIT_FAILURE);
    }

    mIcCsr  = (double*) malloc(sizeof(double)*(nadCsr));
    if( mIcCsr == NULL ){
      printf("erro na alocaco do vetor  mic.\n");
      exit(EXIT_FAILURE);
    }

    mIcCsc  = (double*) malloc(sizeof(double)*(nadCsr));
    if( mIcCsc == NULL ){
      printf("erro na alocaco do vetor  mic.\n");
      exit(EXIT_FAILURE);
    }

    ws1     = (int   *) malloc(sizeof(int   )*nEq   );
    if( ws1 == NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }
    
    ws2     = (int   *) malloc(sizeof(int   )*nEq   );
    if( ws2== NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }
      
    ws3     = (int   *) malloc(sizeof(int   )*nEq   );
    if( ws3== NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }
      
    w1      = (double*) malloc(sizeof(double)*nEq   );
    if( w1 == NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }

    xIcCg = (double*) malloc(sizeof(double)*nEq);
    if(  xIcCg == NULL ){
      printf("Erro na alocaco do vetor  xIcCg.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*... DPCG*/
    deflation = (double*) malloc(sizeof(double)*nEq*nDeflation);
    if( deflation == NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }
    
    az        = (double*) malloc(sizeof(double)*nEq*nDeflation);
    if( az == NULL ){
      printf("Erro na alocaco do vetor  w.\n");
      exit(EXIT_FAILURE);
    }
    
    e         = (double*) malloc(sizeof(double)*nDeflation*nDeflation);
    if( e == NULL ){
      printf("Erro na alocaco do vetor  e.\n");
      exit(EXIT_FAILURE);
    }
    
    wtw       = (double*) malloc(sizeof(double)*nDeflation*nDeflation);
    if( wtw == NULL ){
      printf("Erro na alocaco do vetor  wtw.\n");
      exit(EXIT_FAILURE);
    }
    
    q1        = (double*) malloc(sizeof(double)*nDeflation);
    if( q1 == NULL ){
      printf("Erro na alocaco do vetor  q1.\n");
      exit(EXIT_FAILURE);
    }
    
    q2        = (double*) malloc(sizeof(double)*nDeflation);
    if( q2 == NULL ){
      printf("Erro na alocaco do vetor  q2.\n");
      exit(EXIT_FAILURE);
    }
    
    x0        = (double*) malloc(sizeof(double)*nEq);
    if( x0 == NULL ){
      printf("Erro na alocaco do vetor  x0.\n");
      exit(EXIT_FAILURE);
    }
    
    v         = (double*) malloc(sizeof(double)*nEq);
    if( v  == NULL ){
      printf("Erro na alocaco do vetor  v .\n");
      exit(EXIT_FAILURE);
    }
    
    xDpcg1= (double*) malloc(sizeof(double)*nEq);
    if(  xDpcg1 == NULL ){
      printf("Erro na alocaco do vetor  xDpcg1.\n");
      exit(EXIT_FAILURE);
    }
    
    xDpcg2= (double*) malloc(sizeof(double)*nEq);
    if(  xDpcg2 == NULL ){
      printf("Erro na alocaco do vetor  xDpcg2.\n");
      exit(EXIT_FAILURE);
    }

/*...................................................................*/

/*...................................................................*/

/*...*/
    aux    = (int   *) malloc(sizeof(int)*nEq);
    if(  aux == NULL ){
      printf("Erro na alocaco do vetor  aux.\n");
      exit(EXIT_FAILURE);
    }
    
    e1    = (double*) malloc(sizeof(double)*nEq);
    if(  e1 == NULL ){
      printf("Erro na alocaco do vetor  e1.\n");
      exit(EXIT_FAILURE);
    }
    
    e2    = (double*) malloc(sizeof(double)*nEq);
    if(  e2 == NULL ){
      printf("Erro na alocaco do vetor  e2.\n");
      exit(EXIT_FAILURE);
    }
    
    e3    = (double*) malloc(sizeof(double)*nEq);
    if(  e3 == NULL ){
      printf("Erro na alocaco do vetor  e2.\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*...*/

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
    if(fcg){ 
      alphaProdVector(1.0e0,b ,nEq,f);
      alphaProdVector(1.0e0,b ,nEq,b0);
/*... ad-> m*/
      alphaProdVector(1.0e0,adCsrUp,nEq,m);
//    for(i=0;i<nEq;i++)
//      m[i] = 1.0e0;  
               
/*... pcg*/
     timePCG = getTimeC();
     pcg(nEq    ,nadCsr
//       ,iaCsrUp,jaCsrUp
//       ,aCsrUp ,adCsrUp ,aCsrUp 
         ,iaCsrLo,jaCsrLo
         ,aCsrLo ,adCsrLo ,aCsrLo 
         ,m      ,f ,b0
         ,xPcg   ,v
         ,z      ,r ,tolPcg
         ,maxIt  ,1          
         ,filePcg,1
         ,matVecCsrDSymUpLower,dot);  
      timePCG = getTimeC() - timePCG;
    }
/*
    for(i=0;i<nEq;i++)    
      printf("%9d   %.8f\n",i,x[i]); 
*/
/*...................................................................*/

/*... ic0cg*/
    printf("\n");
    if(fic0cg){
      alphaProdVector(1.0e0,b ,nEq,f);
//    for(i=0;i<nEq;i++)
//      mIcd[i] = 0.0;               
      timeIc0 = getTimeC();
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
         ,f      ,xIcCg
         ,z      ,r
         ,tolPcg
         ,maxIt ,        1          
         ,fLog  ,        0
         ,matVecCsrDSymUpLower,dot);
  
      timeIc0 = getTimeC() -timeIc0;
    }      
/*...*/
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
    if(fdcg) {
      alphaProdVector(1.0e0,adCsrUp,nEq,m);
/*
      for(i=0;i<nEq;i++)
        m[i] = 1.0e0; 
*/           
      alphaProdVector(1.0e0,b ,nEq,f);
      alphaProdVector(1.0e0,b ,nEq,b0);
/*...*/      
      deflatInit(deflation,nEq,nDeflation);
/*...................................................................*/
      timeDpcg1 = getTimeC();
/*  
      dpcg(nEq           ,nadCsr  
           ,iaCsrLo      ,jaCsrLo
           ,aCsrLo       ,adCsrLo ,aCsrLo
           ,m            ,f       ,b0     
           ,xDpcg1
           ,z            ,r       ,v 
           ,tolDpcg 
           ,deflation    
           ,e            ,wtw
           ,az           ,x0             
           ,q1           ,q2             
           ,nDeflation                              
           ,maxIt                   ,1          
           ,fileDpcg1               ,1          
           ,matVecCsrDSymUpLower    ,dot);
      timeDpcg1 = getTimeC() - timeDpcg1;
*/  
           
      alphaProdVector(1.0e0,b ,nEq,f);
/*...*/      
      deflatInit(deflation,nEq,nDeflation);
/*...................................................................*/
      timeDpcg2 = getTimeC() - timeDpcg2;
      dpcg2(nEq            ,nadCsr  
           ,iaCsrLo       ,jaCsrLo
           ,aCsrLo        ,adCsrLo ,aCsrLo
           ,m             ,f       ,b0        
           ,xDpcg2
           ,z             ,r       ,v 
           ,tolDpcg 
           ,deflation     
           ,e             ,wtw
           ,az            ,x0             
           ,q1            ,q2             
           ,nDeflation                              
           ,maxIt                   ,1          
           ,fileDpcg2               ,1          
           ,matVecCsrDSymUpLower    ,dot);
      timeDpcg2 = getTimeC() - timeDpcg2;

    }
    
    printf("TimePcg    %.8f\n",timePCG); 
    printf("TimeIc0    %.8f\n",timeIc0); 
//    printf("TimeDpcg1  %.8f\n",timeDpcg1); 
    printf("TimeDpcg2  %.8f\n",timeDpcg2); 
 
    normInfxPcg = normInf(xPcg,nEq,1);

    for(i=0;i<nEq;i++){
      e1[i] = fabs( (xPcg[i] - xIcCg[i])/xPcg[i]); 
      e2[i] = (xPcg[i] - xDpcg1[i])*(xPcg[i] - xDpcg1[i]); 
      e3[i] = (xPcg[i] - xDpcg2[i])*(xPcg[i] - xDpcg2[i]); 
//      if( fic0cg && e1[i] > tolDx) 
//        printf("%d %.16lf %.16lf %.16lf\n",i,e1[i],xPcg[i],xIcCg[i]);
//      if( fdcg && e2[i] > tolDx) 
//        printf("%d erro %.16lf xPcg %.16e xDpcg %.16e\n"
//              ,i,e2[i],xPcg[i],xDpcg[i]);
    } 
    
//    printf("erro Dpcg1  %.8e\n",sqrt(dot(e2,e2,nEq))/normInfxPcg); 
    printf("erro Dpcg2  %.8e\n",sqrt(dot(e3,e3,nEq))/normInfxPcg); 
    
    return EXIT_SUCCESS; 
}

