#include<Graph.h>
/********************************************************************* 
 * CONVGRAPHCSR: ordena o gafo no formato CSR                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * n      -> numera de linhas                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
/*
void convGraph(INT *xAdj       ,INT *adjncy ,INT const *adj
              ,short const *nViz,short maxViz ,INT numel
              ,bool xAdjFlag    ,bool adjFlag){

  INT nel,kk=0,viz,n;
  int aux = 0;
  
  if(xAdj) xAdj[0] = 0;
  for(nel=0;nel<numel;nel++){
    for(viz=0;viz<nViz[nel];viz++){
      n = MAT2D(nel,viz,adj,maxViz);
      if( n != -1){
        aux++;
        if(adjFlag) 
          adjncy[kk++] = n - 1;  
      }
    } 
    if(xAdjFlag) xAdj[nel+1] = xAdj[nel] + aux;
    aux = 0;
  }

}
*/
/*********************************************************************/ 

/********************************************************************* 
 * SORTGRAPHCSR: ordena o gafo no formato CSR                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * n      -> numera de linhas                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void sortGraphCsr(INT *ia,INT *ja,INT n){

  INT i,nl;
  
  for(i=0;i<n;i++){
    nl = ia[i+1] - ia[i];
    if(nl!=0) bubblesort(&ja[ia[i]],nl);
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * BUBBLESORT :  ordena o arranjo em ordem crescente                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ja     -> arranjo                                                 * 
 * n      -> dimensao do arranjo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> arranjo ordenado                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void bubblesort(INT *ja,INT n){
  
  INT i,j;
  bool itroca;
  
  do{
    itroca=false;
    for(i=1;i<n;i++){
      if(ja[i]   < ja[i-1]){
         j       = ja[i-1];
         ja[i-1] = ja[i];
         ja[i]   = j;
         itroca  = true;
      }
    }
  }while(itroca);

}
/*********************************************************************/ 

