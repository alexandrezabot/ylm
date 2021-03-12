//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Função que faz o Component Labeling de He, Chao & Suzuki, 2011.
//
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef COMPONENT_LABELING
#define COMPONENT_LABELING


#include <iostream>
using namespace std;


#include "meusTipos.hpp"


// Para matrizes 3D com Boost
#include "boost/multi_array.hpp"
using namespace boost;
typedef multi_array<int , 3> matrix;




void component_labeling( matrix &, MTci &, MTci &, MTvi &, 
                         MTvi &, MTvi &, MTci & );
void merge( MTci &, MTci &, MTvi &, MTvi &, MTvi & );
void resolve( MTci &, MTci &, MTvi &, MTvi &, MTvi & );






//------------------------------------------------------------------------------
// DESCRICAO:
//   Identifica as componentes desconexas de uma imagem através do algoritmo
//  de He, Chao & Suzuki, 2011.
//
// ATENCAO:
//   IMG é alterada para conter os labels das regiões desconexas
//
// RECEBE:
//   img          => Imagem para encontrar as regiões desconexas
//   F, B         => Número para Foreground e Background
//   next, tail, rtable => Vetores para trabalhar informações das equivalências
void component_labeling( matrix &IMG, MTci &F, MTci &B, MTvi &next, 
                         MTvi &tail, MTvi &rtable ){


  // Tamanho da imagem  
  MTci nx = IMG.shape()[0];
  MTci ny = IMG.shape()[1];
  MTci nz = IMG.shape()[2];
 
  // Auxiliares
  int lx=0, nl=1;
  MTvi uniq_labels(3);
  int nuniq;
  
  
  for( int x=0; x<nx; x++ ){
  for( int y=0; y<ny; y++ ){
  for( int z=0; z<nz; z++ ){
    if( IMG[x][y][z]==F ){
      
      
      // Pega o label dos vizinhos se eles estão dentro da imagem.
      // Senão, diz que eles são background
      MTci lq = (x>0)?  IMG[x-1][y  ][z  ]:B;
      MTci lp = (y>0)?  IMG[x  ][y-1][z  ]:B;
      MTci lz = (z>0)?  IMG[x  ][y  ][z-1]:B;


      // Cria um vetor com os labels únicos diferentes de background
      // Poderia montar um esquema mais simples, mas quero aproveitar o vetor
      // Senão gasto muito tempo criando e apagando um vetor sempre do mesmo
      // tamanho!
      nuniq=0;
      if( lp!=B ){
        uniq_labels[nuniq] = lp;
        nuniq++;
      }
      if( lq!=B && lq!=lp ){
        uniq_labels[nuniq] = lq;
        nuniq++;
      }
      if( lz!=B && lz!=lp && lz!=lq ){
        uniq_labels[nuniq] = lz;
        nuniq++;
      }


      // Registra equivalências
      switch( nuniq ){
    
        // ---------------------------------------------------------------------
        // Caso 1: Todos os vizinhos são background
        // Cria-se um novo label
        case 0:
          nl++;
          lx = nl;
          
          // He, Chao, Suzuki, 2008, pág 752
          rtable[nl] = nl;
          next[nl]   = -1;
          tail[nl]   = nl;
          break;
      
      
        // ---------------------------------------------------------------------
        // Caso 2: Só há um único label entre os vizinhos
        // Pega este label
        case 1:
          lx=uniq_labels[0];
          break;

      
        // ---------------------------------------------------------------------
        // Caso 3: Há dois labels diferentes entre os vizinhos
        case 2:
          
          // Registra equivalência
          //dual[0] = uniq_labels[0];
          //dual[1] = uniq_labels[1];
          //eq.push_back( dual );
          resolve( uniq_labels[0], uniq_labels[1], next, tail, rtable );
          
          lx = min( uniq_labels[0], uniq_labels[1] );
          break;

      
        // ---------------------------------------------------------------------
        // Caso 4: Há três labels diferentes entre os vizinhos
        case 3:
          
          // Registra equivalência 1
          //dual[0] = uniq_labels[0];
          //dual[1] = uniq_labels[1];
          //eq.push_back( dual );
          
          // Registra equivalência 2
          //dual[0] = uniq_labels[1];
          //dual[1] = uniq_labels[2];
          //eq.push_back( dual );
          resolve( uniq_labels[0], uniq_labels[1], next, tail, rtable );
          resolve( uniq_labels[0], uniq_labels[2], next, tail, rtable );
          resolve( uniq_labels[1], uniq_labels[2], next, tail, rtable );
          
          lx = min( uniq_labels[0], min( uniq_labels[1], uniq_labels[2] ) );
          break;
      }
      
      IMG[x][y][z] = lx;
      
            
    }
  }}}




  // Troco os índices velhos pelos equivalentes
  for( int x=0; x<nx; x++ ){
  for( int y=0; y<ny; y++ ){
  for( int z=0; z<nz; z++ ){
    if( IMG[x][y][z]!=B )
      IMG[x][y][z] = rtable[ IMG[x][y][z] ];
  }}}


}



// He, Chao, Suzuki, 2008, pág 752
void merge( MTci &u, MTci &v, MTvi &next, MTvi &tail, MTvi &rtable ){
  for( int i=v; i!=-1;  ){
    rtable[i] = u;
    i = next[i];
  }
  next[ tail[u] ] = v;
  tail[u] = tail[v];
}

// He, Chao, Suzuki, 2008, pág 752
void resolve( MTci &x, MTci &y, MTvi &next, MTvi &tail, MTvi &rtable ){
  MTci u = rtable[x];
  MTci v = rtable[y];
  if     ( u<v ) merge( u, v, next, tail, rtable );
  else if( v<u ) merge( v, u, next, tail, rtable );
}


#endif // COMPONENT_LABELING





















