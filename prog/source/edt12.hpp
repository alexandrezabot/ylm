//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// 
//  Função que calcula a Transformada Euclidiana de Distâncias
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef EUCLIDIAN_DISTANCE_TRANSFORM
#define EUCLIDIAN_DISTANCE_TRANSFORM


// Padrão do C++
#include <iostream>
using namespace std;


// Minhas bibliotecas
#include "meusTipos.hpp"


// Para matrizes 3D com Boost
#include "boost/multi_array.hpp"
using namespace boost;
typedef multi_array<int , 3> matrix;


// Para operações da Tranformada de Distância Euclidiana
#include "operators10.hpp"










//------------------------------------------------------------------------------
// DESCRICAO:
//   Calcula Transformada de Distância Euclidiana 2D ou 3D
//   Aplico o algoritmo rápido (seção 3.4) de Saito & Toriwaki, 1994
//
//   http://www.sciencedirect.com/science/article/pii/0031320394901333
//
//   Baseado no código sedt.cc de David Coeurjolly
//   http://liris.cnrs.fr/~dcoeurjo
//
// ATENCAO:
//   O resultado é colocado numa matriz passada por referência por questão
//  de eficiência. Senão, eu retornaria, mas isso é lento.
//
// RECEBE:
//   background   => Cor do fundo
//   img          => Imagem para calcular a transformada (sobreescrita com edt)
//   maux         => Matriz auxiliar já inicializada com o tamanho da imagem
template<typename T>
void euclidian_distance_transform( MTci &bg,
                                   multi_array<T , 3> &img,
                                   multi_array<T , 3> &maux ){
  
  // Tamanho da imagem  
  MTci nx = img.shape()[0];
  MTci ny = img.shape()[1];
  MTci nz = img.shape()[2];


  
  // ---------------------------------------------------------------------------
  // First step of  the saito  algorithm
  for( int y=0; y<ny; y++){
  for( int z=0; z<nz; z++){
    img[0][y][z] = (img[0][y][z] == bg)?  0:INFTY;
  
    // Forward scan
    for( int x=1; x<nx; x++){
      img[x][y][z] = (img[x][y][z] == bg)?  0:sum( 1, img[x-1][y][z]); 
    }
  
    //Backward scan
    for(int x=nx-2; x>=0; x--){    
      if( img[x+1][y][z] < img[x][y][z] ) 
        img[x][y][z]=sum(1, img[x+1][y][z]);
    }
  }}




  // ---------------------------------------------------------------------------
  // Second step of the Saito algorithm using the
  //[Meijster/Roerdnik/Hesselink] optimization
  MTvi s(ny); //Center of the upper envelope parabolas
  MTvi t(ny); //Separating index between 2 upper envelope parabolas 
  int q, w;

  for( int x=0; x<nx; x++){
  for( int z=0; z<nz; z++){
    q=0;
    s[0] = 0;
    t[0] = 0;

    //Forward Scan
    for( int u=1; u<ny; u++){
      
      while( (q >= 0) &&
       (F( t[q], s[q], prod(img[x][s[q]][z],img[x][s[q]][z]) ) > 
        F( t[q], u,    prod(img[x][u][z],   img[x][u][z])    ) )
           ){
        q--;
      }
    
      if(q<0){
        q=0;
        s[0]=u;
        
      }else{
        w = 1 + Sep( s[q],
                     u,
                     prod( img[x][s[q]][z],img[x][s[q]][z] ),
                     prod( img[x][u][z]   ,img[x][u][z]    ) );

        if( w < ny ){
          q++;
          s[q]=u;
          t[q]=w;
        }
      }
    }

    //Backward Scan
    for( int u=ny-1; u>=0; --u ){
      maux[x][u][z] = F( u, s[q], prod( img[x][s[q]][z], img[x][s[q]][z]) );	   
      
      if( u==t[q] ) q--;
    }
  }}




  // ---------------------------------------------------------------------------
  // Third step of the Saito algorithm using the
  //[Meijster/Roerdnik/Hesselink] optimization
  s.resize(nz); //Center of the upper envelope parabolas
  t.resize(nz); //Separating index between 2 upper envelope parabolas 


  for( int x=0; x<nx; x++){
  for( int y=0; y<ny; y++){
    q=0;
    s[0] = 0;
    t[0] = 0;

    //Forward Scan
    for( int u=1; u<nz; u++){
      while( (q>=0) &&
             (F(t[q],s[q], maux[x][y][s[q]]) > 
              F(t[q],u,maux[x][y][u]))
           ){
             q--;
      }
      
      if(q<0){
        q=0;
        s[0]=u;
        
      }else{
        w = 1 + Sep( s[q],
                      u,
                      maux[x][y][s[q]],
                      maux[x][y][u] );
  
        if( w<nz ){
          q++;
          s[q]=u;
          t[q]=w;
        }
      }
    }
  
    //Backward Scan
    for( int u=nz-1; u>=0; --u){

      // Alterei o código aqui para ficar mais rápido
      // No original, havia a matriz sdt_xy
      // Mas não preciso dela, posso colocar direto no _img      
      img[x][y][u] = F( u, s[q], maux[x][y][s[q]] );	      

      if( u==t[q] ) 
        q--;
    }
  }}

}







#endif // EUCLIDIAN_DISTANCE_TRANSFORM





















