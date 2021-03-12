//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// OBJETIVO:
//   Young-Laplace Method (scientific).
//________________________________________________________
// RECEBE:
//   cfg => Configuration file
//________________________________________________________
// MODO DE USAR:
//   ylm config.txt
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#include <iomanip>
#include <iostream>
using namespace std;

#include "aborta.hpp"
#include "meusTipos.hpp"

#include "ylm.hpp"








int main( int argc, char *argv[] ){

  if( argc!=2 ) aborta("Número de parâmetros errado.");
  MTcs config( argv[1] );


  Young_Laplace_Method YLM( config, argc, argv );
  YLM.run();
  

  return 0;
}
