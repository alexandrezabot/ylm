//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Funções para ajudar na manipulação de arquivos
//________________________________________________________
//A.Z. - 03/05 => Criacao
//       10/05 => Uso das bibliotecas aborta.h e meusTipos.h
//       01/07 => Funcao outputFileHead
//       03/08 => Funcao outputFileHead imprime data tambem
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef FILE_UTI_H
#define FILE_UTI_H

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "aborta.hpp"
#include "meusTipos.hpp"
using namespace std;


void mymkdir( MTcs & );
template < class T > inline void abriu( T &, MTcs & );
void substitui_caminho( string & );
void outputFileHead( MTci &, char **, ofstream &, MTcvs & );




//------------------------------------------------------------------------------
// DESCRICAO:
//   Cria um diretório
// RECEBE:
//   folder => Nome do diretório
void mymkdir( MTcs &folder ){
  MTcs mkdir( "mkdir -p " + folder );
  int lixo = system( mkdir.c_str() );
  lixo=lixo;
}




//------------------------------------------------------------------------------
// DESCRICAO:
// Verifica se um arquivo foi aberto corretamente
template < class T >
inline void abriu( T &obj, MTcs &nome ){
  if( !obj.is_open() )
    aborta("Nao foi possivel abrir/criar o arquivo "+nome);
}




//------------------------------------------------------------------------------
// DESCRICAO:
// Recebe uma string que representa o nome de um arquivo e quebra ela para
//  descobrir qual eh o caminho ate o arquivo.
// Atencao: A string eh alterada.
//
// Ex.:   $PWD/file_uti.h  -->  /home/vela/zabot/file_uti.h
//
void substitui_caminho( string &arq ){

  // Pega o inicio e o fim da variavel de ambiente
  int inicio  = arq.find("$");
  int fim     = arq.find("/");
  if(fim<0 && inicio>=0) // Se nao achou "/" na string mas achou $
    fim = arq.size();
  int tamanho = fim - inicio - 1;

  // Pega o nome da variavel de ambiente se ela esta na string e substitui
  if( inicio>=0 ){

    // Descobre o que ha na variavel de ambiente
    string variavel          = string( arq, inicio+1, tamanho );
    char *GETENV = getenv(variavel.c_str());

    if(GETENV==NULL){   // Se nao existe a variavel de ambiente
      string msg  = "Erro na funcao \"substitui_caminho(string &)\" da ";
      msg        += "biblioteca file_util.h\nNao foi possivel pegar a variavel";
      msg        += " de ambiente \"" + variavel + "\"";
      aborta(msg);
    }else               // Senao, substitui o valor na string original
      arq.replace(inicio,tamanho+1,string(GETENV));
  }
}



//------------------------------------------------------------------------------
// DESCRICAO:
//   Gera um cabecalho no arquivo de saida.
//   Este cabecalho eh util para facilitar a identificacao de arquivos de dados
//  gerados por programas.
// RECEBE:
// RETORNA:
void outputFileHead( MTci &argc, char *argv[], ofstream &OUT, MTcvs &vec,
                     MTcvs &comment){
  char *dir = getenv("PWD");

  OUT<< "# ##################################################################\n"
     << "# Arquivo criado com o comando:  ";
  for( int i=0;i<argc;i++ ) OUT << argv[i] << " ";

  OUT<< "\n# No diretorio:  ";
  if( dir != NULL ) OUT << dir;
  else              OUT << "(Nao foi possivel identificar o diretorio!)";

  OUT<< "\n# Em:  " << __TIME__ << "  " << __DATE__  << endl;
  //OUT<< "\n# Tempo de Execução:  "
     //<< clock()/static_cast<double>(CLOCKS_PER_SEC) << " s\n";

  OUT<< "#\n";
  if( !comment.empty() ){
    for( unsigned int i=0;i<comment.size();i++ )
      OUT << "#   * " << comment[i] << "\n";
    OUT<< "#\n";
  }

  for( unsigned int i=0;i<vec.size();i++ )
    OUT << "# Coluna " << setw(2) << i+1 << ": " << vec[i] << "\n";
  OUT<< "# ##################################################################"
     << endl;
}



#endif /* FILE_UTI_H */
