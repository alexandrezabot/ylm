#ifndef YOUNG_LAPLACE_METHOD_HPP
#define YOUNG_LAPLACE_METHOD_HPP


#include <fstream>
#include <iomanip>
#include <stdint.h>
using namespace std;



// YLM auxiliary libraries
#include "edt12.hpp"
#include "component_labeling10.hpp"


// Libraries for general use
#include "file_uti.hpp"
#include "meusTipos.hpp"
#include "dataFileReader.hpp"




// -----------------------------------------------------------------------------
// Debugging software options

// Uncomment this to get high performance
#define BOOST_DISABLE_ASSERTS

// When testing the software, it is usefull to print intermediary set images
MTcb print_set_images=false;
// -----------------------------------------------------------------------------



#include "boost/format.hpp"




// 3D Boost matrices
#include "boost/multi_array.hpp"
using namespace boost;
typedef multi_array<int , 3> matrix;







// -----------------------------------------------------------------------------
class Young_Laplace_Method{
public:

  Young_Laplace_Method( MTcs &, MTci &, char *[] );
  void run( void );

  
private:
  int _F, _B; // Foreground and background values


  int _ny, _nx, _nz;  // Rock size
  double _n;
  matrix _mm;         // Rock
  matrix _matrix1, _matrix2;  // Auxiliary matrices
  
  MTvi _d;            // Step diameters, only odd numbers

  double _cost2;   // Cos (contact angle) ^2
                   //   90º: 0 => wet
                   //    0º: 1 => not wet
                   
  bool _wet, _compressible;   // Wettability and compressibility
  

  // Folder and file names
  string _outImgRoot, _whichImg, _outImgDir; 
  ofstream _dat;
  
  
  // Voxel values
  int _I,   // Inlet fluid
      _O,   // Outlet fluid
      _S;   // Solid
    
  
  // Chamber  
  int _source[3] = {0,0,0}; // Voxel of inlet  fluid source
  int _sink[3]   = {0,0,0}; // Voxel of outlet fluid sink

  
  // Vector for the Component Labeling function
  // They are created here to avoid initializatino overhead inside function   
  MTvi _next, _tail, _rtable;

  
  // Euclidian Distance Transform (EDT) of the original image.
  matrix _edt;
  
  
  // Auxiliary function
  template<typename T>
  void _writevtk( MTcs &, const multi_array<T , 3> & );
};







//------------------------------------------------------------------------------
// DESCRIPTION:
//   Constructor
//
// INPUT:
//   config => Configuration file
Young_Laplace_Method::Young_Laplace_Method( MTcs &config, MTci &argc, 
  char *argv[] ){
    
  
  int iaux;
  string saux;
  _B = 0;
  _F = 1;



  DataFileReader dfr(config);
  ++dfr;  dfr>>saux;
  MTcs DigRockFileName(saux);
  
  ++dfr;  dfr>>_outImgRoot;  // Output filenames
  ++dfr;  dfr>>_outImgDir;   // Folder to the output
  
  ++dfr;  dfr>>_whichImg;


  ++dfr;   dfr >> _cost2;
  if( _cost2<0 || _cost2>60 ) aborta("Wrong value for contact angle!");
  _cost2 = pow( cos( _cost2* M_PI/180.0 ), 2 );

  ++dfr;   dfr >> _wet;
  ++dfr;   dfr >> _compressible;
  _compressible = !_compressible;
  

  // Make folder
  _outImgDir += "/";
  mymkdir( _outImgDir );


  ++dfr;  dfr >> _I;
  ++dfr;  dfr >> _O;
  _S=0;
  while( _S==_I || _S==_O ) _S+=1;  
  
  
  // ---------------------------------------------------------------------------
  // Read Digital Rock (vtk)
  DataFileReader digRockFile( DigRockFileName );

  ++digRockFile;
  ++digRockFile; // ASCII
  ++digRockFile; // DATASET STRUCTURED_POINTS
  ++digRockFile; // DIMENSIONS 
  digRockFile >> saux >> _nx >> _ny >> _nz;
  ++digRockFile; // ASPECT_RATIO 1 1 1
  ++digRockFile; // ORIGIN 0 0 0
  ++digRockFile; // POINT_DATA 3330000
  ++digRockFile; // SCALARS geometria int
  ++digRockFile; // LOOKUP_TABLE default
  ++digRockFile; // next line


  // Allocate memory
  matrix::extent_gen extents;
  _mm.resize      ( extents[_nx][_ny][_nz] );
  _edt.resize     ( extents[_nx][_ny][_nz] );
  _matrix1.resize ( extents[_nx][_ny][_nz] );
  _matrix2.resize ( extents[_nx][_ny][_nz] );


  for( int z=0; z<_nz; z++ ){
  for( int y=0; y<_ny; y++ ){
  for( int x=0; x<_nx; x++ ){
    digRockFile >> iaux;
    if( iaux!=_I  &&  iaux!=_O )
      iaux=_S;
    
    _edt[x][y][z] = _mm[x][y][z] = iaux;
  }
  ++digRockFile;
  }}
  digRockFile.close();
  // ---------------------------------------------------------------------------



  // Source and sink of fluids
  ++dfr;  dfr >> _source[0] >> _source[1] >> _source[2];
  ++dfr;  dfr >> _sink[0]   >> _sink[1]   >> _sink[2];
  
  // Python alike: if output fluid voxel index is negative, sum dim length.
  // vec[-1] = vec[n-1]
  if( _sink[0]<0 ) _sink[0] += _nx;
  if( _sink[1]<0 ) _sink[1] += _ny;
  if( _sink[2]<0 ) _sink[2] += _nz;





  // ---------------------------------------------------------------------------
  // Read and evaluate the spheres diameters instructions
  //
  // Configuration file may contain as many lines as the user desires with
  //   dmin dmax dstep
  // whre:
  //   dmin  => Minor diameter
  //   dmax  => Major diameter
  //   dstep => Step in diameter
  //
  // Example:
  //   1 7 2 ->  1,3,5,7
  //   1 4 2 ->  1,3
  //
  // I unify the diameters, eliminating the repeated and the even values.
  // The order of the final vector depends on the inlet fluid wetabillity:
  //   Not wet: increasing order
  //   Wet    : descreasing order

  int dmin, dmax, dstep;
  while( !++dfr ){
    dfr >> dmin >> dmax >> dstep;
    for( int dd=dmin; dd<=dmax; dd+=dstep ){
      if( dd%2==0 ) _d.push_back(dd-1);
      else          _d.push_back(dd);
    }
  }
  dfr.close();     
  
  sort( _d.begin(), _d.end() );
  
  MTvi::iterator it = unique( _d.begin(), _d.end() );
  _d.resize( distance( _d.begin(),it ) );

  if( _d.size()==0 )
  aborta("I couldn't initialize the list of diameters. Check the values and if the config file is according to the template file config.");

  if( !_wet )  reverse( _d.begin(), _d.end() );
  // ---------------------------------------------------------------------------
  
  
  



  // Auxiliary vectors to the Component Labeling Algorithm
  // They are created here to avoid allocation overhead if they were 
  // allocated inside the function
  // For the first neighbours connection, the maximum number of disconnected
  // voxels is half of the total voxels. Think in a 3d chess board.
  _n = static_cast<double>( _nx*_ny*_nz );
  MTci max_eq = static_cast<int> ( ceil( _n/2.0 ) );
  _tail.resize( max_eq );
  _next.resize( max_eq );
  _rtable.resize( max_eq );
  
  
  
  
  euclidian_distance_transform( _S, _edt, _matrix1 );



  // ---------------------------------------------------------------------------
  MTvs vec(6), cmt(1);
  
  vec[0] = "Step";
  vec[1] = "Diameter (vx)";
  vec[2] = "Number of voxels occupied by inlet fluid.";
  vec[3] = "Number of voxels occupied by outlet fluid.";
  vec[4] = "Number of voxels occupied by inlet fluid / number of voxels.";
  vec[5] = "Number of voxels occupied by outlet fluid / number of voxels.";
  
  cmt[0] = "YLMs v3";

  saux = _outImgDir + "/" + _outImgRoot + ".dat";
  _dat.open( saux.c_str() );
  abriu( _dat, saux );
  outputFileHead( argc, argv, _dat, vec, cmt );
}















//------------------------------------------------------------------------------
// DESCRIPTION:
//   Run YLM simulation
void Young_Laplace_Method::run( void ){
  string outfile;


  if( print_set_images){
    mymkdir( _outImgDir + "A" );
    
    if( !_wet ){
      mymkdir( _outImgDir + "N0" ); 
      mymkdir( _outImgDir + "N1" ); 
      mymkdir( _outImgDir + "N2" ); 
      mymkdir( _outImgDir + "N3" ); 
    }else{
      mymkdir( _outImgDir + "W0" );    
      mymkdir( _outImgDir + "W1" );    
      mymkdir( _outImgDir + "W2" );    
      mymkdir( _outImgDir + "W3" );
    }

    mymkdir( _outImgDir + "B" ); 
    
    if( !_compressible ){
      mymkdir( _outImgDir + "C0" );
      mymkdir( _outImgDir + "C1" );
      mymkdir( _outImgDir + "C2" );
      mymkdir( _outImgDir + "C3" );
      mymkdir( _outImgDir + "C"  );
    }
    
    mymkdir( _outImgDir + "omega" );
  }
   

  
  for( int step=1; step<=_d.size(); step++ ){
    MTci D = _d.at(step-1);
    MTcd D24 = D*D/4.0;
    double D24cos2 = D24*_cost2;

    cout << "Step " << step << ", D = " << D << " px." << endl;  

    MTcs stepvtk = boost::str( boost::format( "step%02d.vtk" ) % step );
  
  
    //--------------------------------------------------------------------------
    // Set A
    // Result in _matrix1
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _edt[x][y][z]>=D24cos2 ){
        _matrix1[x][y][z] = _F;
      }else{
         _matrix1[x][y][z] = _B;
      }
    }}}

    if(print_set_images) _writevtk( _outImgDir + "A/" + stepvtk, _matrix1 );
    //--------------------------------------------------------------------------




    //--------------------------------------------------------------------------
    // Stage 2: Sets N or W, whether it's a not wetting or wetting fluid   
    // Results in _matrix1
    
    // Not wetting.
    if( ! _wet ){    


      //-------
      // Set N0
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        _matrix2[x][y][z] = _matrix1[x][y][z];
        if( _mm[x][y][z] == _I ){
          _matrix1[x][y][z] = _F;
        }
      }}}
      if(print_set_images) _writevtk( _outImgDir + "N0/" + stepvtk, _matrix1 );


      
      
      // If D is too great, the Inlet Voxel isn't in the eroded part of the image
      // Therefore, all the image is background and there is no necessity to 
      // do component labeling and dilation
      if( _matrix1[_source[0]][_source[1]][_source[2]] == _F ){
        
        
        //-------
        // Set N1
        component_labeling( _matrix1, _F, _B, _next, _tail, _rtable );
        MTci inlet_chamber_label =  _matrix1[_source[0]][_source[1]][_source[2]];
    
        for( int x=0; x<_nx; x++ ){
        for( int y=0; y<_ny; y++ ){
        for( int z=0; z<_nz; z++ ){
          _matrix1[x][y][z] = (_matrix1[x][y][z] == inlet_chamber_label)? _F:_B;
        }}}
    
        if(print_set_images) _writevtk( _outImgDir + "N1/" + stepvtk, _matrix1 );
    

        //-------
        // Set N2
        for( int x=0; x<_nx; x++ ){
        for( int y=0; y<_ny; y++ ){
        for( int z=0; z<_nz; z++ ){
          if( ! (_matrix1[x][y][z] == _F && _matrix2[x][y][z] == _F ) ){
            _matrix1[x][y][z] = _B;
          }
        }}} 
        if(print_set_images) _writevtk( _outImgDir + "N2/" + stepvtk, _matrix1 );


    
        //-------
        // Set N3
        euclidian_distance_transform( _F, _matrix1, _matrix2  );
    
        for( int x=0; x<_nx; x++ ){
        for( int y=0; y<_ny; y++ ){
        for( int z=0; z<_nz; z++ ){

          if( _mm[x][y][z]==_S ){
            _matrix1[x][y][z] = _B;
          }else{
            _matrix1[x][y][z] = (_matrix1[x][y][z] < D24)? _F:_B;
          } 

        }}}
      
        if(print_set_images) _writevtk( _outImgDir + "N3/" + stepvtk, _matrix1 );
        
        
        
      }else{
        for( int x=0; x<_nx; x++ ){
        for( int y=0; y<_ny; y++ ){
        for( int z=0; z<_nz; z++ ){
          _matrix1[x][y][z] = _B;
        }}}
        if(print_set_images) _writevtk( _outImgDir + "N0/" + stepvtk, _matrix1 );
        if(print_set_images) _writevtk( _outImgDir + "N1/" + stepvtk, _matrix1 );
        if(print_set_images) _writevtk( _outImgDir + "N2/" + stepvtk, _matrix1 );
        if(print_set_images) _writevtk( _outImgDir + "N3/" + stepvtk, _matrix1 );
      }




    
    //--------------------------------------------------------------------------
    // Wetting. 
    }else{
      
      //-------
      // Set W0
      euclidian_distance_transform( _F, _matrix1, _matrix2  );
  
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _mm[x][y][z]==_S ){
          _matrix1[x][y][z] = _B;
        }else{
          _matrix1[x][y][z] = (_matrix1[x][y][z] < D24)? _F:_B;
        } 

      }}}
    
      if(print_set_images) _writevtk( _outImgDir + "W0/" + stepvtk, _matrix1 );


      //-------
      // Set W1
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _mm[x][y][z]!=_S ){
          _matrix1[x][y][z] = (_matrix1[x][y][z] == _F)? _B:_F;
        }
      }}}
    
      if(print_set_images) _writevtk( _outImgDir + "W1/" + stepvtk, _matrix1 );


      //-------
      // Set W2
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _mm[x][y][z]==_I ){
          if( _matrix1[x][y][z]==_B )
            _matrix1[x][y][z] = _F;
        }
      }}}
    
      if(print_set_images) _writevtk( _outImgDir + "W2/" + stepvtk, _matrix1 );


      //-------
      // Set W3
      component_labeling( _matrix1, _F, _B, _next, _tail, _rtable );
      MTci inlet_chamber_label =  _matrix1[_source[0]][_source[1]][_source[2]];
  
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        _matrix1[x][y][z] = (_matrix1[x][y][z] == inlet_chamber_label)? _F:_B;
      }}}
    
      if(print_set_images) _writevtk( _outImgDir + "W3/" + stepvtk, _matrix1 );

    }
    //--------------------------------------------------------------------------




    //-------
    // Set B
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _mm[x][y][z] == _I ){
        if( _matrix1[x][y][z] ==_B ){
          _matrix1[x][y][z] = _F;
        }
      }
    }}}
  
    if(print_set_images) _writevtk( _outImgDir + "B/" + stepvtk, _matrix1 );






    //--------------------------------------------------------------------------
    // Stage 3: If not compressible
    // Results in matrix2 for C0 to C1 and in matrix1 thereafter
    // (it's necessary only one loop to omega, because everthing is in matrix1)
    if( !_compressible ){


      //-------
      // Set C0
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        _matrix2[x][y][z] = (_mm[x][y][z]==_O)? _F:_B;
      }}}
      
      if(print_set_images)
        _writevtk( _outImgDir + "C0/" + stepvtk, _matrix2 );



      //-------
      // Set C1
      component_labeling( _matrix2, _F, _B, _next, _tail, _rtable );
      MTci outlet_chamber_label =  _matrix2[_sink[0]][_sink[1]][_sink[2]];
  
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _matrix2[x][y][z]!=_B ){
          if( _matrix2[x][y][z] != outlet_chamber_label  )
            _matrix2[x][y][z]=_F;
          else
            _matrix2[x][y][z]=_B;
        }
      }}}
      
      if(print_set_images)
        _writevtk( _outImgDir + "C1/" + stepvtk, _matrix2 );



      //-------
      // Set C2
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _matrix1[x][y][z]==_F ){
          if( _matrix2[x][y][z]==_F ){
            _matrix1[x][y][z] = _B;
          }else{
            _matrix1[x][y][z] = _F;
          }
        }
      }}}
      
      if(print_set_images)
        _writevtk( _outImgDir + "C2/" + stepvtk, _matrix1 );



      //-------
      // Set C3
      component_labeling( _matrix1, _F, _B, _next, _tail, _rtable );
      MTci inlet_chamber_label =  _matrix1[_source[0]][_source[1]][_source[2]];
  
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _matrix1[x][y][z]!=_B ){
          if( _matrix1[x][y][z] == inlet_chamber_label  )
            _matrix1[x][y][z]=_F;
          else
            _matrix1[x][y][z]=_B;
        }
      }}}
      
      if(print_set_images)
        _writevtk( _outImgDir + "C3/" + stepvtk, _matrix1 );



      //-------
      // Set C
      for( int x=0; x<_nx; x++ ){
      for( int y=0; y<_ny; y++ ){
      for( int z=0; z<_nz; z++ ){
        if( _matrix1[x][y][z] !=_F ){
          if( _mm[x][y][z] == _I ){
            _matrix1[x][y][z] = _F;
          }
        }
      }}}
    
      if(print_set_images)
        _writevtk( _outImgDir + "C/" + stepvtk, _matrix1 );
      
    }
    

    


    //--------------------------------------------------------------------------
    // Omega in _mm
    int Ninlet=0, Noutlet=0, iaux;
    for( int x=0; x<_nx; x++ ){
    for( int y=0; y<_ny; y++ ){
    for( int z=0; z<_nz; z++ ){
      if( _matrix1[x][y][z] ==_F ){
        _mm[x][y][z] = _I;
      }else{
        if( _mm[x][y][z] != _S ){
          _mm[x][y][z] = _O;
        }
      }
      
      iaux = _mm[x][y][z];
      if     ( iaux==_I ) Ninlet++;
      else if( iaux==_O ) Noutlet++;
      
    }}}
    _dat << setprecision(6) << step << " " << D << " " << Ninlet << " " << Noutlet << " " << Noutlet/_n << " " << Ninlet/_n << endl;




    //--------------------------------------------------------------------------
    // Output files?
    if(print_set_images) _writevtk( _outImgDir + "omega/" + stepvtk, _mm );
    
    if( _whichImg=="all"  || (_whichImg=="last" && step==_d.size()-1) )
        _writevtk( _outImgDir + stepvtk, _mm );
  }
}













//------------------------------------------------------------------------------
// DESCRIPTION:
//   Write an VTK data file with matrix data
// 
// PARAMETERS:
//   name  => File name
//   data  => Data matrix
template<typename T>
void Young_Laplace_Method::_writevtk( MTcs &name, const multi_array<T , 3> &data ){

  ofstream FVTK( name.c_str() );
  abriu( FVTK, name );
  
  FVTK << "# vtk DataFile Version 2.0" << endl;
  FVTK << "Geometria" << endl;
  FVTK << "ASCII" << endl;
  FVTK << "DATASET STRUCTURED_POINTS" << endl;
  FVTK << "DIMENSIONS " << _nx << " " << _ny << " " << _nz << endl;
  FVTK << "ASPECT_RATIO 1 1 1" << endl;
  FVTK << "ORIGIN 0 0 0" << endl;
  FVTK << "POINT_DATA " << _nx * _ny * _nz << endl;
  FVTK << "SCALARS geometria int" << endl;
  FVTK << "LOOKUP_TABLE default" << endl;
    
  for( int z=0; z<_nz; z++ ){
  for( int y=0; y<_ny; y++ ){
  for( int x=0; x<_nx; x++ ){
    FVTK << data[x][y][z] << " ";
  }FVTK << endl;
  }}

  FVTK.close();
}







#endif // YOUNG_LAPLACE_METHOD_HPP

