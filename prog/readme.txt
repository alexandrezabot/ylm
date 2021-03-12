# Name of code: Young-Laplace Method
# developer: Alexandre Zabot
#    contact address: 
#          Universidade Federal
#          de Santa Catarina
#          Campus Joinville
#          R. Dona Francisca, 8300 – Bloco U
#          Zona Industrial Norte
#          Joinville – SC – Brasil
#          89.219-600
#    telephone number:
#          +55 47 9 9606 5819
#    e-mail
#          alexandre.zabot@ufsc.br
# year first available: 2021
# hardware required: personal computer
# software required: Python 3+ and C++ compiler
# program language: C++ and Python
# program size: 50 kb
# how to access the source code: 


# Compile:
  g++ -lm -Wall -W -Ofast -std=c++0x -Wno-sign-compare   source/ylm.cpp -o ylm 




# Run examples:


  cd examples/benchmark


  # For 0° contact angle
  
  cd sim/ang00/nc
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat
  
  cd ../ni
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../wc/
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../wi/
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../../../


  # For 45° contact angle

  cd sim/ang45/nc
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat
  
  cd ../ni
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../wc/
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../wi/
  ../../../../../ylm config.txt
  ../../../../vtk2img.py outimg/*vtk
  ../plot_sat.py outimg/d.dat

  cd ../../../








  cd porosimetry
  
  cd sim/ang00
  ../../../../ylm config.txt
  ../../../vtk2img.py outimg/*vtk
  
  cd ../ang20/
  ../../../../ylm config.txt
  ../../../vtk2img.py outimg/*vtk
  
  cd ../ang40/
  ../../../../ylm config.txt
  ../../../vtk2img.py outimg/*vtk
  
  cd ../
  ./plot_sat.py
  
  
  




















