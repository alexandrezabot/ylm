/**
 * 
 * Biblioteca usada para calcular a Transformada de Distância Euclidiana
 * 
**/


/** operators : Basic arithmetic operation using INFTY numbers
 * 
 * David Coeurjolly (david.coeurjolly@liris.cnrs.fr) - Sept. 2004
 *
**/

#ifndef OPERATORS_HPP
#define OPERATORS_HPP



#include "meusTipos.hpp"






// The original source used a small number as infinity
// Since int range is [-2147483648, 2147483648], I increased it to 1e9
#define INFTY 1000000000
            












/////////Basic functions to handle operations with INFTY

/** 
 **************************************************
 * @b sum
 * @param a Long number with INFTY
 * @param b Long number with INFTY
 * @return The sum of a and b handling INFTY
 **************************************************/
int sum(MTci &a, MTci &b) {
  if ((a==INFTY) || (b==INFTY))     
    return INFTY;    
  else 
    return a+b;
}

/** 
 **************************************************
 * @b prod
 * @param a Long number with INFTY
 * @param b Long number with INFTY
 * @return The product of a and b handling INFTY
 **************************************************/
int prod(MTci &a, MTci &b) {
  if ((a==INFTY) || (b==INFTY)) 
    return INFTY;  
  else 
    return a*b;
}
/** 
 **************************************************
 * @b opp
 * @param a Long number with INFTY
 * @return The opposite of a  handling INFTY
 **************************************************/
int opp (MTci &a) {
  if (a == INFTY) {
    return INFTY;
  }
  else {
    return -a;
  }
}

/** 
 **************************************************
 * @b intdivint
 * @param divid Long number with INFTY
 * @param divis Long number with INFTY
 * @return The division (integer) of divid out of divis handling INFTY
 **************************************************/
int intdivint (MTci &divid, MTci &divis) {
  if (divis == 0) 
    return  INFTY;
  if (divid == INFTY) 
    return  INFTY;
  else 
    return  divid / divis;
}
















////////// Functions F and Sep for the SDT labelling
/** 
 **************************************************
 * @b F
 * @param x 
 * @param i 
 * @param gi2 
 * @return Definition of a parabola
 **************************************************/
int F(MTci &x, MTci &i, MTci &gi2){
  return sum((x-i)*(x-i), gi2);
}

/** 
 **************************************************
 * @b Sep
 * @param i 
 * @param u 
 * @param gi2 
 * @param gu2 
 * @return The abscissa of the intersection point between two parabolas
 **************************************************/
int Sep(MTci &i, MTci &u, MTci &gi2, MTci &gu2) {
  return intdivint(sum( sum( (u*u - i*i),gu2), opp(gi2) ), 2*(u-i));
}
//////////



#endif
