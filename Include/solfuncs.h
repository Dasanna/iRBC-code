/*********************************************************************
**********************************************************************
****      <<< iRBC adhesion Code >>>                       		  ****
****                                                              ****
****        Anil Kumar Dasanna                                    ****
****        Prof. Dr. Ulrich Schwarz Group                        ****
****        Institute of Theoretical Physics                      ****
****        Heidelberg-69115, Germany                             ****
****        contact:anil.dasanna@bioquant.uni-heidelberg.de       ****
**********************************************************************
*********************************************************************/


#ifndef SOLF_H
#define SOLF_H
#include "stocc.h"
#include "randomc.h"

void sortsolvent_initial1() ;
void CollisionStep() ;
void llistupdate() ;
double Solve(double aa, double bb, double cc ) ;
double reflectionFed(int i,double rstep) ;
double bbound(double x,double lx) ;
double bboundmx(int i,double x,double lx) ;
double bboundmy(int i,double y,double ly) ;
double gsl_ran_gamma (double al, double b) ;
void CollisionDyn() ;
void CollisionDynmem() ;
void CollisionDynpara() ;

#endif
