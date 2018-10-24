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


#ifndef MEMF_H
#define MEMF_H

void tricm () ;
int linecount(FILE *fp) ;
void cp(double A[],double B[], double C[]) ;
void centerofmass() ;
void centerofmass2() ;
void Normal(double X1[],double Y1[],double Z1[],double X[]) ;
void Normaldt(int p1,int p2,int p3,double X[]) ;
double Voltri(int p1, int p2, int p3) ;
double Areatri(int p1,int p2,int p3) ;
void TotalVA() ;
void TotalVA0() ;
void force_calculate_dt() ;
void pfdt() ;
double bendenergy() ;
void makeiedge() ;
void filewrite(char filename[]) ;

#endif
