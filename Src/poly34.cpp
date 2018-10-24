// poly.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
//

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

#include "../define.h"
#define	TwoPi  6.28318530717958648
const double eps=1e-14;
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
int SolveP3(double *x,double aa,double b,double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c
	double a2 = aa*aa;
    double q  = (a2 - 3*b)/9; 
	double r  = (aa*(2*a2-9*b) + 27*c)/54;
    double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
    if(r2<q3) {
        double t=r/sqrt(q3);
		if( t<-1) t=-1;
		if( t> 1) t= 1;
        t=acos(t);
        aa/=3; q=-2*sqrt(q);
        x[0]=q*cos(t/3)-aa;
        x[1]=q*cos((t+TwoPi)/3)-aa;
        x[2]=q*cos((t-TwoPi)/3)-aa;
        return(3);
    } else {
        A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
		if( r<0 ) A=-A;
		B = A==0? 0 : B=q/A;

		aa/=3;
		x[0] =(A+B)-aa;
        x[1] =-0.5*(A+B)-aa;
        x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
        return(1);
    }
}// SolveP3(double *x,double a,double b,double c) {	
