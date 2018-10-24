/*********************************************************************
**********************************************************************
****      <<< iRBC adhesion Code >>>                       	      ****
****                                                              ****
****        Anil Kumar Dasanna                                    ****
****        Prof. Dr. Ulrich Schwarz Group                        ****
****        Institute of Theoretical Physics                      ****
****        Heidelberg-69115, Germany                             ****
****        contact:anil.dasanna@bioquant.uni-heidelberg.de       ****
**********************************************************************
*********************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include <iostream>
#include <new>
#include <vector>
#include <algorithm>
#include <fstream>
//Random number generation files
#include "Include/solfuncs.h"
#include "Include/memfuncs.h"
#include "Include/randomc.h"
#include "Include/stocc.h"
#include "Include/poly34.h"





#define Pi 3.141592653
#define Theta 1.0
#define randno 12345678



//for solvent
#define a 1.0
#define Lx 60.0
#define Ly 40.0
#define Lz 40.0
#define Ns 960000
#define dtcd 0.02	//collision time step
#define navg 10.0   //density
#define Theta 1.0	//Teamperature
#define walls 3.0
#define mass 1.0	//mass
//alpha = 135
#define cosalpha -0.70710678118
#define sinalpha 0.70710678118

// for membrane
#define Nt 600  //no of vertices
#define MAXLINE 400   //used for # of lines in file 
#define kappasa 60000.0 
#define kappad 1200.0    
#define kappav 6000.0  
#define kappab 100.0  //150 for HbAS
#define pmass (navg*mass) 	// in units of mass
#define dtmd 0.002	// integration timestep
#define ks 6000.0
#define dthr 0.50
#define FD 1200.0
#define Fg 40.0  
#define el0 0.1  


#define PRmass 10.0
#define paraR 1.7




using namespace std ;

typedef pair<double,int> Pair ;

 struct particle
 {
 	
	double xpos ;
	double ypos ;
	double zpos ;
 	double vx ;
	double vy ;
	double vz ;
 	int pid ;
 	int indio ; //index that tells the particle is inside or outside
 	bool operator < (const particle& rhs) const { return indio < rhs.indio; }
 };
 
 
 struct cell
 {
	double xmax ;
 	double xmin ;
 	double ymax ;
 	double ymin ;
 	double zmax ;
 	double zmin ;
	double vcm[3] ;
	int nc ;
	double ncc ;
	double rand1,rand2 ;
	double dvv ;  //to compute scaling factor
	double ThE ;
	int flgThE ;
	// vertex list for bounce-back collisions
 	vector <int> llist ;
	vector <int> liglist ;
	vector <int> sol ;
	vector <int> mem ;
 };


struct vert
{
	double x ;
	double y ;
	double z ;
	double vx ;
	double vy ;
	double vz ;
	int pid ;
	int ix ;
	int iy ;
	int bondi ; //bond index
	int bondl ; // index of ligand 
	double bondlen ; //initial bond length
	int knob ;
	double poff ;
	double pon ;
	vector < int > liglist ;
	vector < Pair > Pairs ;
};

struct verto
{
	double x ;
	double y ;
	double z ;
};

	
struct linkk
{
	int l1 ;
	int l2 ;
	double l0 ;
	double kappap ;
};

struct linkkk
{
	int l1 ;
	int l2 ;
};


struct triangle
{
	int l1 ;
	int l2 ;
	int l3 ;
	double cmx ;
	double cmy ;
	double cmz ; 
};

struct force
{
	double fx ;
	double fy ;
	double fz ;
};


#ifndef MAIN
extern vector < CRandomMersenne > RanGen ;
#else
vector < CRandomMersenne > RanGen ;
#endif


#ifndef MAIN
extern vector < StochasticLib1 > sto ;
#else
vector < StochasticLib1 > sto ;
#endif


 //1351.351351351 : conversion factor
#ifndef MAIN
extern double onrate ;
#else
double onrate ;
#endif 


#ifndef MAIN
extern double offrate ;
#else
double offrate ;
#endif


#ifndef MAIN
extern vector < vert > vts ;
#else
vector < vert > vts ;
#endif

#ifndef MAIN
extern vector < vert > vtso ;
#else
vector < vert > vtso ;
#endif


#ifndef MAIN
extern vector < linkk > edges ;
#else
vector < linkk > edges ;
#endif

#ifndef MAIN
extern vector < linkkk > edgest ;
#else
vector < linkkk > edgest ;
#endif

#ifndef MAIN
extern vector < int > iedge ;
#else
vector < int > iedge ;
#endif

#ifndef MAIN
extern vector < triangle > tries ;
#else
vector < triangle > tries ;
#endif


#ifndef MAIN
extern vector < double > tarea ;
#else
vector < double > tarea ;
#endif



#ifndef MAIN
extern vector < force > fmapdt ;
#else
vector < force > fmapdt ;
#endif



#ifndef MAIN
extern double Sarea ;
#else
double Sarea ;
#endif



#ifndef MAIN
extern double Volume ;
#else
double Volume ;
#endif


#ifndef MAIN
extern double center[3] ;
#else
double center[3] ;
#endif


#ifndef MAIN
extern double lmax ;
#else
double lmax ;
#endif


#ifndef MAIN
extern double ll0 ;
#else
double ll0 ;
#endif


#ifndef MAIN
extern double Cq ;
#else
double Cq ;
#endif


#ifndef MAIN
extern double facT,factorT ;
#else
double facT,factorT ;
#endif


#ifndef MAIN
extern vector <particle> solvent ; //particles
#else
vector <particle> solvent ; //particles
#endif



#ifndef MAIN
extern vector <cell> Cells ; //Cells
#else
vector <cell> Cells ; //Cells
#endif


#ifndef MAIN
extern int cindex[(int)(Lx+2)][(int)(Ly+2)][(int)(Lz+2)];
#else
int cindex[(int)(Lx+2)][(int)(Ly+2)][(int)(Lz+2)];
#endif


#ifndef MAIN
extern double Sarea0 ;
#else
double Sarea0 ;
#endif


#ifndef MAIN
extern double Volume0 ;
#else
double Volume0 ;
#endif


#ifndef MAIN
extern int nsolb ;
#else
int nsolb ;
#endif



#ifndef MAIN
extern double elap,theta0 ;
#else
double elap,theta0 ;
#endif


#ifndef MAIN
extern vector < double > xlig ;
#else
vector < double > xlig ;
#endif


#ifndef MAIN
extern vector < double > ylig ;
#else
vector < double > ylig ;
#endif

#ifndef MAIN
extern vector < int > ligi ;
#else
vector < int > ligi ;
#endif


#ifndef MAIN
extern double para[3] ;
#else
double para[3] ;
#endif

#ifndef MAIN
extern double paraV[3] ;
#else
double paraV[3] ;
#endif


#ifndef MAIN
extern double paraft[3] ;
#else
double paraft[3] ;
#endif


#ifndef MAIN
extern double parafdt[3],wallspeed ;
#else
double parafdt[3],wallspeed ;
#endif


#ifndef MAIN
extern int paraid ;
#else
int paraid ;
#endif

#ifndef MAIN
extern int paradx ;
#else
int paradx ;
#endif

#ifndef MAIN
extern int parady ;
#else
int parady ;
#endif

