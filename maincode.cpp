/*********************************************************************
**********************************************************************
****      <<< iRBC deformation adhesion Code >>>                  ****
****                                                              ****
****        Anil Kumar Dasanna                                    ****
****        Prof. Dr. Ulrich Schwarz Group                        ****
****        Institute of Theoretical Physics                      ****
****        Heidelberg-69115, Germany                             ****
****        contact:anil.dasanna@bioquant.uni-heidelberg.de       ****
**********************************************************************
*********************************************************************/
#define MAIN
#include "define.h"


void bond_dynamics()
{
	double pon ;
	double B[3],bb ;
	double koff ;	

	for ( int ii = 0 ; ii < vts.size() ; ++ii )
	{	
		if (vts[ii].bondi == 1 )   //already bonded
		{
			double randoff = RanGen[0].Random() ;
			int jj = vts[ii].bondl ;
			int p1modx = vts[ii].ix ;
			int p1mody = vts[ii].iy ; 
			if (vts[ii].poff > randoff)
			{
  	 			vts[ii].bondi = 0 ;
  	 			vts[ii].bondl = randno ;
  	 			ligi[jj] = 0 ; //ligand is freed
			}
		}	
		//for vacant ligands
		else if ( ( vts[ii].z < dthr ) && (vts[ii].bondi == 0) && (vts[ii].knob == 1)  )  
		{
			for ( unsigned int jkl = 0 ; jkl < vts[ii].Pairs.size() ; ++jkl )
			{
					int j = vts[ii].Pairs[jkl].second ;
					double pon = vts[ii].pon ;
					double randon = RanGen[0].Random() ;
					if (( pon > randon ) && (ligi[j] == 0))
					{
						ligi[j] = 1 ; //ligand is filled
						vts[ii].bondi = 1 ;
						vts[ii].bondl = j ;
					}
			 }
		}
	}
}		

int main(int argc, char *argv[])
{
   int seed = (int)time(0) ;
   
   int nthreads = omp_get_max_threads();
   for ( int i = 0 ; i < nthreads ; ++i )
   {
	   int seedn = seed + i*20 ;
	   RanGen.push_back(CRandomMersenne(seedn)) ;
	   sto.push_back(StochasticLib1(seedn)) ;
   }
   
   
   char filename[256];
   int mm ;
   int fx,fy,fz,pi,nn ;
   double A[3],B[3],aa ;
   onrate = atof(argv[1]) ;
   offrate = atof(argv[2]) ;
       
  //ligands positions
  double xdelta = 0.5 ; 
  double ydelta = 0.5 ;
  double xmin = 0.0 ;
  double ymin = 0.0 ;
  
  while ( xmin < Lx )
  {
  	while ( ymin < Ly )
  	{
  		xlig.push_back(xmin) ;
  		ylig.push_back(ymin) ;
  		ligi.push_back(0) ;
  		ymin += ydelta ;
  	}
  	xmin += xdelta ;
  	ymin = 0 ;
  }
   //Solvent initialization
   /*********************************************/
	double vmax = sqrt(3.0*Theta/mass) ;
	int nrt = omp_get_thread_num() ;
 	for ( int j = 0 ; j < Ns ; ++j )
 	{
 	 	solvent.push_back(particle()) ;
		solvent[j].xpos = (Lx*RanGen[nrt].Random()) ;
 	 	solvent[j].ypos = (Ly*RanGen[nrt].Random()) ;
		solvent[j].zpos = (Lz*RanGen[nrt].Random()) ;
		solvent[j].vx = -vmax + 2.0*vmax*RanGen[nrt].Random() ;
		solvent[j].vy = -vmax + 2.0*vmax*RanGen[nrt].Random() ;
		solvent[j].vz = -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	}
	
	/*********************************************/
	
 	 ///////////////////Particles are placed////////////////
 	int count = -1 ;
	int lxi = (int)(Lx/a) ;  
	int lyi = (int)(Ly/a) ; 
	int lzi = (int)(Lz/a) +1 ;
	for ( int i = 0 ; i < lxi ; ++i )
	{
		for ( int j = 0 ; j < lyi ; ++j )
		{
			 for (int k = 0 ; k < lzi ; ++k)
			 {
				   Cells.push_back(cell()) ;
				   count += 1 ;
				   Cells[count].xmin =  (double)i*a ;   
				   Cells[count].xmax =  (double)i*a + a ;
				   Cells[count].ymin =  (double)j*a ;
				   Cells[count].ymax =  (double)j*a + a ;
				   Cells[count].zmin =  -1.0 + (double)k*a ;
				   Cells[count].zmax =  -1.0 + (double)k*a + a ;
				   fx = (int)floor(Cells[count].xmin ) ;
				   fy = (int)floor(Cells[count].ymin ) ;
				   fz = (int)floor(Cells[count].zmin + 1.0) ;
				   cindex[fx][fy][fz] = count ;
			 }
		 }
	 } 
	/**********************************/
	 for ( int i = 0 ; i < Cells.size() ; ++i )
 	 {
		 if ( Cells[i].zmax <= 2.0)
		 {
			 for ( int j = 0 ; j < xlig.size() ; ++j )
			 {
				 double A0 = fabs(((Cells[i].xmin + Cells[i].xmax)*0.5) - xlig[j]) ;
				 double A1 = fabs(((Cells[i].ymin + Cells[i].ymax)*0.5) - ylig[j]) ;
				 A0 -= Lx*round(A0/Lx) ;
				 A1 -= Ly*round(A1/Ly) ;	
				 if ( ( A0 < 1.5 ) && (A1 < 1.5 ) )
				 {
					 Cells[i].liglist.push_back(j) ;
				 }
			 }
		 }
	 }	 
	
   /************************************/
   for ( int i = 0 ; i <  Nt ; ++i )
   {
	   vts.push_back(vert()) ;
	   fmapdt.push_back(force()) ;
   }
   /************************************/
   
   //Reading of Data from file...
   	FILE *fpv ;  //file pointer for membrane files
	double xx, yy, zz ;
	double xxx,yyy,zzz ;
	double Center[3] ;
	int mapi ;
	count = -1 ;
	vmax = sqrt(3.0*Theta/(pmass)) ;
	vector < double > zcrds ;
	 //opening of vertices file..
	fpv = fopen(argv[3],"r" ) ;
	for ( int i = 0 ; i < Nt ; ++i )
	{
		int mmm = fscanf(fpv,"%lf\t%lf\t%lf\n",&xxx,&yyy,&zzz ) ;
		zcrds.push_back(zzz) ;
		vts[i].x = xxx ;
		vts[i].y = yyy ;
		vts[i].z = zzz ;
		vts[i].vx =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
		vts[i].vy =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
		vts[i].vz =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
		fmapdt[i].fx = 0 ;
		fmapdt[i].fy = 0 ;
		fmapdt[i].fz = 0 ;
		vts[i].ix = 0 ;
		vts[i].iy = 0 ;
		vts[i].bondi = 0 ;
		vts[i].bondl = randno ;
		vts[i].bondlen = 0 ;
		vts[i].knob = 0 ;
	}
	fclose(fpv) ; //closing of vertices file
	sort(zcrds.begin(),zcrds.end()) ;
	Center[2] = -zcrds.front() + 0.3 ; //initial height
	Center[0] = Lx/2.0 + 10.0 ;
	Center[1] = Ly/2.0 ;
	
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].x += Center[0] ;
		vts[i].y += Center[1] ;
		vts[i].z += Center[2] ;
	}
	double numer = 1.73205080757*((double)Nt - 2.0) - 5.0*Pi ;
	double denom = 1.73205080757*((double)Nt - 2.0) - 3.0*Pi ;
	theta0 = 180.0*acos(numer/denom)/Pi ;
	//parasite
	/*************************************/
    para[0] = Center[0]-2.0 ;
   	para[1] = Center[1] ;
   	para[2] = Center[2] ;
	vmax = sqrt(3.0*Theta/(PRmass)) ;
	paraV[0] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paraV[1] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paraV[2] =  -vmax + 2.0*vmax*RanGen[nrt].Random() ;
	paradx = 0 ;
	parady = 0 ;
	for ( int i = 0 ; i <  3 ; ++i )
	{
		paraft[i] = 0.0 ;
		parafdt[i] = 0.0 ;
	}
	
	//We now read links(edges) and triangles..
	double l00 ;
    //opening of edges file ..format l1 l2
	fpv = fopen(argv[4],"r") ;
	int nol = linecount(fpv) ;
	fclose(fpv) ;					  
	fpv = fopen(argv[4],"r") ;
	int l1,l2 ;
	double l0 ;
	double meanedge = 0 ;
	double meanedge2 = 0 ;
	//FILE *fuk = fopen("edges.dat","w") ;
	for ( int i = 0 ; i < nol ; ++i )
	{
		int mmm = fscanf(fpv,"%d\t%d\n",&l1,&l2) ;
		//Add this link to edges vector
		edges.push_back(linkk()) ;
		edges[i].l1 = l1-1 ;
		edges[i].l2 = l2-1 ;
		edges[i].kappap = 0 ;
		double distx =  (vts[edges[i].l1].x - vts[edges[i].l2].x)*(vts[edges[i].l1].x - vts[edges[i].l2].x) ;
		double disty =  (vts[edges[i].l1].y - vts[edges[i].l2].y)*(vts[edges[i].l1].y - vts[edges[i].l2].y) ;
		double distz =  (vts[edges[i].l1].z - vts[edges[i].l2].z)*(vts[edges[i].l1].z - vts[edges[i].l2].z) ;
		edges[i].l0 = sqrt( distx + disty + distz ) ;
		//fprintf(fuk,"%lf\n",edges[i].l0 ) ; 
		meanedge += edges[i].l0 ;
		meanedge2 += edges[i].l0 * edges[i].l0 ;
	}
	//fclose(fuk) ;
	cout << "Mean edge is:\t"<< meanedge/(double)nol << endl ;
	cout << "Standard deviation is:\t"<< sqrt((meanedge2/(double)nol) - ((meanedge/(double)nol)*(meanedge/(double)nol))) << endl ;
	fclose(fpv) ;  //closing of links file..
	fpv = fopen(argv[5],"r") ;
	nol = linecount(fpv) ;
	fclose(fpv) ;
	fpv = fopen(argv[5],"r") ;
	int t1,t2,t3 ;
	for ( int i = 0 ; i < nol ; ++i )
	{
			int mmm  = fscanf(fpv,"%d\t%d\t%d\n",&t1,&t2,&t3) ;
			tries.push_back(triangle()) ;
			tries[i].l1 = t1-1 ;
			tries[i].l2 = t2-1 ;
			tries[i].l3 = t3-1 ;
	}
	fclose(fpv) ;
	/*****************************************/
			
	//Now we make other list "edgest" for bending energy
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		edgest.push_back(linkkk()) ;
	}
	int flg1 = 0 ; 
	int flg2 = 0 ;
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		int l1 = edges[i].l1 ;
		int l2 = edges[i].l2 ;
		flg1 = 0 ; 
		flg2 = 0 ;
		for ( int j = 0 ; j < tries.size() ; ++j )
		{
			flg1 = 0 ;
			flg2 = 0 ;
			flg1 += (l1 == tries[j].l1 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l2 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l3 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l1 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l2 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l3 ) ? 1 : 0 ;
			if (( flg1 > 0 ) && (flg2 > 0 ) )
			{
				if (( tries[j].l1 != l1 ) && (tries[j].l1 != l2 ) )
				{
					edgest[i].l1 = tries[j].l1 ;
				}
				if (( tries[j].l2 != l1 ) && (tries[j].l2 != l2 ) )
				{
					edgest[i].l1 = tries[j].l2 ;
				}
				if (( tries[j].l3 != l1 ) && (tries[j].l3 != l2 ) )
				{
					edgest[i].l1 = tries[j].l3 ;
				}
				break ;
			}
		}
		for ( int j = tries.size()-1 ; j >= 0 ; j-- )
		{
			flg1 = 0 ;
			flg2 = 0 ; 
			flg1 += (l1 == tries[j].l1 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l2 ) ? 1 : 0 ;
			flg1 += (l1 == tries[j].l3 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l1 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l2 ) ? 1 : 0 ;
			flg2 += (l2 == tries[j].l3 ) ? 1 : 0 ;
			if (( flg1 > 0 ) && (flg2 > 0 ) )
			{
				if (( tries[j].l1 != l1 ) && (tries[j].l1 != l2 ) )
				{
					edgest[i].l2 = tries[j].l1 ;
				}
				if (( tries[j].l2 != l1 ) && (tries[j].l2 != l2 ) )
				{
					edgest[i].l2 = tries[j].l2 ;
				}
				if (( tries[j].l3 != l1 ) && (tries[j].l3 != l2 ) )
				{
					edgest[i].l2 = tries[j].l3 ;
				}
				break ;
			}
		}
	}
	/////////////////////////////////////////////
	TotalVA0() ;  //initial S.area and volume
	cout << "Sarea:\t" << Sarea0 << endl ;
	cout << "Volume:\t" <<  Volume0 << endl ;
	elap = 0.0005 ; 
	double shear = 1200.0 ;	
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		ll0 = edges[i].l0 ;
		lmax = 2.5*ll0 ;
		double x0 = ll0/lmax ;
		double onemx0 = 1.0 - x0 ;
		double shear1 = 1.0/(4.0*onemx0*onemx0) ;
		double shear2 = x0/(2.0*onemx0*onemx0*onemx0) ;
		double shear4 = (1.732 /(4.0*elap*lmax*x0) ) ;
		double sh = shear4*shear2 - shear4*shear1 + shear4*0.25 ;
		edges[i].kappap = (shear-sh)*4.0*ll0*ll0*ll0/(3.0*1.732) ;
		
	}

	double areaC_modulus = 2.0*shear + kappasa + kappad ;	
	double young_modulus = (4.0*areaC_modulus*shear) / (areaC_modulus + shear) ;
	double poissonR = (areaC_modulus - shear) / (areaC_modulus + shear) ;
	cout << "Shear modulus::\t" << shear << endl ;
	cout << "Young's modulus::\t" << young_modulus << endl ;
	cout << "Poisson's ratio::\t" << poissonR << endl ;


	for ( int i = 0 ; i < tries.size() ; ++i )
	{
		double are = Areatri(tries[i].l1,tries[i].l2,tries[i].l3) ;
		tarea.push_back(are) ; 
	}
	
	/*********************************************/
	//Indexing of edgest
	for ( int i = 0 ; i <  edges.size() ; ++i )
	{
		iedge.push_back(int()) ;
	}
	centerofmass() ;
	makeiedge() ;
	TotalVA() ;
	pfdt() ;
	force_calculate_dt() ;
	
	tricm() ;
	sortsolvent_initial1() ;
	llistupdate() ;
	
	////////////////////////////
    int knob_size = atoi(argv[6]) ;
	vector < double > parea ;
	for ( int i = 0 ; i < tarea.size() ; ++i )
	{
		parea.push_back(tarea[i]/Sarea0) ;
	}
	double aaa = 0 ;
	int knob_count = 0 ;
	
	while ( knob_count < knob_size )
	{
		double rad = RanGen[0].Random() ;
		aaa = 0 ;
		for ( int j = 0 ; j < parea.size() ; ++j )
		{
			aaa += parea[j] ;
			if ( aaa > rad )
			{
				int radi = RanGen[0].IRandom(1,3) ;
				if (( radi == 1 ) && (vts[tries[j].l1].knob == 0 ) )
				{
					vts[tries[j].l1].knob = 1 ;
					knob_count += 1 ;
					break ;
				}
				else if (( radi == 2 ) && (vts[tries[j].l2].knob == 0 ) )
				{
					vts[tries[j].l2].knob = 1 ;
					knob_count += 1 ;
					break ;
				}
				else if (( radi == 3 ) && (vts[tries[j].l3].knob == 0 ) )
				{
					vts[tries[j].l3].knob = 1 ;
					knob_count += 1 ;
					break ;
				}
			}
		}
	}
	
	int cnt = 0 ;
	FILE *fknob = fopen("knobs.dat","w") ;
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		if ( vts[i].knob == 1 )
		{
			fprintf(fknob,"%lf\t%lf\t%lf\n",vts[i].x,vts[i].y,vts[i].z) ;
			cnt += 1 ;
		}
	}
	fclose(fknob) ;
	cout << "no of vertices selected:\t" << cnt << endl ;
	
	//We start iterations here ...

	int filec = 0 ;
	int colsteps = -1 ;
	int filecc = 0 ;
	int bbi = 0 ;
	int upd = 0 ;
	int csteps = (int)(dtcd/dtmd) ;
	FILE *fdat = fopen(argv[7],"w");
	int stl_flag = atoi(argv[8]) ;
	
	//Now we create shear flow...
	int uplist = (int)(0.5/(walls*dtmd)) ;
	for ( int mpcd_steps = 0 ; mpcd_steps < 4000000 ; ++mpcd_steps )
	{
		if ( mpcd_steps < 250000 )
		{
			wallspeed = 0.1 ;
		}
		else
		{
			wallspeed = walls ;
		}
		
		colsteps += 1 ;
		bbi += 1 ;
		upd += 1 ;
		if (mpcd_steps > 0 )	
		{
			filecc += 1 ;
		}
		///////////////////////////////////////////	
		if (colsteps == csteps )
		{
			CollisionStep() ;
			colsteps = 0 ;
	}//End of solvent loop
	
	/***********************************************************/
	#pragma omp  parallel for			
	for ( unsigned int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].vx += (fmapdt[i].fx*dtmd/(2.0*pmass) ) ;
		vts[i].vy += (fmapdt[i].fy*dtmd/(2.0*pmass) ) ;
		vts[i].vz += (fmapdt[i].fz*dtmd/(2.0*pmass) ) ;
		fmapdt[i].fx = 0.0 ;
		fmapdt[i].fy = 0.0 ;
		fmapdt[i].fz = 0.0 ;
	}	
	
	
	//Membrane part
	//centerofmass() ;
	if ( upd == uplist)
	{
		llistupdate() ;
		upd = 0 ;
	}

	//Solvent
	//Streaming + Bounce back .. For Poisellue flow
	//Change Bounce Back boundary conditions
 #pragma omp parallel
 {
	 #pragma omp for private(A,B)
	for ( int i = 0 ; i < Ns ; ++i )
	{
		double fullts = reflectionFed(i,0.0) ;
		solvent[i].xpos += solvent[i].vx*(dtmd-fullts) ;
		solvent[i].ypos += solvent[i].vy*(dtmd-fullts) ;
		solvent[i].zpos += solvent[i].vz*(dtmd-fullts) ;
		double ul = Lz ;
		double ll = 0 ; 
		if (solvent[i].zpos > ul)
		{
			A[0] = solvent[i].xpos - solvent[i].vx*(dtmd-fullts) ;
			A[1] = solvent[i].ypos - solvent[i].vy*(dtmd-fullts) ;
			A[2] = solvent[i].zpos - solvent[i].vz*(dtmd-fullts) ;
			double s = (ul-A[2])/(solvent[i].zpos-A[2]) ; //parameter of line equation
			B[0] = A[0] + s*(solvent[i].xpos-A[0]) ; // point of crossing
			B[1] = A[1] + s*(solvent[i].ypos-A[1]) ;
			B[2] = A[2] + s*(solvent[i].zpos-A[2]) ;
			double dist = sqrt( ((B[0]-A[0])*(B[0]-A[0])) + ((B[1]-A[1])*(B[1]-A[1])) + ((B[2]-A[2])*(B[2]-A[2])) ) ;
			double dtt = dist/sqrt( (solvent[i].vx*solvent[i].vx) + (solvent[i].vy*solvent[i].vy) + (solvent[i].vz*solvent[i].vz)) ;
			solvent[i].vx = wallspeed ;
			solvent[i].vy = 0 ;
			solvent[i].vz = 0 ;

			solvent[i].xpos = B[0] + solvent[i].vx*((dtmd-fullts)-dtt) ;
			solvent[i].ypos = B[1] + solvent[i].vy*((dtmd-fullts)-dtt) ;
			solvent[i].zpos = B[2] + solvent[i].vz*((dtmd-fullts)-dtt) ;

		}
		else if(solvent[i].zpos < ll)
		{
			A[0] = solvent[i].xpos - solvent[i].vx*(dtmd-fullts) ;
			A[1] = solvent[i].ypos - solvent[i].vy*(dtmd-fullts) ;
			A[2] = solvent[i].zpos - solvent[i].vz*(dtmd-fullts) ;
			double s = (ll-A[2])/(solvent[i].zpos-A[2]) ; //parameter of line equation
			B[0] = A[0] + s*(solvent[i].xpos-A[0]) ; // point of crossing
			B[1] = A[1] + s*(solvent[i].ypos-A[1]) ;
			B[2] = A[2] + s*(solvent[i].zpos-A[2]) ;
			double dist = sqrt( ((B[0]-A[0])*(B[0]-A[0])) + ((B[1]-A[1])*(B[1]-A[1])) + ((B[2]-A[2])*(B[2]-A[2])) ) ;
			double dtt = dist/sqrt( (solvent[i].vx*solvent[i].vx) + (solvent[i].vy*solvent[i].vy) + (solvent[i].vz*solvent[i].vz)) ;
			solvent[i].vx *= -1.0 ;
			solvent[i].vy *= -1.0 ;
			solvent[i].vz *= -1.0 ;
			solvent[i].xpos = B[0] + solvent[i].vx*((dtmd-fullts)-dtt) ;
			solvent[i].ypos = B[1] + solvent[i].vy*((dtmd-fullts)-dtt) ;
			solvent[i].zpos = B[2] + solvent[i].vz*((dtmd-fullts)-dtt) ;
		}
		solvent[i].xpos = bbound(solvent[i].xpos,Lx) ;
		solvent[i].ypos = bbound(solvent[i].ypos,Ly) ;
	}
 }
		
	
	/*********************************************************/
	#pragma omp parallel for
	for ( unsigned int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].x += (vts[i].vx*dtmd)  ;
		vts[i].y += (vts[i].vy*dtmd)  ;
        vts[i].z += (vts[i].vz*dtmd)  ;
		vts[i].x = bboundmx(i,vts[i].x,Lx) ;  
		vts[i].y = bboundmy(i,vts[i].y,Ly) ;
	}
 //Parasite position update
  para[0] += (paraV[0]*dtmd) + (paraft[0]*dtmd*dtmd/(2.0*PRmass)) ;
  para[1] += (paraV[1]*dtmd) + (paraft[1]*dtmd*dtmd/(2.0*PRmass)) ;
  para[2] += (paraV[2]*dtmd) + (paraft[2]*dtmd*dtmd/(2.0*PRmass)) ; 	
    if ( para[0] > Lx )
    {
    	para[0] -= Lx ;
		paradx += 1 ; 	
    }
    else if ( para[0] < 0 )
    {
    	para[0] += Lx ;
    	paradx -= 1 ;
    }	
    if ( para[1] > Ly )
    {
    	para[1] -= Ly ;
		parady += 1 ; 	
    }
    else if ( para[1] < 0 )
    {
    	para[1] += Ly ;
    	parady -= 1 ;
    }

 
 	//Computation of force of dt and updating velocity	
	TotalVA() ;
	pfdt() ;  //for parasite
	force_calculate_dt() ;
   // Gravitational force
   #pragma omp parallel for
   for ( int i = 0 ; i < vts.size() ; ++i )
   {
   	fmapdt[i].fz += -Fg/(double)vts.size() ;
   }
   //Bond forces
   int bonds = 0 ;
   double B1[3] ;
   #pragma omp parallel for default(shared) private(A,B1) reduction(+:bonds)
	   for ( int i = 0 ; i < vts.size() ; ++i )
	   {
		if ( vts[i].bondi == 1 )
		{
				bonds += 1 ;
				int jj = vts[i].bondl ;
				A[0] = vts[i].x - xlig[jj] ;
				A[1] = vts[i].y - ylig[jj];
				A[2] = vts[i].z ;
				A[0] -= Lx*round(A[0]/Lx) ;
				A[1] -= Ly*round(A[1]/Ly) ;		
				double aa = sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
				A[0] /= aa ;
				A[1] /= aa ;
				A[2] /= aa ;
				if ( aa > el0 )
				{
					B1[0] = -ks * (aa - el0) * A[0] ;
					B1[1] = -ks * (aa - el0) * A[1] ;
					B1[2] = -ks * (aa - el0) * A[2] ;
					fmapdt[i].fx += B1[0] ;
					fmapdt[i].fy += B1[1] ;
					fmapdt[i].fz += B1[2] ;
				}
			}
		}
	
  /**************************************************/
  #pragma omp parallel for
	  for ( int i = 0 ; i < vts.size() ; ++i )
	  {  	  	   
		  double c = vts[i].z ;
		  double x2 = 0.0 ;
		  if ( c < 0.5 )
		  {
			x2 = 300.0*exp(-50.0*c) ;
			fmapdt[i].fz += x2 ;
		  }
	   }
   
   #pragma omp parallel for
	  for ( int i = 0 ; i < vts.size() ; ++i )
	  {  	  	   
		  double c = Lz - vts[i].z ;
		  double x2 = 0.0 ;
		  if ( c < 0.5 )
		  {
			x2 = 300.0*exp(-50.0*c) ;
			fmapdt[i].fz -= x2 ;
		  }
	   }		
	/*********************************************************/
	//Preparation to bond dynamics
	#pragma omp parallel for private(A,B)
	for ( unsigned int ii = 0 ; ii < vts.size() ; ++ii )
	{
		if ( vts[ii].bondi == 1)
		{
			int p1modx = vts[ii].ix ;
			int p1mody = vts[ii].iy ; 
			int jj = vts[ii].bondl ;
			A[0] = (vts[ii].x + p1modx*Lx) - (xlig[jj] + p1modx*Lx) ;
			A[1] = (vts[ii].y + p1mody*Ly) - (ylig[jj] + p1mody*Ly) ;
			A[2] = vts[ii].z - 0.0 ;				
			double xx = sqrt ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
			if ( xx > el0)
			{
				
				B[0] = -ks * (xx - el0) * A[0]/xx ;
				B[1] = -ks * (xx - el0) * A[1]/xx ;
				B[2] = -ks * (xx - el0) * A[2]/xx ;
			}
			else
			{
				B[0] = 0.0 ;
				B[1] = 0.0 ;
				B[2] = 0.0 ;
			}
			double bb = sqrt( B[0]*B[0] + B[1]*B[1] + B[2]*B[2] ) ;
			double koff = offrate*exp(bb/FD) ;
			vts[ii].poff = 1.0 - exp(-koff*dtmd) ;
		}
		
		else if ( ( vts[ii].z < dthr ) && (vts[ii].bondi == 0) && (vts[ii].knob == 1)  )  
		{
			vts[ii].Pairs.clear() ;
			int pid = vts[ii].pid ;
			double kon = onrate ;
			vts[ii].pon = 1.0 - exp(-kon*dtmd) ;
			for ( unsigned int jkl = 0 ; jkl < Cells[pid].liglist.size() ; ++jkl )
			{
					int j = Cells[pid].liglist[jkl] ;
					int p1modx = vts[ii].ix ;
					int p1mody = vts[ii].iy ;
					Pair apair ; 
					A[0] = (vts[ii].x + p1modx*Lx) - (xlig[j] + p1modx*Lx) ;
					A[1] = (vts[ii].y + p1mody*Ly) - (ylig[j] + p1mody*Ly) ;
					A[2] = vts[ii].z - 0.0 ;
					double aaa =  ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
					if (aaa <= dthr*dthr )
					{
						apair.first = aaa ;
						apair.second = j ;
						vts[ii].Pairs.push_back(apair) ;
					}
			 }
			 sort(vts[ii].Pairs.begin(),vts[ii].Pairs.end()) ;
		}
	}
	//Bond dynamics
	bond_dynamics() ;
 
	 #pragma omp parallel for			
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		vts[i].vx += (fmapdt[i].fx*dtmd/(2.0*pmass) ) ;
		vts[i].vy += (fmapdt[i].fy*dtmd/(2.0*pmass) ) ;
		vts[i].vz += (fmapdt[i].fz*dtmd/(2.0*pmass) ) ;
	}
 		
    //for parasite
	paraV[0] += ((paraft[0]+parafdt[0])*dtmd/(2.0*PRmass) ) ;
	paraV[1] += ((paraft[1]+parafdt[1])*dtmd/(2.0*PRmass) ) ;
	paraV[2] += ((paraft[2]+parafdt[2])*dtmd/(2.0*PRmass) ) ; 
	paraft[0] = parafdt[0] ;
	paraft[1] = parafdt[1] ;
	paraft[2] = parafdt[2] ;
	parafdt[0] = 0.0 ;
	parafdt[1] = 0.0 ;
	parafdt[2] = 0.0 ;
 
	if ( filecc == 5000 ) 
	{
		filec += 1 ;
		fprintf(fdat,"#%d\n",filec) ;		
		for ( int i = 0 ; i < Nt ; ++i )
		{
			fprintf(fdat,"%lf\t%lf\t%lf\t%d\t%d\t%d\n",vts[i].x,vts[i].y,vts[i].z,vts[i].ix,vts[i].iy,vts[i].bondi) ;
		}
		fflush(fdat) ;
		//writing an stl file
		if ( stl_flag == 1 )
		{
			int nm = sprintf(filename,"rbc%d.stl",filec) ;
			filewrite(filename) ;
		}
		filecc = 0 ;
	}
	
}

fclose(fdat) ;

return 0 ;

}
