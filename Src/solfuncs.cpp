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

#include "../define.h"


//sorting solvent initially inside and outside
void sortsolvent_initial1()
{
	int cflag = 0 ;
	double A[3] , A1[3], A2[3],X[3];
	double Q1[3],Q2[3] ;
	int count = 0 ;
	double Rd[3] , P0[3];
	//ray direction
	Rd[0] = 1.0 ;
	Rd[1] = 0.0 ;
	Rd[2] = 0.0 ;
	double P2P1[3],P3P2[3],P1P3[3] ;
	double P0P1[3],P0P2[3],P0P3[3] ;
	double f1[3],f2[3],f3[3] ;
	for ( unsigned int i = 0 ; i < solvent.size() ; ++i )
	{
		count = 0 ; 
		Q1[0] = solvent[i].xpos + Rd[0]*0.3 ;
		Q1[1] = solvent[i].ypos + Rd[1]*0.3 ;
		Q1[2] = solvent[i].zpos + Rd[2]*0.3 ;
		Q2[0] = solvent[i].xpos + Rd[0]*0.7 ;
		Q2[1] = solvent[i].ypos + Rd[1]*0.7 ;
		Q2[2] = solvent[i].zpos + Rd[2]*0.7 ;		
		for ( unsigned int j = 0 ; j <  tries.size() ; ++j )
		{
			int p1 = tries[j].l1 ;
			int p2 = tries[j].l2 ;
			int p3 = tries[j].l3 ;
			Normaldt(p1,p2,p3,X) ; 
			double dotP1Q1 = (vts[p1].x-Q1[0])*X[0] + (vts[p1].y-Q1[1])*X[1] + (vts[p1].z-Q1[2])*X[2] ;
			double dotQ2Q1 = (Q2[0]-Q1[0])*X[0] + (Q2[1]-Q1[1])*X[1] + (Q2[2]-Q1[2])*X[2] ;
			P0[0] = Q1[0] + ((dotP1Q1/dotQ2Q1)*(Q2[0]-Q1[0]) ) ;
			P0[1] = Q1[1] + ((dotP1Q1/dotQ2Q1)*(Q2[1]-Q1[1]) ) ;
			P0[2] = Q1[2] + ((dotP1Q1/dotQ2Q1)*(Q2[2]-Q1[2]) ) ;
			P2P1[0] = vts[p2].x - vts[p1].x ;
			P2P1[1] = vts[p2].y - vts[p1].y ;
			P2P1[2] = vts[p2].z - vts[p1].z ;
			P3P2[0] = vts[p3].x - vts[p2].x ;
			P3P2[1] = vts[p3].y - vts[p2].y ;
			P3P2[2] = vts[p3].z - vts[p2].z ;
			P1P3[0] = vts[p1].x - vts[p3].x ;
			P1P3[1] = vts[p1].y - vts[p3].y ;
			P1P3[2] = vts[p1].z - vts[p3].z ;
			
			P0P1[0] = P0[0] - vts[p1].x ;
			P0P1[1] = P0[1] - vts[p1].y ;
			P0P1[2] = P0[2] - vts[p1].z ;
			P0P2[0] = P0[0] - vts[p2].x ;
			P0P2[1] = P0[1] - vts[p2].y ;
			P0P2[2] = P0[2] - vts[p2].z ;
			P0P3[0] = P0[0] - vts[p3].x ;
			P0P3[1] = P0[1] - vts[p3].y ;
			P0P3[2] = P0[2] - vts[p3].z ;
			
			cp(P2P1,P0P1,f1) ;
			cp(P3P2,P0P2,f2) ;
			cp(P1P3,P0P3,f3) ;
			double ch1 = f1[0]*X[0] + f1[1]*X[1] + f1[2]*X[2] ;
			double ch2 = f2[0]*X[0] + f2[1]*X[1] + f2[2]*X[2] ;
			double ch3 = f3[0]*X[0] + f3[1]*X[1] + f3[2]*X[2] ;
			if (( ch1 >= 0 ) && (ch2 >= 0) && (ch3 >= 0) && (tries[j].cmx >= solvent[i].xpos) )
			{
				count += 1 ;
			}			
		}
		if ( (count % 2) == 0 )
		{
			solvent[i].indio = 1 ; // outside
		}
		else if ( (count % 2) == 1 )
		{
			solvent[i].indio = 0 ; // inside
		}
		else
		{
			cout << "error:solcreate" << endl ;
		}
	} // Done with creating indices "inside" or "outside"
	//Now we should sort out solvent array -- to have clear inside first and outside second
	sort(solvent.begin(),solvent.end()) ; //sorting full solvent array
	for ( unsigned int i = 0 ; i < solvent.size() ; ++i)
	{
		if ( ( solvent[i].indio == 0 ) && ( solvent[i+1].indio == 1 ) )
		{
			nsolb = i ;
			break ;
		}
	}
	//cout << "solvent index that divides inside and outside:\t" << nsolb << endl ;	
}


void CollisionStep()
{
	int pi,nn ;
	double A[3],mv[3],aa ;
	double rpos[3] ;
	int fx,fy,fz ;
	//Grid Shifting     
	rpos[0] = a*RanGen[0].Random() ;
	rpos[1] = a*RanGen[0].Random() ;
	rpos[2] = a*RanGen[0].Random() ;
	
	#pragma omp parallel
	{
			#pragma omp for		   
		   for ( unsigned int i = 0 ; i < Cells.size() ; ++i )
		   {
			   Cells[i].vcm[0] = 0.0 ; 
			   Cells[i].vcm[1] = 0.0 ;
			   Cells[i].vcm[2] = 0.0 ;
			   Cells[i].nc = 0 ;
			   Cells[i].ncc = 0 ;
			   Cells[i].dvv = 0 ;
			   Cells[i].ThE = 0 ;
			   Cells[i].flgThE = 0 ;
			   Cells[i].sol.clear() ;
			   Cells[i].mem.clear() ;
		   }
	 }
	//We sort particles in shifted grid
	/**********************************/				
		for ( int i = 0 ; i < solvent.size() ; i++ )
		{
				fx = (int)floor(bbound(solvent[i].xpos - rpos[0],Lx)) ;
				fy = (int)floor(bbound(solvent[i].ypos - rpos[1],Ly)) ;
				fz = (int)floor(solvent[i].zpos - rpos[2]+1.0) ;
				solvent[i].pid = cindex[fx][fy][fz];
				int ii = solvent[i].pid ;
				Cells[ii].sol.push_back(i) ;
		}
			
		for ( int  i = 0 ; i < vts.size() ; ++i )
		{
			fx = (int)floor(bbound(vts[i].x - rpos[0],Lx)) ;
			fy = (int)floor(bbound(vts[i].y - rpos[1],Ly)) ;
			fz = (int)floor(vts[i].z - rpos[2]+1.0) ;
			vts[i].pid = cindex[fx][fy][fz];
			int ii = vts[i].pid ;
			Cells[ii].mem.push_back(i) ;
		}
	
	fx = (int)floor(bbound(para[0] - rpos[0],Lx)) ;
	fy = (int)floor(bbound(para[1] - rpos[1],Ly)) ;
	fz = (int)floor(para[2] - rpos[2]+1.0) ;
	paraid = cindex[fx][fy][fz];
			
   #pragma omp parallel for private(pi)
   for ( int ii = 0 ; ii < Cells.size() ; ++ii )
   { 			   
	   for ( int i = 0 ; i < Cells[ii].sol.size() ; ++i )
	   {
			pi = Cells[ii].sol[i] ;
			Cells[ii].vcm[0] += mass*solvent[pi].vx ;
			Cells[ii].vcm[1] += mass*solvent[pi].vy ;
			Cells[ii].vcm[2] += mass*solvent[pi].vz ;
			Cells[ii].ncc += (1*mass) ;
			Cells[ii].nc += 1 ;
	   }
   }
				
	#pragma omp parallel private(nn)
	{
		int nrrt = omp_get_thread_num() ;
		#pragma omp for
		for ( unsigned int i = 0 ; i <  Cells.size() ; ++i )
		{
				nn = Cells[i].nc ;
				if ( (Cells[i].zmax >= Lz))
				{
					if ( nn < navg )
					{
						Cells[i].vcm[0] += mass*sto[nrrt].Normal((navg-nn)*wallspeed,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].vcm[1] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].vcm[2] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].ncc = navg*mass ;
					}
				}
				else if( Cells[i].zmin <= 0 )
				{
					if ( nn < navg )
					{
						Cells[i].vcm[0] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].vcm[1] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].vcm[2] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
						Cells[i].ncc = navg*mass ;
					}
				}	
			}
		}
	   
	   #pragma omp parallel for private(pi)
	   for ( int ii = 0 ; ii < Cells.size() ; ++ii )
	   {	   
		   for ( int i = 0 ; i < Cells[ii].mem.size() ; ++i )
		   {
			   pi = Cells[ii].mem[i] ;
			   Cells[ii].vcm[0] += pmass*vts[pi].vx ;
			   Cells[ii].vcm[1] += pmass*vts[pi].vy ;
			   Cells[ii].vcm[2] += pmass*vts[pi].vz ;
			   Cells[ii].ncc += (1*pmass) ;
			   Cells[ii].nc += 1 ;
		   }
	   }
		   
	   
		pi = paraid ;
		Cells[pi].vcm[0] += PRmass*paraV[0] ;
		Cells[pi].vcm[1] += PRmass*paraV[1] ;
		Cells[pi].vcm[2] += PRmass*paraV[2] ;
		Cells[pi].ncc += (1*PRmass) ;
		Cells[pi].nc += 1 ;
	   
	   #pragma omp parallel
	   {
		   int nrrt = omp_get_thread_num() ;
		   #pragma omp for private(nn)
		   for ( int ii = 0 ; ii < Cells.size() ; ++ii )
		   {
				nn = Cells[ii].nc ;
				if ( ( nn < navg ) && (Cells[ii].mem.size() > 0 ) )
				{	
					Cells[ii].vcm[0] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
					Cells[ii].vcm[1] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
					Cells[ii].vcm[2] += mass*sto[nrrt].Normal(0,sqrt((navg-nn)*Theta/mass)) ;
					Cells[ii].ncc += (navg-nn)*mass ;
				}
			}
		}    					   
	   
	   
		#pragma omp parallel private(nn)
		{
				#pragma omp for		   
			   for ( unsigned int ii = 0 ; ii < Cells.size() ; ++ii )
			   {
					nn = Cells[ii].ncc ;
					if ( nn != 0 )
					{
						Cells[ii].vcm[0] /= (double)nn ;
						Cells[ii].vcm[1] /= (double)nn ;
						Cells[ii].vcm[2] /= (double)nn ;
					}
			   }
		  }
	   #pragma omp parallel
	   {
			int nrrt = omp_get_thread_num() ;
			 #pragma omp for      
			for ( unsigned int i = 0 ; i < Cells.size() ; ++i )
			{
					Cells[i].rand1 =  2.0*Pi*RanGen[nrrt].Random() ;
					Cells[i].rand2 = -1.0 + (2.0*RanGen[nrrt].Random()) ;
			}
		}
	  		   
	   #pragma omp parallel for private(pi,aa,A)
	   for ( int ii = 0 ; ii < Cells.size() ; ++ii )
	   {
		   for ( int i = 0 ; i <  Cells[ii].sol.size() ; ++i )
		   {
				pi = Cells[ii].sol[i] ;
			   A[0] = solvent[pi].vx - Cells[ii].vcm[0] ;
			   A[1] = solvent[pi].vy - Cells[ii].vcm[1] ;
			   A[2] = solvent[pi].vz - Cells[ii].vcm[2] ;
			   aa = ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
			   Cells[ii].dvv += (mass*aa) ;
		   }
	   }
	   
	   	#pragma omp parallel for private(pi,aa,A)
	   for ( int ii = 0 ; ii < Cells.size() ; ++ii )
	   {	   
		   for ( int i = 0 ; i < Cells[ii].mem.size() ; ++i )
		   {
			   pi = Cells[ii].mem[i] ;
			   A[0] = vts[pi].vx - Cells[ii].vcm[0] ;
			   A[1] = vts[pi].vy - Cells[ii].vcm[1] ;
			   A[2] = vts[pi].vz - Cells[ii].vcm[2] ;
			   aa = ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
			   Cells[ii].dvv += (pmass*aa) ;
		   }
	   }
	   
		
		pi = paraid ;
		A[0] = paraV[0] - Cells[pi].vcm[0] ;
		A[1] = paraV[1] - Cells[pi].vcm[1] ;
		A[2] = paraV[2] - Cells[pi].vcm[2] ;
		aa = ( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
		Cells[pi].dvv += (PRmass*aa) ;				   
	   
		//Now Collision Dynamics ........... 	
		#pragma omp parallel for
		for ( unsigned int i = 0 ; i < Cells.size() ; ++i )
		{
			Cells[i].ThE = 1.0 ;
			if (( Cells[i].nc > 0 ) && (Cells[i].dvv != 0 ))
			{
				Cells[i].ThE = sqrt(3.0*((double)Cells[i].nc - 1.0)*Theta/(Cells[i].dvv) ) ;
			}
		 }					
				
	CollisionDyn() ;
	CollisionDynmem() ;
	CollisionDynpara() ;
}
		 

void llistupdate()
{
  #pragma omp parallel
  {
	 #pragma omp for	
	for ( unsigned int i = 0 ; i < Cells.size() ; ++i )
	{ 
		int count = 0 ;
		int nt = Cells[i].llist.size() ;	
		for ( unsigned int j = 0 ; j <  tries.size() ; ++j )
		{
			int p1 = tries[j].l1 ;
			int p2 = tries[j].l2 ;
			int p3 = tries[j].l3 ;
			double cellx = ( Cells[i].xmin + Cells[i].xmax ) * 0.5 ;
			double celly = ( Cells[i].ymin + Cells[i].ymax ) * 0.5 ;
			double cellz = ( Cells[i].zmin + Cells[i].zmax ) * 0.5 ;								
			
			double xxx1 =  fmod((((vts[p1].x + vts[p1].ix*Lx) +  (vts[p2].x  + vts[p2].ix*Lx)  + (vts[p3].x + vts[p3].ix*Lx))/3.0),Lx) - cellx ;
			double yyy1 =  fmod((((vts[p1].y + vts[p1].iy*Ly) +  (vts[p2].y  + vts[p2].iy*Ly)  + (vts[p3].y + vts[p3].iy*Ly))/3.0),Ly) - celly ;
			xxx1 -= Lx*round(xxx1/Lx) ;
			yyy1 -= Ly*round(yyy1/Ly) ;
			double zzz =  fabs( ((vts[p1].z  +   vts[p2].z   +  vts[p3].z )/3.0) - cellz) ;
			double A1aa = xxx1*xxx1 + yyy1*yyy1 + zzz*zzz  ;
			if ( A1aa < 4.0*a*a )
			{
				count += 1 ;
				if ( count > nt )
				{
					Cells[i].llist.push_back(j) ;
				}
				else
				{
					Cells[i].llist[count-1] = j ;
				}
			}
		}
	}
 }
} 			

			

double Solve(double aa, double bb, double cc )
{
	double *dt ;
	dt = new double[3] ;
	double tt = 0 ;
	int n ;
	n = SolveP3(dt,aa,bb,cc) ;
	if (( n == 3 ) || (n == 2 ) )
	{
		sort(dt,dt+3);
		tt = dt[0] ;
		if (( tt < 0 ) || ( tt > dtmd ) )
		{
			tt = dt[1] ;
			if (( tt < 0 ) || ( tt > dtmd ) )
			{
				tt = dt[2] ;
			}
		}	
	}
	else
	{
		tt = dt[0] ;
	}
	delete [] dt ;
	return tt ;
	
}
		
					

double reflectionFed(int i,double rstep)
{
	double X[3],Xdt[3] ;
	double PSdt[3] ;
	double PSt[3] ;
	double d1[3],d2[3] ;
	double a1[3],a2[3] ;
	double b1[3],b2[3],b3[3] ;
	double a10[3],a20[3] ;
	double d11[3],d22[3] ;
	double vd[3],pd[3] ;
	double vtsp1x1,vtsp1y1,vtsp1z1 ;
	double vtsp2x1,vtsp2y1,vtsp2z1 ;
	double vtsp3x1,vtsp3y1,vtsp3z1 ;
	double vtsp1x2,vtsp1y2,vtsp1z2 ;
	double vtsp2x2,vtsp2y2,vtsp2z2 ;
	double vtsp3x2,vtsp3y2,vtsp3z2 ;	
	double fulldt = rstep;
	double A1[3],B1[3],C1[3] ;
	double A2[3],B2[3],C2[3] ;
	int flg = 0 ;
	int p1modx = 0, p2modx = 0, p3modx = 0 ;	
	int p1mody = 0, p2mody = 0, p3mody = 0 ;
	int fx = (int)floor(bbound(solvent[i].xpos,Lx)) ;
	int fy = (int)floor(bbound(solvent[i].ypos,Ly)) ; 
	int fz = (int)floor(solvent[i].zpos+1.0) ;
	int pid = cindex[fx][fy][fz];	
	for ( unsigned int j = 0 ; j < Cells[pid].llist.size() ; ++j )
	{
			flg = 0 ;
			int jj = Cells[pid].llist[j] ;
			int p1 = tries[jj].l1 ;
			int p2 = tries[jj].l2 ;
			int p3 = tries[jj].l3 ;

			/*p1modx = vts[p1].ix ;
			p2modx = vts[p2].ix ;
			p3modx = vts[p3].ix ;*/
			p1modx = 0 ;
			p2modx = 0 ;
			p3modx = 0 ;
			int xx1 = min(p1modx,p2modx) ;
			int xx = min(xx1,p3modx) ;
				
			/*p1mody = vts[p1].iy ;
			p2mody = vts[p2].iy ;
			p3mody = vts[p3].iy ;*/
			p1mody = 0 ;
			p2mody = 0 ;
			p3mody = 0 ;
			xx1 = min(p1mody,p2mody) ;
			int yy = min(xx1,p3mody) ;
			double solvx = solvent[i].xpos + xx*Lx ;
			double solvy = solvent[i].ypos + yy*Ly ;
			double solvz = solvent[i].zpos ;
			
			
			A1[0] = vtsp1x1 = (vts[p1].x + p1modx*Lx) + vts[p1].vx*rstep;
			A1[1] = vtsp1y1 = (vts[p1].y + p1mody*Ly) + vts[p1].vy*rstep;
			A1[2] = vtsp1z1 = vts[p1].z + vts[p1].vz*rstep;
			B1[0] = vtsp2x1 = (vts[p2].x + p2modx*Lx) + vts[p2].vx*rstep;
			B1[1] = vtsp2y1 = (vts[p2].y + p2mody*Ly) + vts[p2].vy*rstep;
			B1[2] = vtsp2z1 = vts[p2].z + vts[p2].vz*rstep;
			C1[0] = vtsp3x1 = (vts[p3].x + p3modx*Lx) + vts[p3].vx*rstep;
			C1[1] = vtsp3y1 = (vts[p3].y + p3mody*Ly) + vts[p3].vy*rstep;
			C1[2] = vtsp3z1 = vts[p3].z + vts[p3].vz*rstep;
			
			/**************************************************/
			double rstencilx = Lx - (3.0*a) ;
			double lstencilx = 3.0*a ;
			double rstencily = Ly - (3.0*a) ;
			double lstencily = 3.0*a ;
						
			if ( solvx < lstencilx )
			{
				A1[0] = vtsp1x1 -= (vtsp1x1 > rstencilx)?Lx:0 ;
				B1[0] = vtsp2x1 -= (vtsp2x1 > rstencilx)?Lx:0 ;
				C1[0] = vtsp3x1 -= (vtsp3x1 > rstencilx)?Lx:0 ;
			}
			else if ( solvx > rstencilx)
			{
				A1[0] = vtsp1x1 += (vtsp1x1 < lstencilx)?Lx:0 ;
				B1[0] = vtsp2x1 += (vtsp2x1 < lstencilx)?Lx:0 ;
				C1[0] = vtsp3x1 += (vtsp3x1 < lstencilx)?Lx:0 ;
			}
			if ( solvy < lstencily )
			{
				A1[1] = vtsp1y1 -= (vtsp1y1 > rstencily)?Ly:0 ;
				B1[1] = vtsp2y1 -= (vtsp2y1 > rstencily)?Ly:0 ;
				C1[1] = vtsp3y1 -= (vtsp3y1 > rstencily)?Ly:0 ;
			}
			else if ( solvy > rstencily)
			{
				A1[1] = vtsp1y1 += (vtsp1y1 < lstencily)?Ly:0 ;
				B1[1] = vtsp2y1 += (vtsp2y1 < lstencily)?Ly:0 ;
				C1[1] = vtsp3y1 += (vtsp3y1 < lstencily)?Ly:0 ;
			}
			/**************************************************/			
			A2[0] = vtsp1x2 = vtsp1x1 + vts[p1].vx*(dtmd-rstep);
			A2[1] = vtsp1y2 = vtsp1y1 + vts[p1].vy*(dtmd-rstep);
			A2[2] = vtsp1z2 = vtsp1z1 + vts[p1].vz*(dtmd-rstep);
			B2[0] = vtsp2x2 = vtsp2x1 + vts[p2].vx*(dtmd-rstep);
			B2[1] = vtsp2y2 = vtsp2y1 + vts[p2].vy*(dtmd-rstep);
			B2[2] = vtsp2z2 = vtsp2z1 + vts[p2].vz*(dtmd-rstep);
			C2[0] = vtsp3x2 = vtsp3x1 + vts[p3].vx*(dtmd-rstep);
			C2[1] = vtsp3y2 = vtsp3y1 + vts[p3].vy*(dtmd-rstep);
			C2[2] = vtsp3z2 = vtsp3z1 + vts[p3].vz*(dtmd-rstep);	
			
			Normal(A1,B1,C1,X) ;
			Normal(A2,B2,C2,Xdt) ;
			PSt[0] = solvx - vtsp1x1 ;
			PSt[1] = solvy - vtsp1y1 ;
			PSt[2] = solvz - vtsp1z1 ;
			PSdt[0] = (solvx  + solvent[i].vx*(dtmd-rstep)) -  vtsp1x2;
			PSdt[1] = (solvy  + solvent[i].vy*(dtmd-rstep)) -  vtsp1y2;
			PSdt[2] = (solvz  + solvent[i].vz*(dtmd-rstep)) -  vtsp1z2 ;
			double bt = X[0]*PSt[0] + X[1]*PSt[1] + X[2]*PSt[2] ;
			double bdt = Xdt[0]*PSdt[0] + Xdt[1]*PSdt[1] + Xdt[2]*PSdt[2] ;
			if (( solvent[i].indio == 1 ) && ( bt < 0 ) )
			{
				flg = 1 ;
			}
			if (( solvent[i].indio == 0 ) && ( bt > 0 ) )
			{
				flg = 1 ;
			}			
			if ( ( bt*bdt <= 0 ) && (flg == 0) )
			{
				/****************************/
				a10[0] = vtsp2x1 - vtsp1x1 ;
				a10[1] = vtsp2y1 - vtsp1y1 ;
				a10[2] = vtsp2z1 - vtsp1z1 ;
				a20[0] = vtsp3x1 - vtsp1x1 ;
				a20[1] = vtsp3y1 - vtsp1y1 ;
				a20[2] = vtsp3z1 - vtsp1z1 ;
				d1[0] = vts[p2].vx - vts[p1].vx ;
				d1[1] = vts[p2].vy - vts[p1].vy ;
				d1[2] = vts[p2].vz - vts[p1].vz ;
				d2[0] = vts[p3].vx - vts[p1].vx ;
				d2[1] = vts[p3].vy - vts[p1].vy ;
				d2[2] = vts[p3].vz - vts[p1].vz ;						
				/****************************/						
				a1[0] = PSt[0] ;
				a1[1] = PSt[1] ;
				a1[2] = PSt[2] ;
				a2[0] = solvent[i].vx - vts[p1].vx ;
				a2[1] = solvent[i].vy - vts[p1].vy ;
				a2[2] = solvent[i].vz - vts[p1].vz ;
				cp(a10,a20,b1) ;
				cp(a10,d2,d11) ;
				cp(d1,a20,d22) ;
				cp(d1,d2,b3) ;
				b2[0] = d11[0] + d22[0] ;
				b2[1] = d11[1] + d22[1] ;
				b2[2] = d11[2] + d22[2] ;
				double a2b3 = a2[0]*b3[0] + a2[1]*b3[1] + a2[2]*b3[2] ;
				double a1b3 = a1[0]*b3[0] + a1[1]*b3[1] + a1[2]*b3[2] ;
				double a2b2 = a2[0]*b2[0] + a2[1]*b2[1] + a2[2]*b2[2] ;
				double a1b2 = a1[0]*b2[0] + a1[1]*b2[1] + a1[2]*b2[2] ;
				double a2b1 = a2[0]*b1[0] + a2[1]*b1[1] + a2[2]*b1[2] ;
				double a1b1 = a1[0]*b1[0] + a1[1]*b1[1] + a1[2]*b1[2] ;
				double tt = Solve((a1b3+a2b2)/a2b3,(a1b2+a2b1)/a2b3,a1b1/a2b3) ;
				PSdt[0] = (solvx  + solvent[i].vx*tt) - (vtsp1x1 + vts[p1].vx*tt) ;
				PSdt[1] = (solvy  + solvent[i].vy*tt) - (vtsp1y1 + vts[p1].vy*tt) ;
				PSdt[2] = (solvz  + solvent[i].vz*tt) - (vtsp1z1 + vts[p1].vz*tt) ;
				pd[0] = (solvx  + solvent[i].vx*tt) ;
				pd[1] = (solvy  + solvent[i].vy*tt) ;
				pd[2] = (solvz  + solvent[i].vz*tt) ;
				a10[0] = (vtsp2x1 + vts[p2].vx*tt)  - (vtsp1x1 + vts[p1].vx*tt) ;
				a10[1] = (vtsp2y1 + vts[p2].vy*tt)  - (vtsp1y1 + vts[p1].vy*tt) ;
				a10[2] = (vtsp2z1 + vts[p2].vz*tt)  - (vtsp1z1 + vts[p1].vz*tt) ;
				a20[0] = (vtsp3x1 + vts[p3].vx*tt)  - (vtsp1x1 + vts[p1].vx*tt) ;
				a20[1] = (vtsp3y1 + vts[p3].vy*tt)  - (vtsp1y1 + vts[p1].vy*tt) ;
				a20[2] = (vtsp3z1 + vts[p3].vz*tt)  - (vtsp1z1 + vts[p1].vz*tt) ;
				double pa1 = PSdt[0]*a10[0] + PSdt[1]*a10[1] + PSdt[2]*a10[2] ;
				double pa2 = PSdt[0]*a20[0] + PSdt[1]*a20[1] + PSdt[2]*a20[2] ;
				double a1a2 = a10[0]*a20[0] + a10[1]*a20[1] + a10[2]*a20[2] ; 
				double a12 = a10[0]*a10[0] + a10[1]*a10[1] + a10[2]*a10[2] ;
				double a22 = a20[0]*a20[0] + a20[1]*a20[1] + a20[2]*a20[2] ;
				double denom = (a12*a22) - (a1a2*a1a2) ;
				double eta = ((pa1*a22)-(pa2*a1a2) ) / denom ;
				double zeta = ((pa2*a12)-(pa1*a1a2) ) / denom ;
				if (( eta >= 0 ) && (zeta >= 0 ) && (eta+zeta <= 1) && (tt < (dtmd-rstep)) )
				{
					fulldt += tt ;
					vd[0] = (1.0-eta-zeta)*vts[p1].vx + eta*vts[p2].vx + zeta*vts[p3].vx ;
					vd[1] = (1.0-eta-zeta)*vts[p1].vy + eta*vts[p2].vy + zeta*vts[p3].vy ;
					vd[2] = (1.0-eta-zeta)*vts[p1].vz + eta*vts[p2].vz + zeta*vts[p3].vz ;
					solvent[i].vx = 2.0*vd[0] - solvent[i].vx ;
					solvent[i].vy = 2.0*vd[1] - solvent[i].vy ;
					solvent[i].vz = 2.0*vd[2] - solvent[i].vz ;
					solvent[i].xpos = pd[0] ;
					solvent[i].ypos = pd[1] ;
					solvent[i].zpos = pd[2] ;
				}
			}
		}
		return fulldt ;
}
				
				

//used for periodic boundaries (wrapping)
double bbound(double x,double lx)
{
    double r = 0.0 ;
    if ( x > lx )
    {
    	while(x > lx)
    	{
    		x -= lx ;
    	}
    }
    if ( x < 0 )
    {
    	while ( x < 0)
    	{
    		x += lx ;
    	}
    }
    r = x ;
    return r ;
} 

//used for periodic boundaries (wrapping) membrane x
double bboundmx(int i,double x,double lx)
{
    double r = 0.0 ;
    if ( x > lx )
    {
    	x -= lx ;
		vts[i].ix += 1 ; 	
    }
    if ( x < 0 )
    {
    	x += lx ;
    	vts[i].ix -= 1 ;
    }
    r = x ;
    return r ;
} 


//used for periodic boundaries (wrapping) membrane y
double bboundmy(int i,double y,double ly)
{
    double r = 0.0 ;
    if ( y > ly )
    {
    	y -= ly ;
		vts[i].iy += 1 ; 	
    }
    if ( y < 0 )
    {
    	y += ly ;
    	vts[i].iy -= 1 ;
    }
    r = y ;
    return r ;
} 

//Random no generation from Gamma distribution
double gsl_ran_gamma (double al, double b) //a :: alpha > 1.0 b :: is  always 1.0 
{
 
  if (al < 1)
    {
      double u = RanGen[0].Random();
      return gsl_ran_gamma (1.0 + al, b) * pow (u, 1.0 / al);
    }
 
  {
    double x, v, u;
    double d = al - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while (1)
      {
        do
          {
            x = sto[0].Normal(0.0, 1.0);
            v = 1.0 + c * x;
          }
        while (v <= 0);
 
        v = v * v * v;
        u = RanGen[0].Random();
 
        if (u < 1 - 0.0331 * x * x * x * x) 
          break;
 
        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
 
    return b * d * v;
  }
}


void CollisionDyn()
{
	#pragma omp parallel
	{
		#pragma omp for
		for ( unsigned int i = 0  ; i < solvent.size() ; ++i )
		{
			int pi = solvent[i].pid ;
			double dvv[3] ;
			double dvpn[3],dvp[3],dvn[3],unit[3] ;
			double A[3] ;			
			double thet = Cells[pi].rand1 ;
			double zz = Cells[pi].rand2 ;
			double unitx = sqrt(1.0 - (zz*zz))*cos(thet) ;
			double unity = sqrt(1.0 - (zz*zz))*sin(thet) ;
			double unitz = zz ;	
			double scalev = Cells[pi].ThE ;
			dvv[0] = solvent[i].vx - Cells[pi].vcm[0] ;
			dvv[1] = solvent[i].vy - Cells[pi].vcm[1] ;
			dvv[2] = solvent[i].vz - Cells[pi].vcm[2] ;
			unit[0] = unitx ;
			unit[1] = unity ;
			unit[2] = unitz ;
			double xxx = dvv[0]*unit[0] + dvv[1]*unit[1] + dvv[2]*unit[2] ;
			dvp[0] = xxx*unit[0] ;
			dvp[1] = xxx*unit[1] ;
			dvp[2] = xxx*unit[2] ;
			dvn[0] = dvv[0] - dvp[0] ;
			dvn[1] = dvv[1] - dvp[1] ;
			dvn[2] = dvv[2] - dvp[2] ;
			cp(dvn,unit,dvpn) ;
			solvent[i].vx = Cells[pi].vcm[0] + scalev*(cosalpha*dvn[0] + sinalpha*dvpn[0] + dvp[0]) ;
			solvent[i].vy = Cells[pi].vcm[1] + scalev*(cosalpha*dvn[1] + sinalpha*dvpn[1] + dvp[1]) ;
			solvent[i].vz = Cells[pi].vcm[2] + scalev*(cosalpha*dvn[2] + sinalpha*dvpn[2] + dvp[2]) ;
		}
	}
}


void CollisionDynmem()
{
	#pragma omp parallel
	{
		#pragma omp for
		for ( unsigned int i = 0  ; i < vts.size() ; ++i )
		{
			double dvv[3] ;
			double dvpn[3],dvp[3],dvn[3],unit[3] ;
			double A[3] ;
			double alpha = Pi / 2.0 ;
			int pi = vts[i].pid ;
			double thet = Cells[pi].rand1 ;
			double zz = Cells[pi].rand2 ;
			double unitx = sqrt(1.0 - (zz*zz))*cos(thet) ;
			double unity = sqrt(1.0 - (zz*zz))*sin(thet) ;
			double unitz = zz ;	
			double scalev = Cells[pi].ThE ;
			dvv[0] = vts[i].vx - Cells[pi].vcm[0] ;
			dvv[1] = vts[i].vy - Cells[pi].vcm[1] ;
			dvv[2] = vts[i].vz - Cells[pi].vcm[2] ;
			unit[0] = unitx ;
			unit[1] = unity ;
			unit[2] = unitz ;
			double xxx = dvv[0]*unit[0] + dvv[1]*unit[1] + dvv[2]*unit[2] ;
			dvp[0] = xxx*unit[0] ;
			dvp[1] = xxx*unit[1] ;
			dvp[2] = xxx*unit[2] ;
			dvn[0] = dvv[0] - dvp[0] ;
			dvn[1] = dvv[1] - dvp[1] ;
			dvn[2] = dvv[2] - dvp[2] ;
			cp(dvn,unit,dvpn) ;
			vts[i].vx = Cells[pi].vcm[0] + scalev*(cosalpha*dvn[0] + sinalpha*dvpn[0] + dvp[0]) ;
			vts[i].vy = Cells[pi].vcm[1] + scalev*(cosalpha*dvn[1] + sinalpha*dvpn[1] + dvp[1]) ;
			vts[i].vz = Cells[pi].vcm[2] + scalev*(cosalpha*dvn[2] + sinalpha*dvpn[2] + dvp[2]) ;
		}
	}
}


void CollisionDynpara()
{
	double dvv[3],dv,dvvv;
	double dvpn[3],dvp[3],dvn[3],unit[3] ;
	double A[3] ;
	
		double alpha = Pi / 2.0 ;
		int pi = paraid ;
		double thet = Cells[pi].rand1 ;
		double zz = Cells[pi].rand2 ;
		double unitx = sqrt(1.0 - (zz*zz))*cos(thet) ;
		double unity = sqrt(1.0 - (zz*zz))*sin(thet) ;
		double unitz = zz ;	
		double scalev = Cells[pi].ThE ;
		dvv[0] = paraV[0] - Cells[pi].vcm[0] ;
		dvv[1] = paraV[1] - Cells[pi].vcm[1] ;
		dvv[2] = paraV[2] - Cells[pi].vcm[2] ;
		unit[0] = unitx ;
		unit[1] = unity ;
		unit[2] = unitz ;
		double xxx = dvv[0]*unit[0] + dvv[1]*unit[1] + dvv[2]*unit[2] ;
		dvp[0] = xxx*unit[0] ;
		dvp[1] = xxx*unit[1] ;
		dvp[2] = xxx*unit[2] ;
		dvn[0] = dvv[0] - dvp[0] ;
		dvn[1] = dvv[1] - dvp[1] ;
		dvn[2] = dvv[2] - dvp[2] ;
		cp(dvn,unit,dvpn) ;
		paraV[0] = Cells[pi].vcm[0] + scalev*(cosalpha*dvn[0] + sinalpha*dvpn[0] + dvp[0]) ;
		paraV[1] = Cells[pi].vcm[1] + scalev*(cosalpha*dvn[1] + sinalpha*dvpn[1] + dvp[1]) ;
		paraV[2] = Cells[pi].vcm[2] + scalev*(cosalpha*dvn[2] + sinalpha*dvpn[2] + dvp[2]) ;
}
