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



//computation of center of mass of triangles
void tricm ()
{
	for ( unsigned int i = 0 ; i < tries.size() ; ++i )
	{
		int p1 = tries[i].l1 ;
		int p2 = tries[i].l2 ;
		int p3 = tries[i].l3 ;
		tries[i].cmx = ( (vts[p1].x + vts[p1].ix*Lx) + (vts[p2].x + vts[p2].ix*Lx) + (vts[p3].x + vts[p3].ix*Lx) ) / 3.0 ;
		tries[i].cmy = ( (vts[p1].y + vts[p1].iy*Ly) + (vts[p2].y + vts[p2].iy*Ly) + (vts[p3].y + vts[p3].iy*Ly) ) / 3.0 ;
		tries[i].cmz = ( vts[p1].z  + vts[p2].z  + vts[p3].z ) / 3.0 ; 
		tries[i].cmx = fmod(tries[i].cmx,Lx ) ;
		tries[i].cmy = fmod(tries[i].cmy,Ly ) ;
	}
}




int linecount(FILE *fp)
{
    char buff[MAXLINE];
    int count = 0;
    while(fgets(buff,MAXLINE,fp) != NULL)
    {
        count++;
    }
    return count;
}

void cp(double A[],double B[], double C[])
{
	C[0] = (A[1] * B[2]) - (B[1] * A[2]) ;
	C[1] = (A[2] * B[0]) - (B[2] * A[0]) ;
	C[2] = (A[0] * B[1]) - (B[0] * A[1]) ;
}


void centerofmass()
{
	double rr[3] ;
	rr[0] = 0 ;
	rr[1] = 0 ;
	rr[2] = 0 ;
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		rr[0] += vts[i].x + vts[i].ix*Lx ;
		rr[1] += vts[i].y + vts[i].iy*Ly ;
		rr[2] += vts[i].z ;
	}
	
	center[0] = rr[0]/(double)vts.size() ;
	center[1] = rr[1]/(double)vts.size() ;
	center[2] = rr[2]/(double)vts.size() ;
	center[0] = fmod(center[0],Lx) ;
	center[1] = fmod(center[1],Ly) ;
	center[2] = fmod(center[2],Lz) ;
}

void centerofmass2()
{
	double rr[3] ;
	rr[0] = 0 ;
	rr[1] = 0 ;
	rr[2] = 0 ;
	for ( int i = 0 ; i < vts.size() ; ++i )
	{
		rr[0] += vts[i].x + vts[i].ix*Lx ;
		rr[1] += vts[i].y + vts[i].iy*Ly ;
		rr[2] += vts[i].z ;
	}
	
	center[0] = rr[0]/(double)vts.size() ;
	center[1] = rr[1]/(double)vts.size() ;
	center[2] = rr[2]/(double)vts.size() ;
}	



void Normal(double X1[],double Y1[],double Z1[],double X[])
{
		double A[3],B[3] ;
		X[0] = X[1] = X[2] = 0.0 ;
		A[0] = Y1[0] - X1[0] ;
		A[1] = Y1[1] - X1[1] ;
		A[2] = Y1[2] - X1[2] ;
		B[0] = Z1[0] - X1[0] ;
		B[1] = Z1[1] - X1[1] ;
		B[2] = Z1[2] - X1[2] ;
		cp(A,B,X) ;

}

void Normaldt(int p1,int p2,int p3,double X[])
{
		double A[3],B[3],D[3],C[3] ;
		X[0] = X[1] = X[2] = 0.0 ;
		A[0] = (vts[p2].x + vts[p2].ix*Lx) - (vts[p1].x + vts[p1].ix*Lx) ;
		A[1] = (vts[p2].y + vts[p2].iy*Ly) - (vts[p1].y + vts[p1].iy*Ly) ;
		A[2] = vts[p2].z - vts[p1].z ;
		B[0] = (vts[p3].x + vts[p3].ix*Lx) - (vts[p1].x + vts[p1].ix*Lx) ;
		B[1] = (vts[p3].y + vts[p3].iy*Ly) - (vts[p1].y + vts[p1].iy*Ly) ;
		B[2] = vts[p3].z - vts[p1].z ;
		cp(A,B,X) ;

}

double Voltri(int p1, int p2, int p3) 
{
	double A[3],B[3],C[3],D[3],Cn[3] ;
	A[0] = vts[p1].x + vts[p1].ix*Lx ;
	A[1] = vts[p1].y + vts[p1].iy*Ly ;
	A[2] = vts[p1].z ;
	B[0] = vts[p2].x + vts[p2].ix*Lx  ;
	B[1] = vts[p2].y + vts[p2].iy*Ly ;
	B[2] = vts[p2].z ;
	C[0] = vts[p3].x + vts[p3].ix*Lx ;
	C[1] = vts[p3].y + vts[p3].iy*Ly ;
	C[2] = vts[p3].z ;
	Cn[0] = (A[0] + B[0] + C[0] ) / 3.0 ;
	Cn[1] = (A[1] + B[1] + C[1] ) / 3.0 ;
	Cn[2] = (A[2] + B[2] + C[2] ) / 3.0 ;
	Normaldt(p1,p2,p3,D) ;
	double vol = (1.0/6.0)*( (D[0]*Cn[0]) + (D[1]*Cn[1]) + (D[2]*Cn[2]) ) ;
	return vol ;
}

double Areatri(int p1,int p2,int p3)
{
		double A[3],B[3],eta[3] ;
		double area ;
		A[0] = (vts[p2].x + vts[p2].ix*Lx) - (vts[p1].x + vts[p1].ix*Lx) ;
		A[1] = (vts[p2].y + vts[p2].iy*Ly) - (vts[p1].y + vts[p1].iy*Ly) ;
		A[2] = vts[p2].z - vts[p1].z ;
		B[0] = (vts[p3].x + vts[p3].ix*Lx) - (vts[p1].x + vts[p1].ix*Lx) ;
		B[1] = (vts[p3].y + vts[p3].iy*Ly) - (vts[p1].y + vts[p1].iy*Ly) ;
		B[2] = vts[p3].z - vts[p1].z ;
		cp(A,B,eta) ;
		area = 0.5*sqrt(eta[0]*eta[0] + eta[1]*eta[1] + eta[2]*eta[2] ) ; //area of this triangle
		return area ;
}



void TotalVA()
{
	double totV = 0 ;
	double totA = 0 ;
	double A[3],B[3],eta[3] ;
	double factor = 1.0 ;
	int flag = 0 ;
	
	#pragma omp parallel for reduction(+:totA)
	for ( int i = 0 ; i < tries.size() ; ++i )
	{
		int l1 = tries[i].l1 ;
		int l2 = tries[i].l2 ;
		int l3 = tries[i].l3 ;	
		//Area
		totA += Areatri(l1,l2,l3) ;
	}
	Sarea = totA ;
	
	#pragma omp parallel for reduction(+:totV)
	for ( int i = 0 ; i < tries.size() ; ++i )
	{
		int l1 = tries[i].l1 ;
		int l2 = tries[i].l2 ;
		int l3 = tries[i].l3 ;	
		//Volume
		totV += Voltri(l1,l2,l3) ;
	}
	Volume = totV ;
}

void TotalVA0()
{
	double totV = 0 ;
	double totA = 0 ;
	double A[3],B[3],eta[3] ;
	double factor = 1.0 ;
	int flag = 0 ;
	for ( int i = 0 ; i < tries.size() ; ++i )
	{
		int l1 = tries[i].l1 ;
		int l2 = tries[i].l2 ;
		int l3 = tries[i].l3 ;	
		//Area
		totA += Areatri(l1,l2,l3) ;
		//Volume
		totV += Voltri(l1,l2,l3) ;
	}
	Sarea0 = totA ;
	Volume0 = totV ;
}


void force_calculate_dt()
{
	double *ffx,*ffy,*ffz ;
  #pragma omp parallel private(ffx,ffy,ffz) default(shared)
  {
	double A[3],A21[3],A13[3],A32[3],B[3],eta[3],C[3] ;
	double a21[3],a31[3],a24[3],a34[3],zeta[3],t1[3],t2[3],a32[3],D[3],fact ;
	double b11,b12,b22,E[3],F[3] ;
	double DD1[3],DD2[3],dd1,dd2 ;
	int l2,l3 ;
	ffx = (double *)calloc(Nt,sizeof(double)) ;
	ffy = (double *)calloc(Nt,sizeof(double)) ;
	ffz = (double *)calloc(Nt,sizeof(double)) ;
	for ( int i = 0 ; i < Nt ; ++i )
	{
		ffx[i] = 0 ;
		ffy[i] = 0 ;
		ffz[i] = 0 ;
	}
	#pragma omp for schedule(dynamic,5)
	for ( int i = 0 ; i < edges.size() ; ++i )
	{
				int l1 = edges[i].l1 ;
				int l2 = edges[i].l2 ;
				double l0 = edges[i].l0 ;
				A[0] = (vts[l1].x + vts[l1].ix*Lx) - (vts[l2].x + vts[l2].ix*Lx) ;
				A[1] = (vts[l1].y + vts[l1].iy*Ly) - (vts[l2].y + vts[l2].iy*Ly) ;
				A[2] = vts[l1].z - vts[l2].z ;
				double aa = sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] ) ;
				A[0] /= aa ;
				A[1] /= aa ;
				A[2] /= aa ;
				double xx = aa/(2.5*l0) ;
				double xxx = 1.0 / (4.0*(1.0 - xx)*(1.0 - xx)) ;
				ffx[l1] += (-1.0/(elap))*( xxx - 0.25 + xx )*A[0] ;
				ffy[l1] += (-1.0/(elap))*( xxx - 0.25 + xx )*A[1] ;
				ffz[l1] += (-1.0/(elap))*( xxx - 0.25 + xx )*A[2] ;
				ffx[l2] += (1.0/(elap))*( xxx - 0.25 + xx )*A[0] ;
				ffy[l2] += (1.0/(elap))*( xxx - 0.25 + xx )*A[1] ;
				ffz[l2] += (1.0/(elap))*( xxx - 0.25 + xx )*A[2] ;
					
				ffx[l1] += edges[i].kappap*A[0]/(aa*aa) ;
				ffy[l1] += edges[i].kappap*A[1]/(aa*aa) ;
				ffz[l1] += edges[i].kappap*A[2]/(aa*aa) ;
				ffx[l2] += -edges[i].kappap*A[0]/(aa*aa) ;
				ffy[l2] += -edges[i].kappap*A[1]/(aa*aa) ;
				ffz[l2] += -edges[i].kappap*A[2]/(aa*aa) ;										
	}
	//Surface area and Volume conservation
	#pragma omp for schedule(dynamic,5)
	for ( int i = 0 ; i <  tries.size() ; ++i )
	{
				int l1 = tries[i].l1 ;
				int l2 = tries[i].l2 ;
				int l3 = tries[i].l3 ;
				double beta = -kappasa*(Sarea-Sarea0)/Sarea0 ;
				double betav = -kappav*(Volume - Volume0)/(Volume0) ;  				 
				double area = Areatri(l1,l2,l3) ;
				double alpha = beta/(4.0*area) ;  
				double alphad = -kappad*(area - tarea[i])/(4.0*tarea[i]*area) ;			   
									
				//first forces due to S.area constraint
				A21[0] = (vts[l2].x + vts[l2].ix*Lx) - (vts[l1].x + vts[l1].ix*Lx) ;
				A21[1] = (vts[l2].y + vts[l2].iy*Ly) - (vts[l1].y + vts[l1].iy*Ly) ;
				A21[2] = vts[l2].z - vts[l1].z ;
				A13[0] = (vts[l1].x + vts[l1].ix*Lx) - (vts[l3].x + vts[l3].ix*Lx) ;
				A13[1] = (vts[l1].y + vts[l1].iy*Ly) - (vts[l3].y + vts[l3].iy*Ly) ;
				A13[2] = vts[l1].z - vts[l3].z ;
				A32[0] = (vts[l3].x + vts[l3].ix*Lx) - (vts[l2].x + vts[l2].ix*Lx) ;
				A32[1] = (vts[l3].y + vts[l3].iy*Ly) - (vts[l2].y + vts[l2].iy*Ly) ;
				A32[2] = vts[l3].z - vts[l2].z ;				
				cp(A13,A21,eta) ; 
				//force on first vertex
				cp(eta,A32,C) ;
				ffx[l1] += alpha*C[0] ;
				ffy[l1] += alpha*C[1] ;
				ffz[l1] += alpha*C[2] ;
				ffx[l1] += alphad*C[0] ;
				ffy[l1] += alphad*C[1] ;
				ffz[l1] += alphad*C[2] ;
				cp(eta,A13,C) ;
				ffx[l2] += alpha*C[0] ;
				ffy[l2] += alpha*C[1] ;
				ffz[l2] += alpha*C[2] ;
				ffx[l2] += alphad*C[0] ;
				ffy[l2] += alphad*C[1] ;
				ffz[l2] += alphad*C[2] ;
				cp(eta,A21,C) ;
				ffx[l3] += alpha*C[0] ;
				ffy[l3] += alpha*C[1] ;
				ffz[l3] += alpha*C[2] ;
				ffx[l3] += alphad*C[0] ;
				ffy[l3] += alphad*C[1] ;
				ffz[l3] += alphad*C[2] ;
				B[0] = ((vts[l1].x + vts[l1].ix*Lx) + (vts[l2].x + vts[l2].ix*Lx) + (vts[l3].x + vts[l3].ix*Lx))/3.0  ;
				B[1] = ((vts[l1].y + vts[l1].iy*Ly) + (vts[l2].y + vts[l2].iy*Ly) + (vts[l3].y + vts[l3].iy*Ly))/3.0  ;
				B[2] = (vts[l1].z + vts[l2].z + vts[l3].z)/3.0 ;					
				cp(B,A32,C) ;
				ffx[l1] += betav*((eta[0]/3.0) + C[0])/6.0 ;
				ffy[l1] += betav*((eta[1]/3.0) + C[1])/6.0 ;
				ffz[l1] += betav*((eta[2]/3.0) + C[2])/6.0 ;
				cp(B,A13,C) ;
				ffx[l2] += betav*((eta[0]/3.0) + C[0])/6.0 ;
				ffy[l2] += betav*((eta[1]/3.0) + C[1])/6.0 ;
				ffz[l2] += betav*((eta[2]/3.0) + C[2])/6.0 ;
				cp(B,A21,C) ;
				ffx[l3] += betav*((eta[0]/3.0) + C[0])/6.0 ;
				ffy[l3] += betav*((eta[1]/3.0) + C[1])/6.0 ;
				ffz[l3] += betav*((eta[2]/3.0) + C[2])/6.0 ;			
	}
			
			//Bending force
			#pragma omp for schedule(dynamic,5)
			for (int i = 0 ; i <  edges.size() ; ++i )
			{
				if ( iedge[i] == 1 )
				{
					l2 = edges[i].l1 ;
					l3 = edges[i].l2 ;
				}
				if ( iedge[i] == -1 )
				{
					l3 = edges[i].l1 ;
					l2 = edges[i].l2 ;
				}
				int l1 = edgest[i].l1 ; //These pairs are used for bending energy calculation
				int l4 = edgest[i].l2 ; //
				a21[0] = (vts[l2].x + vts[l2].ix*Lx) - (vts[l1].x + vts[l1].ix*Lx) ;
				a21[1] = (vts[l2].y + vts[l2].iy*Ly) - (vts[l1].y + vts[l1].iy*Ly) ;
				a21[2] = vts[l2].z - vts[l1].z ;
				a31[0] = (vts[l3].x + vts[l3].ix*Lx) - (vts[l1].x + vts[l1].ix*Lx) ;
				a31[1] = (vts[l3].y + vts[l3].iy*Ly) - (vts[l1].y + vts[l1].iy*Ly) ;
				a31[2] = vts[l3].z - vts[l1].z ;
				a24[0] = (vts[l2].x + vts[l2].ix*Lx) - (vts[l4].x + vts[l4].ix*Lx) ;
				a24[1] = (vts[l2].y + vts[l2].iy*Ly) - (vts[l4].y + vts[l4].iy*Ly) ;
				a24[2] = vts[l2].z - vts[l4].z ;
				a34[0] = (vts[l3].x + vts[l3].ix*Lx) - (vts[l4].x + vts[l4].ix*Lx) ;
				a34[1] = (vts[l3].y + vts[l3].iy*Ly) - (vts[l4].y + vts[l4].iy*Ly) ;
				a34[2] = vts[l3].z - vts[l4].z ;
				a32[0] = (vts[l3].x + vts[l3].ix*Lx) - (vts[l2].x + vts[l2].ix*Lx) ;
				a32[1] = (vts[l3].y + vts[l3].iy*Ly) - (vts[l2].y + vts[l2].iy*Ly) ;
				a32[2] = vts[l3].z - vts[l2].z ;
				
				t1[0] = ((vts[l1].x + vts[l1].ix*Lx) + (vts[l2].x + vts[l2].ix*Lx) + (vts[l3].x + vts[l3].ix*Lx))/3.0  ;
				t1[1] = ((vts[l1].y + vts[l1].iy*Ly) + (vts[l2].y + vts[l2].iy*Ly) + (vts[l3].y + vts[l3].iy*Ly))/3.0  ;
				t1[2] = (vts[l1].z + vts[l2].z + vts[l3].z)/3.0 ;
				t2[0] = ((vts[l4].x + vts[l4].ix*Lx) + (vts[l2].x + vts[l2].ix*Lx) + (vts[l3].x + vts[l3].ix*Lx))/3.0 ;
				t2[1] = ((vts[l4].y + vts[l4].iy*Ly) + (vts[l2].y + vts[l2].iy*Ly) + (vts[l3].y + vts[l3].iy*Ly))/3.0  ;
				t2[2] = (vts[l4].z + vts[l2].z + vts[l3].z)/3.0 ;
				cp(a21,a31,eta) ;
				cp(a34,a24,zeta) ;
				double etam = sqrt( eta[0]*eta[0] + eta[1]*eta[1] + eta[2]*eta[2] ) ;
				double zetam = sqrt( zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2] ) ;
				double costh = (eta[0]*zeta[0] + eta[1]*zeta[1] + eta[2]*zeta[2] ) / (etam * zetam) ;
				if ( costh > 1.0 )
				{
					costh = 1.0 ;
				}
				if (costh < -1.0 )
				{
					costh = -1.0 ;
				}		
				if ( costh == 1.0 )
				{
					costh -= 0.0005 ;
				}
				if ( costh == -1.0 )
				{
					costh += 0.0005 ;
				}	
				double sinth = sqrt ( 1.0 - costh*costh ) ;
				double tt0 = t1[0]-t2[0] ;
				double tt1 = t1[1]-t2[1] ;
				double factor  = ( ((eta[0]-zeta[0])*tt0) + ((eta[1]-zeta[1])*tt1) + ((eta[2]-zeta[2])*(t1[2]-t2[2])) ); 
				if ( factor >= 0.0 )
				{
					fact = 1.0 ;
				}
				else
				{
					fact = -1.0 ; 
				}
				sinth = fact*sinth ;
				if ( sinth == 0.0 )
				{
					cout << "error bend " << endl ;
				}
				double betab = kappab*(sinth*cos(theta0*Pi/180.0) - costh*sin(theta0*Pi/180.0))/sinth ;
				b11 = -(betab*costh)/(etam*etam) ;
				b22 = -(betab*costh)/(zetam*zetam) ;
				b12 = betab/(etam*zetam) ;
				cp(eta,a32,C) ;
				cp(zeta,a32,D) ;
				ffx[l1] += b11*C[0] + b12*D[0] ;
				ffy[l1] += b11*C[1] + b12*D[1] ;
				ffz[l1] += b11*C[2] + b12*D[2] ;
				ffx[l4] += -b12*C[0] - b22*D[0] ;
				ffy[l4] += -b12*C[1] - b22*D[1] ;
				ffz[l4] += -b12*C[2] - b22*D[2] ;
				cp(eta,a31,C) ;
				cp(zeta,a31,D) ;
				cp(eta,a34,E) ;
				cp(zeta,a34,F) ;
				ffx[l2] += -b11*C[0] + b12*(E[0]-D[0]) + b22*F[0] ;
				ffy[l2] += -b11*C[1] + b12*(E[1]-D[1]) + b22*F[1] ;
				ffz[l2] += -b11*C[2] + b12*(E[2]-D[2]) + b22*F[2] ;
				cp(eta,a21,C) ;
				cp(zeta,a21,D) ;
				cp(eta,a24,E) ;
				cp(zeta,a24,F) ;
				ffx[l3] += b11*C[0] + b12*(-E[0]+D[0]) - b22*F[0] ;
				ffy[l3] += b11*C[1] + b12*(-E[1]+D[1]) - b22*F[1] ;
				ffz[l3] += b11*C[2] + b12*(-E[2]+D[2]) - b22*F[2] ;
			}
			#pragma omp critical
			{
				for(int i = 0 ; i < Nt ; ++i )
				{
					fmapdt[i].fx += ffx[i] ;
					fmapdt[i].fy += ffy[i] ;
					fmapdt[i].fz += ffz[i] ;
				}
			}
			free(ffx) ;
			free(ffy) ;
			free(ffz) ;
   }
} 

void pfdt()
{
	double cutoff = paraR ;
	double epsilon = 0.1 ;
	double sigma = 0.890898718*cutoff ;
	#pragma omp parallel default(shared)
	{
		double fx,fy,fz ;
		double A[3],X[3] ;
		X[0] = X[1] = X[2] = 0.0 ;
		fx = 0 ;
		fy = 0 ;
		fz = 0 ;
		#pragma omp for
		for ( int i = 0 ; i < vts.size() ; ++i )
		{
			A[0] = (para[0] + paradx*Lx) - (vts[i].x + vts[i].ix*Lx) ;
			A[1] = (para[1] + parady*Ly) - (vts[i].y + vts[i].iy*Ly) ;
			A[2] = para[2] - vts[i].z ;
			double aa =  A[0]*A[0] + A[1]*A[1] + A[2]*A[2]  ;
			if ( aa <= cutoff*cutoff)
			{
				double ex = (sigma * sigma * sigma * sigma * sigma * sigma)/(aa * aa * aa);
				 X[0] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[0];
				 X[1] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[1]; 
				 X[2] = (24.0 * epsilon/aa) * ((2.0 * ex * ex) - ex ) * A[2];
				 fx += X[0] ;
				 fy += X[1] ;
				 fz += X[2] ;
				 fmapdt[i].fx -= X[0] ;
				 fmapdt[i].fy -= X[1] ;
				 fmapdt[i].fz -= X[2] ;
			}	
		}
		#pragma omp critical
		{
			parafdt[0] += fx ;
			parafdt[1] += fy ;
			parafdt[2] += fz ;
		}
	} //End of parallel region
}
	


double bendenergy()
{
	double energy = 0.0 ;
	double eta[3],a21[3],a31[3],a24[3],a34[3],zeta[3],t1[3],t2[3],a32[3] ;
	for ( int i = 0 ; i < edges.size() ; ++i )
	{
			int l2 = edges[i].l1 ;
			int l3 = edges[i].l2 ;
			int l1 = edgest[i].l1 ; //These pairs are used for bending energy calculation
			int l4 = edgest[i].l2 ; //
			a21[0] = vts[l2].x - vts[l1].x ;
			a21[1] = vts[l2].y - vts[l1].y ;
			a21[2] = vts[l2].z - vts[l1].z ;
			a31[0] = vts[l3].x - vts[l1].x ;
			a31[1] = vts[l3].y - vts[l1].y ;
			a31[2] = vts[l3].z - vts[l1].z ;
			a24[0] = vts[l2].x - vts[l4].x ;
			a24[1] = vts[l2].y - vts[l4].y ;
			a24[2] = vts[l2].z - vts[l4].z ;
			a34[0] = vts[l3].x - vts[l4].x ;
			a34[1] = vts[l3].y - vts[l4].y ;
			a34[2] = vts[l3].z - vts[l4].z ;
			a32[0] = vts[l3].x - vts[l2].x ;
			a32[1] = vts[l3].y - vts[l2].y ;
			a32[2] = vts[l3].z - vts[l2].z ;
			t1[0] = (vts[l1].x + vts[l2].x + vts[l3].x)/3.0 ;
			t1[1] = (vts[l1].y + vts[l2].y + vts[l3].y)/3.0 ;
			t1[2] = (vts[l1].z + vts[l2].z + vts[l3].z)/3.0 ;
			t2[0] = (vts[l4].x + vts[l2].x + vts[l3].x)/3.0 ;
			t2[1] = (vts[l4].y + vts[l2].y + vts[l3].y)/3.0 ;
			t2[2] = (vts[l4].z + vts[l2].z + vts[l3].z)/3.0 ;
			cp(a21,a31,eta) ;
			cp(a34,a24,zeta) ;
			double etam = sqrt( eta[0]*eta[0] + eta[1]*eta[1] + eta[2]*eta[2] ) ;
			double zetam = sqrt( zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2] ) ;
			double costh = (eta[0]*zeta[0] + eta[1]*zeta[1] + eta[2]*zeta[2] ) / (etam * zetam) ;
			energy += kappab*(1.0 - costh*costh) ;
	}
	return energy ;
}
		




void makeiedge()
{

	double a21[3],a31[3],a34[3],a24[3],t1[3],t2[3],tt1[3],tt2[3],eta[3],zeta[3] ;
	for (int i = 0 ; i <  edges.size() ; ++i )
	{
	    
		int l2 = edges[i].l1 ;
		int l3 = edges[i].l2 ;
		int l1 = edgest[i].l1 ; //These pairs are used for bending energy calculation
		int l4 = edgest[i].l2 ; //
		a21[0] = vts[l2].x - vts[l1].x ;
		a21[1] = vts[l2].y - vts[l1].y ;
		a21[2] = vts[l2].z - vts[l1].z ;
		a31[0] = vts[l3].x - vts[l1].x ;
		a31[1] = vts[l3].y - vts[l1].y ;
		a31[2] = vts[l3].z - vts[l1].z ;
		a24[0] = vts[l2].x - vts[l4].x ;
		a24[1] = vts[l2].y - vts[l4].y ;
		a24[2] = vts[l2].z - vts[l4].z ;
		a34[0] = vts[l3].x - vts[l4].x ;
		a34[1] = vts[l3].y - vts[l4].y ;
		a34[2] = vts[l3].z - vts[l4].z ;
		t1[0] = (vts[l1].x + vts[l2].x + vts[l3].x)/3.0 ;
		t1[1] = (vts[l1].y + vts[l2].y + vts[l3].y)/3.0 ;
		t1[2] = (vts[l1].z + vts[l2].z + vts[l3].z)/3.0 ;
		t2[0] = (vts[l4].x + vts[l2].x + vts[l3].x)/3.0 ;
		t2[1] = (vts[l4].y + vts[l2].y + vts[l3].y)/3.0 ;
		t2[2] = (vts[l4].z + vts[l2].z + vts[l3].z)/3.0 ;
		tt1[0] = t1[0] - center[0] ;
		tt1[1] = t1[1] - center[1] ;
		tt1[2] = t1[2] - center[2] ;
		tt2[0] = t2[0] - center[0] ;
		tt2[1] = t2[1] - center[1] ;
		tt2[2] = t2[2] - center[2] ;
		cp(a21,a31,eta) ;
		cp(a34,a24,zeta) ;
		double che1 = eta[0]*tt1[0] + eta[1]*tt1[1] + eta[2]*tt1[2] ;
		double che2 = zeta[0]*tt2[0] + zeta[1]*tt2[1] + zeta[2]*tt2[2] ;
		if (( che1 >= 0) && (che2 >= 0 ))
		{
			 iedge[i] = 1 ;
		}
		else if (( che1 < 0) && (che2 < 0 ))
		{
			 iedge[i] = -1 ;
		}
		else
		{
		    cout << "error making iedge: Aborting.." << endl ;
		    //abort() ;
		    iedge[i] = 1 ;
		}		

	}
}

void filewrite(char filename[])
{
			double N[3] ;
			string s0 = "facet normal" ;
			string s1 = "outer loop" ;
			string s2 = "vertex" ;
			string s3 = "endloop";
			string s4 = "endfacet" ;
			ofstream myfile ;
			myfile.open(filename) ;
			myfile << "solid ascii" << endl ;
			for ( unsigned int kk = 0 ; kk < tries.size() ; ++kk )
			{
				int l1 = tries[kk].l1 ;
				int l2 = tries[kk].l2 ;
				int l3 = tries[kk].l3 ;
				Normaldt(l1,l2,l3,N) ;
				double f0 = vts[l1].x + vts[l1].ix*Lx ;
				double f1 = vts[l1].y + vts[l1].iy*Ly;
				double f2 = vts[l1].z ;
				double f3 = vts[l2].x + vts[l2].ix*Lx ;
				double f4 = vts[l2].y + vts[l2].iy*Ly;
				double f5 = vts[l2].z ;
				double f6 = vts[l3].x + vts[l3].ix*Lx;
				double f7 = vts[l3].y + vts[l3].iy*Ly;
				double f8 = vts[l3].z ;
				myfile << s0 << "\t" << N[0] << "\t" << N[1] << "\t" << N[2] << endl ;
				myfile << s1 << endl ;
				myfile << s2 << "\t" << f0 << "\t" << f1 << "\t" << f2 << endl ;
				myfile << s2 << "\t" << f3 << "\t" << f4 << "\t" << f5 << endl ;
				myfile << s2 << "\t" << f6 << "\t" << f7 << "\t" << f8 << endl ;
				myfile << s3 << endl ;
				myfile << s4 << endl ;
			}
			myfile << "endsolid" << endl ;
			myfile.close() ;
}


