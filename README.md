This is adhesive dynamics code for simulating infected red blood cells (iRBCs) in shear flow
written by Dr. Anil Kumar Dasanna when working with Prof. Dr. Ulrich Schwarz at the 
Institute for Theoretical Physics at Heidelberg, Germany.

The solvent is modelled using Multi-Particle Collision Dynamics (MPCD) and iRBCs are modelled
with triangulated surfaces with elasticity as well as total surface area and total volume constraints.  

Usage:
make clean ; make
./MAIN on-rate off-rate vertices_file edges_list triangle_list Number_knobs rbc_data_file stl_flag

Arguments:
onrate --> onrate of bond dynamics between rbc and substrate
offrate --> offrate of bond dynamics between rbc and substrate
vertices_file --> vertices list of triangulated surface
edges_list --> edge list of triangulated surface
triangle_list --> triangles list of triangulated surface
Number_knobs --> Number of adhesive patches or knobs or receptors(must not exceed the total number of vertices)
rbc_data_file --> file to write RBC vertex coordinates for every few thousand timesteps, format: x y z ix iy bond_id (bond_id indicates whether the vertex is bonded or not)
stl_flag --> either 0 (NO) or 1 (YES) to produce stl files of RBC for every few thousand timesteps

Parallelization:
The code also uses openmp threading. So provide an appropriate statement regarding the number of threads before running. For example:
export OMP_NUM_THREADS=4.

Other parameters:
To change the shear rate, change "walls" in define.h. The shear rate would be walls/Lz.
If you want to use different shape such as discocyte, you need to have appropriate vertices, edges, triangles list files and you need to change 
"Nt" value in define.h which indicates the number of vertices. For the format of these files, please have a look at the files provided in the folder.

Example: 
For the data files provided here, below is the list of commands to simulate rolling adhesion of trophozoite cell in shear flow:

make clean; make
export OMP_NUM_THREADS=8
./MAIN 0.005 0.015 TrophV600.dat TrophE600.dat TrophT600.dat 500 RBCdata.dat 1
