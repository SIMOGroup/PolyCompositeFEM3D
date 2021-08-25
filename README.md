# Introduction
A Matlab-based PolyCompositeFEM3D for three-dimensional solid mechanics probmems was released. We proposed a new formulation 
based on arbitrary polytopes, which we name as Polytopal Composite finite Element (PCE). The key idea is to 
devise a polynomial projection of compatible strain fields through the least-squares approximation. 
For nearly incompressible problems, we designed a pair of projection operators of volumetric and deviatoric 
strains resulting in stability of pressure solution. The present approach satisfies the inf–sup stability. 
Through several numerical investigations, the proposed method reached the theoretical convergence rate 
and significantly improved the solution accuracy of polygonal finite elements. 

This Matlab codes can be extended to a wide range of engineering problems. 

1. Structure of PolyCompositeFEM3D package: 
- mainMainPCEnPrismaticBeamFanTri.m & MainPCEnPrismaticBeamSundarMesh.m: the main function for running three-dimensional solid mechanics probmems.
- Other functions are given in subfolders.
2. How to run PolyCompositeFEM3D: 
- You need to define input of new problem following code structures.  
- Run main....m and then use Paraview to show results in displacements and stresses.
- Get output

# Contributors
- Hung Nguyen-Xuan
- Khanh N. Chau
- Khai N. Chau

# References
H. Nguyen-Xuan, Khanh N. Chau, Khai N. Chau, Polytopal composite finite elements, Computer Methods in Applied Mechanics and Engineering, 355, 405-437, 2019 https://doi.org/10.1016/j.cma.2019.06.030.
