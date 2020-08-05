# Introduction
We provide here a Matlab-based PolyCompositeFEM3D for three-dimensional solid mechanics probmems. We propose a new formulation 
based on arbitrary polytopes, which we name as Polytopal Composite finite Element (PCE). The key idea is to 
devise a polynomial projection of compatible strain fields through the least-squares approximation. 
For nearly incompressible problems, we design a pair of projection operators of volumetric and deviatoric 
strains resulting in stability of pressure solution. The present approach fulfills a patch test and satisfies 
the inf–sup stability. Through several numerical investigations, we show that the proposed method reaches 
the theoretical convergence rate and significantly improves the accuracy of polygonal element based solutions. 

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

# References:
H. Nguyen-Xuan, Khanh N. Chau, Khai N. Chau, Polytopal composite finite elements, Computer Methods in Applied Mechanics and Engineering, 355, 405-437, 2019. https://doi.org/10.1016/j.cma.2019.06.030
