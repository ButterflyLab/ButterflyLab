MultiDimJacobi
==============
1 INTRODUCTION
-----------------
MultiDimJacobi is a package in Matlab and Fortran and fastly computes multidimensional Jacobi polynomial 
transform based on non-oscillatory phase function method [[1]](https://arxiv.org/abs/1803.03889) (for 1D transform) and a modified one [[2]](https://arxiv.org/abs/1901.07275) (for 2D and 3D transform).
Three kinds of the transform are covered in the package, including inverse and forward uniform 
transform, and forward nonuniform transform. The fast inverse in the nonuniform case is still an open problem. It involves a highly illconditioned linear system of equations. We are working on an efficient preconditioner for this
linear system.

2 CONTENTS
----------------------
This section introduces the following codes computing 1D, 2Dand 3D transform.

- 1D/:  contains three files, JPT1D.m, invJPT1D.m and NJPT1D.m, which compute 1D forward 
uniform, inverse uniform and forward nonuniform, respectively, based on the non-oscillatory phase 
function method in [[1]](https://arxiv.org/abs/1803.03889). 

- 2D/:  contains three files, JPT2D.m, invJPT2D.m and NJPT2D.m, which compute 2D forward
 uniform, inverse uniform and forward nonuniform, respectively. These codes are based on the modified non-oscillatory phase function method in [[2]](https://arxiv.org/abs/1901.07275).

- 3D/:  contains three files, JPT3D.m, invJPT3D.m and NJPT3D.m, which compute 3D forward
 uniform, inverse uniform and forward nonuniform, respectively.These codes are based on the modified non-oscillatory phase function method in [[2]](https://arxiv.org/abs/1901.07275).

- src/: contains the core files for this package 

- test/: contains test codes and comparison codes  

         - test/1D/: contains testJPT1D.m, testinvJPT1D.m, and NJPT1D.m to test JPT1D.m, invJPT1D.m, and NJPT1D.m respectively.
         - test/2D/: contains testJPT2D.m, testinvJPT2D.m, and NJPT2D.m to test JPT2D.m, invJPT2D.m, and NJPT2D.m respectively.
         - test/3D/: contains testJPT3D.m, testinvJPT3D.m, and NJPT3D.m to test JPT3D.m, invJPT3D.m, and NJPT3D.m respectively.
         - test/Comp/: contains three files, Comp_RSFFT_CHEBNF_1D.m, Comp_RSFFT_CHEBNF_2D.m, and Comp_RSFFT_CHEBNF_3D.m to compare the methods in  [1] and [2] for 1D, 2D, and 3D transforms respectively. For 1D transform, we recommend the non-oscillatory phase function method in [1], and for 2D and 3D transforms, we recommend the modified method in [2].

- MultiDimJacobi_startup.m: startup file

3 ARGUMENT SETTING
-------
For polynomial parameters da and db in all files mentioned above, we recommend users to set them in the open interval (-1,1) except -0.5. Our package would break down if da or db was set to -0.5. Values beyond (-1,1) are 
not be recommended either, since our package would be less accurate and even fail when da or db goes beyond (-1,1).
For other arguments, please refer to the annotation in the codes.

4 INSTALLATION
---------------------
Before using the package, you must compile the mex files in /src first. To do this, you can enter
 /src and then type "make" in your cmd to compile all the mex files.
The flags in /src/make.inc are set as follows by default:

###### #define MEX, in mac mex may also refer to pdftex, a command conflict    
MEX = /usr/local/MATLAB/R2018a/bin/mex    
###### #Tested the flag on R2018a    
FLAGS =    
FORTRAN   = gfortran   
OPTS      = -O3   
FPPFLAGS  = -DLAP -DSLAP -fPIC   
FOPTS     = -w -Ofast -cpp   
LIBNAME = -lfftw3 -lblas   

You can change the flags if necessary.
Several test files are provided in  /test for you to test whether you install the package successfully
 and serve as demos to show the usage of our codes.

5 REFERENCES
----
[1]. James Bremer and Haizhao Yang. Fast algorithms for Jacobi expansions via nonoscillatory 
phase functions. arXiv:1803.03889 [math.NA], 2018. [pdf](https://arxiv.org/abs/1803.03889)

[2]. James Bremer, Qiyuan Pang, and Haizhao Yang. Fast Algorithms for the Multi-dimensional 
Jacobi Polynomial Transform. arXiv:1901.07275 [math.NA], 2019. [pdf](https://arxiv.org/abs/1901.07275)

6 COPYRIGHT
----
MultiDimJacobi is copyright reserved. For further information, please contact main contributors below.
 Contact information 
James Bremer at [bremer@math.ucdavis.edu](bremer@math.ucdavis.edu)  Qiyuan Pang at [ppangqqyz@foxmail.com ](ppangqqyz@foxmail.com)  Haizhao Yang at [haizhao@nus.edu.sg](haizhao@nus.edu.sg)