# ButterflyLab

1 OVERVIEW

ButterflyLab is a software package in MATLAB and c++ containing various algorithms for nearly optimal fast matvec and dense linear system solvers for (hierarchical) complementary low-rank matrices (see the definition in [5]) based on the butterfly algorithm [1-4], (hierarchical) butterfly factorization [5-9], and triangular solvers [10]. Famous examples of complementary low-rank matrices are the matrix representation of Fourier integral operators (FIOs) for wave equations and a class of imaging methods, a class of transforms for harmonic analysis and special functions, high-frequency integral equations for scattering problems.

2 INTRODUCTION

2.1 The first version of the ButterflyLab contains the following codes and examples in [1-7].

2.1.1 MBA.m

Multiscale Butterfly Algorithm is a algorithm designed for fast evaluation of Fourier Integral Operators. Both 2D and 3D code are provided. Examples are in test folders.

2.1.2 BF.m

Contains various examples of the butterfly factorization in 1D and 2D as follows.

1D Butterfly Factorization

Butterfly factorization for all 1D problems.

2D Butterfly Factorization

General Butterfly Factorization (GBF)

Butterfly factorization for 2D problems without singularity at origin.

Polar Butterfly Factorization (PBF)

Butterfly factorization for 2D problems with singularity at origin.

Multiscale Butterfly Factorization (MBF)

Butterfly factorization for 2D problems with singularity at origin.

2.1.3 IBF.m

Contains various examples of the Interpolative Butterfly Factorization [7].

The current implementation supports any dimension
Fourier integral operators with/without singularity at the origin.

Fourier integral operatiors without singularity could be solved via `IBF` in the `src` folder, whereas the Fourier integral operators with singularity at the origin could be solved via either `fastPBF` with polar transform or `fastMBF`. The former adopts fast butterfly factorization with the idea given in [3] and the later adopts the idea of multiscale domain decomposition [4] together
with the interpolative butterfly factorization.

1) The example for `fastBF` can be found in `IBF/test/test_fastbf_1D` and `IBF/test/test_fastbf_2D`,

2) The example for `fastPBF` with polar transform can be found in `IBF/test/test_fastpbf_2D`,

3) The example for `fastMBF` can be found in `IBF/test/test_fastmbf_2D`.

More examples of special function transforms can be found in `test` folder as well.

2.2 The second version of the ButterflyLab reorganize codes in the first version and contains more codes and examples in [8-10].

3 REFERENCES:

[1] Eric Michielssen and Amir Boag, MULTILEVEL EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994. [PDF][ba]

[2] E. Candes, L. Demanet and L. Ying. A fast butterfly algorithm for the computation of Fourier integral operators. SIAM Multiscale Modeling and Simulation 7 (2009). [PDF][pbf]

[3] L. Ying. Sparse Fourier transform via butterfly algorithm. SIAM J. Sci. Comput., 31(3):1678–1694, Feb. 2009. [PDF][spF]

[4] Y. Li, H. Yang, and L. Ying. A multiscale butterfly aglorithm for Fourier integral operators. Multiscale Modeling and Simulation, 13(2):614–631, 2015. [PDF][mba]

[5] Y. Li, H. Yang, E. R. Martin, K. L. Ho, and L. Ying. Butterfly Factorization. Multiscale Modeling & Simulation, 13(2):714–732, 2015. [PDF][bf]

[6] Y. Li, H. Yang, and L. Ying. Multidimensional butterfly factorization. Applied and Computational Harmonic Analysis, Volume 44, Issue 3, May 2018, Pages 737-758. [PDF][mbf]

[7] Y. Li and H. Yang. Interpolative Butterfly Factorization. SIAM J. Sci. Comput., 39(2), A503–A531, 2017. [PDF][ibf]

[8] H. Yang. A Unified Framework for Oscillatory Integral Transform: When to Use NUFFT or Butterfly factorization? Journal of Computational Physics, 2019. [PDF][uni]

[9] Q. Pang, K. L. Ho, and H. Yang. Interpolative Decomposition Butterfly Factorization. SIAM Journal of Scientific Computing, 2020. [PDF][idbf]

[10] Y. Liu, and H. Yang. A Hierarchical Butterfly LU Preconditioner for Two-Dimensional Electromagnetic Scattering Problems Involving Open Surfaces. Journal of Computational Physics, 2020. [PDF][HIBLU]

[11] Z. Chen, J. Zhang, K. Ho, and H. Yang. Multidimensional Phase Recovery and Interpolative Decomposition Butterfly Factorization. Journal of Computational Physics, 2020. [PDF][MIDBF]

[12] J. Bremer, Z. Chen, and H. Yang. Rapid Application of the Spherical Harmonic Transform via Interpolative Decomposition Butterfly Factorization. Submitted. [PDF][MIDBF]

[ba]: https://onlinelibrary.wiley.com/doi/abs/10.1002/mop.4650071707
[pbf]: http://epubs.siam.org/doi/abs/10.1137/080734339
[spF]: http://epubs.siam.org/doi/abs/10.1137/08071291X
[mba]: http://epubs.siam.org/doi/abs/10.1137/140997658
[bf]: http://epubs.siam.org/doi/abs/10.1137/15M1007173
[ibf]: https://epubs.siam.org/doi/abs/10.1137/16M1074941
[mbf]: https://www.sciencedirect.com/science/article/abs/pii/S1063520317300271
[uni]: https://www.sciencedirect.com/science/article/pii/S0021999119301615
[idbf]: https://epubs.siam.org/doi/abs/10.1137/19M1294873?mobileUi=0
[HIBLU]: https://www.sciencedirect.com/science/article/pii/S002199911930720X
[MIDBF]: https://www.sciencedirect.com/science/article/pii/S0021999120302011
[SHT]: https://www.researchgate.net/publication/340859249_Rapid_Application_of_the_Spherical_Harmonic_Transform_via_Interpolative_Decomposition_Butterfly_Factorization?channel=doi&linkId=5ea1028c458515ec3aff1ed2&showFulltext=true

4 COPY RIGHT

ButterflyLab is copyright reserved. For further information, please contact main contributors below.
Contact information
Ze Chen at chenze@u.nus.edu
Kenneth L. Ho at klho@alumni.caltech.edu
Yingzhou Li at yingzhou.li@duke.edu
Qiyuan Pang at ppangqqyz@foxmail.com
Haizhao Yang at haizhao@purdue.edu
Lexing Ying at lexing@math.stanford.edu
