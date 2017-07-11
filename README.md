# ButterflyLab

1 OVERVIEW

ButterflyLab is a software package in MATLAB and c++ containing various algorithms for fast matvec for complementary low-rank matrices (see the definition in [5]) based on the butterfly algorithm [1-4] and butterfly factorization [5-7]. Famous examples of complementary low-rank matrices are the matrix representation of Fourier integral operators (FIOs) for wave equations and a class of imaging methods, and a class of transforms for harmonic analysis and special functions. 

2 INTRODUCTION

2.1 MBA.m

Multiscale Butterfly Algorithm is a algorithm designed for fast evaluation of Fourier Integral Operators. Both 2D and 3D code are provided. Examples are in test folders.

2.2 BF.m

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

2.3 IBF.m

Contains various examples of the Interpolative Butterfly Factorization [7].

The current implementation supports any dimension
Fourier integral operators with/without singularity at the origin.

Fourier integral operatiors without singularity could be solved via `IBF` in the `src` folder, whereas the Fourier integral operators with singularity at the origin could be solved via either `fastPBF` with polar transform or `fastMBF`. The former adopts fast butterfly factorization with the idea given in [3] and the later adopts the idea of multiscale domain decomposition [4] together
with the interpolative butterfly factorization.

1) The example for `fastBF` can be found in `IBF/test/test_fastbf_1D` and `IBF/test/test_fastbf_2D`,

2) The example for `fastPBF` with polar transform can be found in `IBF/test/test_fastpbf_2D`,

3) The example for `fastMBF` can be found in `IBF/test/test_fastmbf_2D`.

More examples of special function transforms can be found in `test` folder as well.


3 REFERENCES:

[1] E. Michielssen and A. Boag. A multilevel matrix decomposition algorithm for analyzing scattering from large structures. Antennas and Propagation, IEEE Transactions on, 44(8):1086–1093, Aug 1996.. [PDF][ba]

[2] E. Candes, L. Demanet and L. Ying. A fast butterfly algorithm for the computation of Fourier integral operators. SIAM Multiscale Modeling and Simulation 7 (2009). [PDF][pbf]

[3] L. Ying. Sparse Fourier transform via butterfly algorithm. SIAM J. Sci. Comput., 31(3):1678–1694, Feb. 2009. [PDF][spF]

[4] Y. Li, H. Yang, and L. Ying. A multiscale butterfly aglorithm for Fourier integral operators. Multiscale Modeling and Simulation, 13(2):614–631, 2015. [PDF][mba]

[5] Y. Li, H. Yang, E. R. Martin, K. L. Ho, and L. Ying. Butterfly Factorization. Multiscale Modeling & Simulation, 13(2):714–732, 2015. [PDF][bf]

[6] Y. Li, H. Yang, and L. Ying. Multidimensional butterfly factorization. Applied and Computational Harmonic Analysis, in press. [PDF][mbf]

[7] Y. Li and H. Yang. Interpolative Butterfly Factorization. SIAM J. Sci. Comput., 39(2), A503–A531, 2017. [PDF][ibf]

[ba]: http://ieeexplore.ieee.org/document/511816/
[pbf]: http://epubs.siam.org/doi/abs/10.1137/080734339
[spF]: http://epubs.siam.org/doi/abs/10.1137/08071291X
[mba]: http://epubs.siam.org/doi/abs/10.1137/140997658
[bf]: http://epubs.siam.org/doi/abs/10.1137/15M1007173
[ibf]: http://arxiv.org/abs/1605.03616
[mbf]: http://arxiv.org/abs/1509.07925

4 COPY RIGHT

ButterflyLab is copyright reserved. For further information, please contact 
Contact information
Yingzhou Li at ryanlee@stanford.edu
Haizhao Yang at haizhao@math.duke.edu
Lexing Ying at lexing@math.stanford.edu
