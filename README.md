# log-det

Matlab code that can be used to reproduce the results in the paper "Preconditioned log-determinant approximation: one Gaussian probe vector is almost always enough" by Alice Cortinovis and Daniele Toni. https://arxiv.org/abs/2601.05778

- ''tools'' contains all the technical tools introduced in the paper: all the matrices described in Section 5, an implementation of the Nyström approximation, of the Lanczos method for the quadratic form involving the logarithm and of the Stochastic Lanczos quadrature, and the log-det-ective (Algorithm 4.1)

- any x_comparison file computes the experiments explained in the corresponding Section

- any x_comparison folder contains the saved datas of the corresponding plotted experiment
