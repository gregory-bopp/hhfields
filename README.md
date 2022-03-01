# Helmholtz-Hodge Decomposition

This repository is an R implementation of the Helmhotz-Hodge decomposition of a vector field defined on a bounded domain described in [Fuselier and Wright, 2016](https://arxiv.org/abs/1502.01575). 

[Fuselier and Wright, 2016](https://arxiv.org/abs/1502.01575) propose a radial basis function method using matrix-valued kernels to compute the Helmholtz-Hodge decomposition on a bounded domain using samples of the field at a collection of nodes.

This repository applies their technique to bounded domains constructed from intersecting polygons. This approach may be useful for doing spatial inference on irregularly sampled solenoidal vector fields (e.g. magnetic fields or the velocity field of an incompressible fluid) with boundary constraints.

General Helmholtz-Hodge background [here](./doc/details.pdf)

# Example Simulation

![](analysis/fig/hh_multipoly.png)