# Helmholtz-Hodge Decomposition

The Helmholtz-Hodge decomposition says that a twice-differentiable vector field **F(r)** on a bounded domain $V \subset \mathbb{R}^3$ can be expressed as the sum of a curl-free component $\mathbf{F_l(r)}$ and  divergence-free component $\mathbf{F_t(r)}$ component:

$$
\mathbf{F(r)} = \mathbf{F_l(r)} + \mathbf{F_t(r)}
$$
Namely, $\mathbf{F_l(r)}$ is the longitudinal or irrotational (curl-free) part of vector:

$$
\nabla \times \mathbf{F_l(r) = 0}.
$$

and $\mathbf{F_t(r)}$ is the solenoidal part (divergence-free, i.e. no sources or sinks):

$$
\nabla \cdot \mathbf{F_t(r)} = 0.
$$


The irrotational part can be expressed as the negative gradient of a scalar potential function $\phi(\mathbf{r})$:

$$
\mathbf{F_l(r) = -\nabla \phi(\mathbf{r})}.
$$


The solenoidal part can be expressed as the curl of a vector potential $\mathbf{a(r)}$

$$
\mathbf{F_t(r)} = \nabla \times \mathbf{a(r)}
$$


[Fuselier and Wright, 2016](https://arxiv.org/abs/1502.01575) propose a radial basis function method using matrix-valued kernels to compute the Helmholtz-Hodge decomposition on a bounded domain using samples of the field at a collection of nodes.

This repository applies their technique to bounded domains constructed from intersecting polygons in 2-D. This approach may be useful for doing spatial inference on irregularly sampled solenoidal vector fields (e.g. magnetic fields or the velocity field of an incompressible fluid) with boundary constraints

# Example Simulation

![](../analysis/fig/hh_multipoly.png)