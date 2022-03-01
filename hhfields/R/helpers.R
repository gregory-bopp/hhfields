#' Matern kernel generator
#' @param ep scale parameter
#'
#' @export
matern_kerngen <- function(ep){
  rbf <- function(ep, r){
    x <- ep*r
    exp(-x)*(x^5 + 15*x^4 + 105*x^3 + 420*x^2 + 945*x + 945)/945
  }

  drbf <- function(xp, r){
    x <- ep*r
    -(ep^2)*exp(-x)*(x^4 + 10*x^3 + 45*x^2 + 105*x + 105)/945
  }

  d2rbf <- function(xp, r){
    x <- ep*r
    (ep^4)*exp(-x)*(x^3 + 6*x^2 +15*x +15)/945
  }
  kerngen <- list(rbf = rbf,
                  drbf = drbf,
                  d2rbf = d2rbf,
                  ep = ep)
}

#' Simulate Domain 
#' @description Simulate a domain over which to do a curl-free and divergence-
#' free dcomposition of a vector field. The domain is simulated by removing 
#' polygon subsets of the interior of an outer polygon. A number of interior
#' squares are generated and intersected with each other to define the removal 
#' region from the outer polygon.
#' @param outer_poly (sf polygon) outer polygon. 
#'
#' @param n_isq number of inner squares to use
#' @param n_bndry Number of outer polygon boundary nodes to sample 
#' @param n_inner Number of inner polygon boundary nodes to sample
#'
#' @export
sim_domain <- function(outer_poly, n_isq, n_bndry, n_inner){
  # Simmulate inner squares
  inner_centers <- st_sample(outer_poly, n_isq, type = "random",
                             exact = T)
  max_dist <- max(calc_distC(st_coordinates(outer_poly)))
  ipolys <- map(1:n_isq,
                ~{
                  make_sq(runif(1, 0, max_dist/10))*rot(runif(1, 0, pi/2)) +
                    st_point(runif(2, -max_dist/2, max_dist/2))
                }) %>%
    st_sfc() %>%
    st_intersection(outer_poly,.x) %>%
    st_union

  # Make exterior square ----------------------------------------------------
  poly_list <- st_geometry(ipolys) %>%
                  st_cast("POLYGON")
  # Remove holes that may have formed from intersecting polygons
  poly_list <- map(poly_list, ~ st_polygon(lapply(st_geometry(.x),
                                                  function(x) x[1])[[1]]))%>%
    st_sfc()

  bndry_nodes <- map(poly_list, ~sample_poly_points(n_bndry, .x, remove_vert = T))
  norm_pts <- map2(poly_list, bndry_nodes,
                   ~get_poly_normals(.x, .y, interior = F)$norm %>%
                     st_sfc())

  trans_norm <- map2(poly_list, bndry_nodes,
                     ~get_poly_normals(.x, .y, interior = F)$trans_norm %>%
                       st_sfc())
  bndry_nodes <- map(bndry_nodes, ~ .x[[1]] %>%  st_sfc())
  trans_norm <- Reduce(c,trans_norm)
  norm <- Reduce(c, norm_pts)
  bndry_nodes <- Reduce(c, bndry_nodes)


# Sample inner nodes ------------------------------------------------------
  inner_diff <-st_difference(outer_poly, ipolys)
  inner_nodes <-  st_sample(inner_diff, n_inner, type = "hexagonal",
                    exact = T)

  return(list(norm = norm,
              trans_norm = trans_norm,
              bndry_nodes = bndry_nodes,
              inner_poly = poly_list,
              inner_nodes = inner_nodes))
}
