###########################################################
# Description: Create an sf exterior square, with normals 
#   pointed towards square interior and an interior triangle
#   with normals pointed away from triangle interior
##########################################################

library(sf)
library(tidyverse)
library(here)
devtools::load_all(here('hhfields'))


# Parameters --------------------------------------------------------------
npt <- 11

# Make exterior square ----------------------------------------------------
sq <- make_sq(4)
sq_pts <- sample_poly_points(npt, sq, remove_vert = T)
sq_norm_ls <- get_poly_normals(sq, sq_pts,interior = T)$trans_norm


# Make interior polygon ---------------------------------------------------
tri <- st_polygon(x = list(cbind(x = c(1,3,3,1),
                                 y = c(1,3,1,1)))) - st_point(c(3,3))
tri_points <- sample_poly_points(npt, tri, remove_vert = T)
norm_ls <- get_poly_normals(tri, tri_points, interior = F)$trans_norm


# Plot polygons, points, and normals --------------------------------------
plot(sq)
plot(sq_pts, add = T)
plot(sq_norm_ls, add = T)
plot(tri, add = T)
plot(tri_points, add = T)
plot(norm_ls, add = T)

