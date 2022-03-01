#############################################################################
# Description: Example domain generation with normals calculated at a set of 
#    interior and boundary nodes. Domain consists of an outer square that has
#    had a set of interior polygons removed.
############################################################################

library(sf)
library(tidyverse)
library(here)
devtools::load_all(here('hhfields'))


set.seed(1)
# Outer square with side length = 8 
outer_poly <- make_sq(4)
# Maximum number of interior squares to generate
n_isq <- 6
# Simulate domain
dom <- sim_domain(outer_poly,
                  n_isq =  n_isq, 
                  n_bndry = 10,
                  n_inner =  101)
# Plot result
ggplot() +
  geom_sf(data = outer_poly) +
  geom_sf(data = dom$inner_poly) +
  geom_sf(data = dom$bndry_pts) +
  geom_sf(data = dom$trans_norm) +
  geom_sf(data = dom$inner_nodes, col = 2)
