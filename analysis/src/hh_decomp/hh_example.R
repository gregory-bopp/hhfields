#######################################################################
#
# Description: Perform a Helmholtz-Hodge decomposition of a vector field
#   with a bounded domain given the normal component of the divergence
#   free component on the boundary of the domain. For this example the
#   boundary of the domain is taken to be an outer square with an 
#   interior square removed. The target vector field is taken to be 
#   the gradient of a bi-variate normal density centered within the 
#   reconstruction domain.
#
#######################################################################
library(sf)
library(tidyverse)
library(here)
devtools::load_all(here('hhfields'))

set.seed(1)

# Parameters --------------------------------------------------------------
n_nodes <- 401
n_bndry <- 61

shift <- c(4,4)
# Make exterior square ----------------------------------------------------
sq <- make_sq(2) + st_point(shift)
sq_pts <- sample_poly_points(n_bndry, sq, remove_vert = T)
sq_norm <-  get_poly_normals(sq, sq_pts,  interior = T)$norm

# Make interior polygon ----------------------------------------------------
isq <- make_sq(0.4)*rot(pi/4) +st_point(c(0.5, 0.5))  +  st_point(shift)
isq_pts <- sample_poly_points(n_bndry, isq, remove_vert = T)
isq_norm <-  get_poly_normals(isq, isq_pts, interior = F)$norm

# Kernel knots/nodes ------------------------------------------------------
inner_nodes <- st_sample(sf::st_difference(sq, isq),
                         n_nodes, type = "hexagonal", exact = T)



# Cast to matrix ----------------------------------------------------------
boundary <- rbind(
                  # st_coordinates(sq_pts)[,1:2],
                  st_coordinates(isq_pts)[,1:2])

normals <- rbind(
                  # st_coordinates(sq_norm)[,1:2],
                 st_coordinates(isq_norm)[,1:2])
interior <- st_coordinates(inner_nodes)[,1:2]

all_nodes <- rbind(interior, boundary)

# True vector field -------------------------------------------------------
mu <-c(0,0) + shift
sigma <- diag(2)
dxdy <- -gauss_grad(all_nodes, mu = mu, sigma = sigma)
dx <- dxdy[,1]
dy <- dxdy[,2]


z <- mvtnorm::dmvnorm(all_nodes, mu, sigma)
df <- data.frame(x = all_nodes[,1],
                 y = all_nodes[,2],
                 dx = dx,
                 dy = dy,
                 z= z)

# HH Projection--------------------------------------------------------------
kg <- matern_kerngen(5)
hh <- hh_proj(dx, dy, Xi = all_nodes, Xb = boundary,
              Nb = normals, Xe = all_nodes, kerngen = kg)
df$vx <- hh$divfree[,1]
df$vy <- hh$divfree[,2]
df$wx <- hh$curlfree[,1]
df$wy <- hh$curlfree[,2]

print(paste0("Max Recon error:",max(abs(c(range(df$dx - (df$vx + df$wx)),
range(df$dy - (df$vy + df$wy)))))))


# Plot projections --------------------------------------------------------
arrow_head <-0.02
base_plot <- ggplot()+
  geom_sf(data = sq) +
  geom_sf(data = isq, fill = 4) +
  theme_bw()
p1 <- base_plot +
  geom_segment(data = df,
               aes(x = x, y = y, xend = x + dx, yend = y + dy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Target Field")
p2 <- base_plot +
  geom_segment(data = df,
               aes(x = x, y = y, xend = x + vx, yend = y + vy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Div. Free Projection")
p3 <- base_plot +
  geom_segment(data = df,
               aes(x = x, y = y, xend = x + wx, yend = y + wy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Curl Free Projection")
p4 <- base_plot +
  geom_segment(data = df,
               aes(x = x, y = y, xend = x + vx + wx, yend = y + vy + wy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Div. + Curl.")
cowplot::plot_grid(plotlist = list(p1, p4, p2, p3), nrow = 2, ncol = 2)
