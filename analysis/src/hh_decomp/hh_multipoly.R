#######################################################################
#
# Description: Perform a Helmholtz-Hodge decomposition of a vector field
#   with a bounded domain given the normal component of the divergence
#   free component on the boundary of the domain. For this example the
#   boundary of the domain is taken to be an outer square with a set of
#   interior polygons removed. The target vector field is taken to be 
#   the gradient of a bi-variate normal density centered within the 
#   reconstruction domain.
#
#######################################################################
library(sf)
library(tidyverse)
library(here)
devtools::load_all(here('hhfields'))

set.seed(1)
shift <- c(0,0)
outer_poly <- make_sq(4)
n_isq <- 30

# Simulate Domain ---------------------------------------------------------
dom <- sim_domain(outer_poly, n_isq, n_bndry = 10, 201)

ggplot() +
  geom_sf(data = outer_poly) +
  geom_sf(data = dom$inner_poly)+
  geom_sf(data = dom$bndry_nodes) +
  geom_sf(data = dom$trans_norm) +
  geom_sf(data = dom$inner_nodes, col = 2)


# Convert to coordinates --------------------------------------------------
coords  <- map(dom, ~rbind(st_coordinates(.x)[,1:2]))
coords$all_nodes <- rbind(coords$inner_nodes, coords$bndry_nodes)



# True vector field -------------------------------------------------------
mu <-c(-4,2) + shift
sigma <- diag(2)
dxdy <- -gauss_grad(coords$all_nodes, mu = mu, sigma = sigma)
dx <- dxdy[,1]
dy <- dxdy[,2]


z <- mvtnorm::dmvnorm(coords$all_nodes, mu, sigma)
df <- data.frame(x = coords$all_nodes[,1],
                 y = coords$all_nodes[,2],
                 dx = dx,
                 dy = dy,
                 z= z)

# HH Projection--------------------------------------------------------------
kg <- matern_kerngen(5)
hh <- hh_proj(dx, dy, Xi = coords$all_nodes, Xb = coords$bndry_nodes,
              Nb = coords$norm, Xe = coords$all_nodes, kerngen = kg)
df$vx <- hh$divfree[,1]
df$vy <- hh$divfree[,2]
df$wx <- hh$curlfree[,1]
df$wy <- hh$curlfree[,2]

print(paste0("Max Recon error:",max(abs(c(range(df$dx - (df$vx + df$wx)),
                                          range(df$dy - (df$vy + df$wy)))))))


# Plot projections --------------------------------------------------------
arrow_head <-0.02
base_plot <- ggplot()+
  geom_sf(data = outer_poly) +
  geom_sf(data = dom$inner_poly, fill = 4) +
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
plt <- cowplot::plot_grid(plotlist = list(p1, p4, p2, p3), nrow = 2, ncol = 2)

fig_dir <- here('analysis', 'fig')
cowplot::save_plot(filename = file.path(fig_dir, 'hh_multipoly.png'), 
       plot = plt, 
       base_height = 9)
