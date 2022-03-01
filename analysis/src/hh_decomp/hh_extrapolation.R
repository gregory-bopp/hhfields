library(sf)
library(tidyverse)
library(here)
devtools::load_all(here('hhfields'))

# Functions ---------------------------------------------------------------
build_target <- function(coords, mu, sigma){
  dxdy <- -gauss_grad(coords, mu = mu, sigma = sigma)
  dx <- dxdy[,1]
  dy <- dxdy[,2]
  dx <- dnorm(coords[,1], mu[1])
  dx <- dnorm(coords[,2], mu[2])

  z <- mvtnorm::dmvnorm(coords, mu, sigma)
  data.frame(x = coords[,1],
             y = coords[,2],
             dx = dx,
             dy = dy,
             z= z)
}


# Parameters --------------------------------------------------------------
n_nodes <- 401
n_bndry <- 8


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
  st_coordinates(sq_pts)[,1:2],
  st_coordinates(isq_pts)[,1:2])

normals <- rbind(
  st_coordinates(sq_norm)[,1:2],
  st_coordinates(isq_norm)[,1:2])
interior <- st_coordinates(inner_nodes)[,1:2]


train_nodes <- boundary
test_nodes <- interior

all_nodes <-  rbind(train_nodes,
                    test_nodes)


# True vector field -------------------------------------------------------
shift <- c(0,0)
mu <-c(0,0) + shift
sigma <- diag(2)
train <- build_target(train_nodes, mu, sigma)
test <- build_target(test_nodes, mu, sigma)


# HH Projection--------------------------------------------------------------
# Train
hh <- hh_proj(train$dx, train$dy, Xi = train_nodes, Xb = boundary,
              Nb = normals, Xe = train_nodes, kerngen = matern_kerngen(0.2))
train$vx <- hh$divfree[,1]
train$vy <- hh$divfree[,2]
train$wx <- hh$curlfree[,1]
train$wy <- hh$curlfree[,2]
# Test
hh <- hh_proj(train$dx, train$dy, Xi = train_nodes, Xb = boundary,
              Nb = normals, Xe = test_nodes, kerngen = matern_kerngen(0.2))
test$vx <- hh$divfree[,1]
test$vy <- hh$divfree[,2]
test$wx <- hh$curlfree[,1]
test$wy <- hh$curlfree[,2]

print(paste0("Max Recon error:",max(abs(c(range(test$dx - (test$vx + test$wx)),
                                          range(test$dy - (test$vy + test$wy)))))))

# Plot projections --------------------------------------------------------
arrow_head <-0.02
base_plot <- ggplot()+
  geom_sf(data = sq) +
  geom_sf(data = isq, fill = 4) +
  theme_bw()
p1 <- base_plot +
  geom_segment(data = train,
               aes(x = x, y = y, xend = x + dx, yend = y + dy),
               arrow = arrow(length = unit(arrow_head, "npc")),
               color = 'red') +
  geom_segment(data = test,
               aes(x = x, y = y, xend = x + dx, yend = y + dy),
               arrow = arrow(length = unit(arrow_head, "npc")))+
  ggtitle("Target Field")
p2 <- base_plot +
  geom_segment(data = train,
               aes(x = x, y = y, xend = x + vx, yend = y + vy),
               arrow = arrow(length = unit(arrow_head, "npc")),
               color = 'red') +
  geom_segment(data = test,
               aes(x = x, y = y, xend = x + vx, yend = y + vy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Div. Free Projection")
p3 <- base_plot +
  geom_segment(data = train,
               aes(x = x, y = y, xend = x + wx, yend = y + wy),
               arrow = arrow(length = unit(arrow_head, "npc")),
               color = 'red') +
  geom_segment(data = test,
               aes(x = x, y = y, xend = x + wx, yend = y + wy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Curl Free Projection")
p4 <- base_plot +
  geom_segment(data = train,
               aes(x = x, y = y, xend = x + vx+wx, yend = y + vy + wy),
               arrow = arrow(length = unit(arrow_head, "npc")),
               color = 'red') +
  geom_segment(data = test,
               aes(x = x, y = y, xend = x + vx + wx, yend = y + vy + wy),
               arrow = arrow(length = unit(arrow_head, "npc"))) +
  ggtitle("Div. + Curl.")
cowplot::plot_grid(plotlist = list(p1, p4, p2, p3), nrow = 2, ncol = 2)
