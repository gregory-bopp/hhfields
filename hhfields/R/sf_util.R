#' Rotate matrix counter clockwise through angle alpha
#'
#' @param angle of rotation
#'
#' @return rotation matrix
#' @export
rot <- function(angle){
  matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2,2)
}

#' Make square 
#' 
#' @description Make an sf square of side length 2*rad centered at (0,0)
#' @param rad half-length of square side
#'
#' @export
make_sq <- function(rad){
  sf::st_polygon(list(rbind(c(-rad, -rad),
                            c(rad, -rad),
                            c(rad,rad),
                            c(-rad, rad),
                            c(-rad, -rad))))
}


#' Sample points on boundary of polygon
#' @param n number of polygon points to sample
#'
#' @param poly polygon list
#' @param remove_vert should vertices be included in sample set
#'
#' @export
sample_poly_points <- function(n, poly, remove_vert = F){
  # Calculate perimeter
  perim <- poly %>% st_cast("MULTILINESTRING") %>% st_length
  nvert <- poly %>% st_coordinates() %>% nrow()
  if(remove_vert){
   # Calculate number of points to sample
    dfMax <- perim/(n - 2 + nvert)
  }
  else{
    dfMax <- perim/(n-1)
  }
poly_points <- st_segmentize(poly, dfMaxLength = dfMax) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::select(X, Y) %>%
  st_as_sf(coords = c("X", "Y"))
if(remove_vert){
  poly_points <- poly_points %>%
    st_difference(poly %>% st_cast("MULTIPOINT"))
}
return(poly_points)
}

#' Get linestrings from linestring or polygon
#' @param x linestring or polygon
#'
#' @export
get_linestrings <- function(x){
  in_geom <- st_geometry(x)

  if(!unique(st_geometry_type(in_geom)) %in% c("POLYGON", "LINESTRING")){
    stop("Input should be LINESTRING or POLYGON")
  }
  for(k in 1:length(st_geometry(in_geom))){
    sub <- in_geom[k]
    coords <- st_coordinates(sub)
    geom <- lapply(
      1:(nrow(coords) - 1),
      function(i)
        rbind(
          as.numeric(coords[i,   1:2]),
          as.numeric(coords[i+1, 1:2])
        ) %>% st_linestring
    ) %>%
      st_sfc()
  }

  if(k ==1){
    endgeom <- geom
  }
  else{
    endgeom <- rbind(endgeom, geom)
  }
  endgeom <- endgeom %>% st_sfc(crs = st_crs(x))
  return(endgeom)
}

#' Get normal vectors to a polygon
#' @param poly polygon
#'
#' @param poly_points points along polygon at which to find normal unit vectors
#'        to the polygon boundary
#' @param interior (logical; default=T) should the normal vectors point to the 
#' interior or the exterior of the polygon
#'
#' @export
get_poly_normals <- function(poly, poly_points, interior = F){
# Convert polygon into component linestrings
  poly_multi <- get_linestrings(poly)
  # Find linestring touching point
  line_idx <- st_is_within_distance(poly_points, poly_multi,
                                    sparse = T, dist = 1e-15) %>%
    unlist()
  # Get perpendicular lines rotated about points on poly
  perp_ls <- purrr::map2(st_geometry(poly_points), 1:length(line_idx),  ~{
    perp_l <- poly_multi[[line_idx[[.y]]]]*rot(pi/2)/1e14
    perp_l <- perp_l - st_centroid(perp_l) + .x  }) %>%
    st_sfc()

  normal_ls <- map2(perp_ls, st_geometry(poly_points), ~{
    X <- st_geometry(.x)
    X <- st_coordinates(X)[,1:2]
    p1 <- st_point(X[1,])
    p2 <- st_point(X[2,])
    if(!interior){
      id <- min(which(c(!st_contains(poly, p1, sparse = F),
                        !st_contains(poly, p2, sparse = F))))
    }
    else{
      id <- min(which(c(st_contains(poly, p1, sparse = F),
                        st_contains(poly, p2, sparse = F))))
    }
    out_pt <- list(p1,p2)[[id]]
    out_ls <- (st_linestring(rbind(st_coordinates(.y),
                                   st_coordinates(out_pt))) -
                 .y)*(2/st_length(.x))
    out_ls
  }) %>%
    st_sfc()

  # Translate linestrings to points where normals were calculated
  trans_norm <- map2(normal_ls, st_geometry(poly_points), ~{
  .x + .y
  }) %>%
    st_sfc()
  # Convert (un-translated() normal linestrings into points (representing
  # normal vectors)
  norm <- map(normal_ls, ~{
                    st_coordinates(.x)[2,1:2] %>%
                      st_point()}) %>%
                    st_sfc()
  return(list(norm = norm,
              trans_norm = trans_norm))
}










