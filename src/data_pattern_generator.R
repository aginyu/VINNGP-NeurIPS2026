
home.wd <- getwd()

source('./src/dual_mesh_loc.R')

Data_pattern_generator = function(beta_0, beta_1, nu, kappa,sigma2, random_seed = 1, max_edge = 0.5){
  set.seed(random_seed)
  
  
  # The factor log(2) is used to double the intensity function
  lg.s.c = rLGCP('exponential', im(beta_0 + log(2) + beta_1*gridcov, xcol=x0, yrow=y0),
                 var=sigma2, scale=kappa, win=win, saveLambda=TRUE)
  n = lg.s.c$n
  indices <- 1:n
  
  # Separate the test set. The test set contains a realization of Poisson process with the same intensity as the training set.
  keep <- rbinom(n, size = 1, prob = 0.5) == 1
  train_indices <- indices[keep]
  test_indices <- indices[(!keep)]
  xy.test <- cbind(lg.s.c$x[test_indices],lg.s.c$y[test_indices])
  
  
  lg.s.c$n <- length(train_indices)
  lg.s.c$x <- lg.s.c$x[train_indices]
  lg.s.c$y <- lg.s.c$y[train_indices]
  attr(lg.s.c, 'Lambda')$v <- attr(lg.s.c, 'Lambda')$v / 2
  
  

  # # Domain boundary
  # pl.dom = cbind(c(0,domain_max,domain_max,0), c(0,0,domain_max,domain_max))
  # bnd1 = Polygon(pl.dom, hole = FALSE)
  # bnd_hole = Polygon(pl.dom_hole, hole = TRUE)
  # bnd2 = Polygon(pl.dom_hole, hole = FALSE)
  # sp = SpatialPolygons(list(Polygons(list(bnd1, bnd_hole), "0"), 
  #                           Polygons(list(bnd2), "1")))
  # sp_hole = SpatialPolygons(list(Polygons(list(bnd_hole), "hole")))
  # bound = inla.sp2segment(sp)
  # xy.c = cbind(lg.s.c$x, lg.s.c$y)
  # mesh = inla.mesh.2d(xy.c, pl.dom, boundary = bnd1, max.edge = max_edge, cutoff = max_edge/5)
  # result_list = list(lg.s.c, xy.c, xy.test, mesh)
  # names(result_list) = c("lg_s_c", "xy_c_observed", "xy_c_test", "mesh")
  # return(result_list)
  
  # Domain boundary
  pl.dom = cbind(c(0, domain_max, domain_max, 0),
                 c(0, 0, domain_max, domain_max))
  
  bnd1 = Polygon(pl.dom, hole = FALSE)
  
  xy.c = cbind(lg.s.c$x, lg.s.c$y)
  
  mesh = inla.mesh.2d(
    xy.c, pl.dom,
    boundary = bnd1,
    max.edge = max_edge,
    cutoff = max_edge / 5
  )
  
  result_list = list(lg.s.c, xy.c, xy.test, mesh)
  names(result_list) = c("lg_s_c", "xy_c_observed", "xy_c_test", "mesh")
  return(result_list)
}

