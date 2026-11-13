# rm(list = ls())
#setwd("C:/Users/lys/Desktop/supplement")
getwd()

library("INLA")
library("Matrix")
library("sp")
library("parallel")
library("lattice")

library(spatstat)
#library(RandomFields)
library(fields)
library(inlabru)


# data pattern matrix B
B = function(loc, mesh){
  B = inla.spde.make.A(mesh, loc)
  return(B)
}

C = function(mesh){
  fem = inla.mesh.fem(mesh, order=2)
  C = as(fem$c0, 'dgTMatrix')
  return(C)
}

G = function(mesh){
  fem = inla.mesh.fem(mesh, order = 2)
  G = as(fem$g1, 'dgTMatrix')
  return(G)
}

# Dual mesh
inla.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

# Dual mesh weight
DualMesh_parameter = function(mesh_spde){
  # Dual weight
  dual_weight = c()
  dual_mesh = inla.mesh.dual(mesh_spde)
  num_dualloc = length(dual_mesh)
  dual_loc = matrix(0, nrow = num_dualloc, ncol = 2)
  
  for (i in 1:num_dualloc) {
    dual_weight[i] = dual_mesh@polygons[[i]]@area
    dual_loc[i, ] = dual_mesh@polygons[[i]]@labpt
  }
  
  # dual mesh location
  A_dual = inla.spde.make.A(mesh_spde, dual_loc)
  
  # mesh location
  A_mesh = inla.spde.make.A(mesh_spde, mesh_spde$loc[ ,1:2])
  
  return(list(dual_weight, dual_loc, A_dual, A_mesh))
}
