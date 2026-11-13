logintensity_vis = function(mesh, field_value_list, case = c("hole", "full"), method = c("VINNGP", "INLA", "VIFRK"), contour_step,ref_values){
  values = unlist(ref_values)
  # log method
  # values = log(values + 1 - min(values))
  val_min = min(values)
  val_max = max(values)
  # contour_levels = seq(val_min, val_max, by = contour_step)
  field_value_list[field_value_list < val_min] = val_min
  field_value_list[field_value_list > val_max] = val_max
  
  # quantile method
  probs = seq(0, 1, length.out = 1/contour_step)
  contour_levels = quantile(values, probs)
  
  contour_levels = round(contour_levels, digits = 3)
  proj = inla.mesh.projector(mesh, dims=c(100,100))
  
  if(case == "hole"){
    if(method == "VINNGP"){
      field = field_value_list$hole$VINNGP
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }
    if(method == "INLA"){
      field = field_value_list$hole$INLA
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }
    if(method == "VIFRK"){
      field = field_value_list$hole$VIFRK
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }
  }
  if(case == "full"){
    if(method == "VINNGP"){
      field = field_value_list$full$VINNGP
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }
    if(method == "INLA"){
      field = field_value_list$full$INLA
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }
    if(method == "VIFRK"){
      field = field_value_list$full$VIFRK
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = viridis(256), zlim = c(val_min,val_max), xlab="", ylab="")
      contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }    
    if(method == "TRUE"){
      field = field_value_list
      image(proj$x, proj$y, inla.mesh.project(proj, field), col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(256)), zlim = c(val_min,val_max), xlab="", ylab="")
      #contour(proj$x, proj$y, matrix(inla.mesh.project(proj, field), 100), levels = contour_levels, add = TRUE, col = "black")
    }  
  }
}
