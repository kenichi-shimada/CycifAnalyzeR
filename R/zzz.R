.onLoad <- function(libname, pkgname) {
  # use super-assignment to update global reference to scipy
  if(reticulate::py_module_available("scimap")){
    sm <<- reticulate::import("scimap", delay_load = TRUE)
  }else{
    stop("scimap is not avaialble in the current environemnt: ",reticulate::py_config()$python)
  }
}
