# global reference to umap (will be initialized in .onLoad)
#' @importFrom reticulate import_builtins
python_builtins <- NULL
.onLoad <- function(libname, pkgname) {
  try(
    {
      python_builtins <<- reticulate::import_builtins()
    },
    silent = TRUE
  )
}
