.onLoad <- function(...) {
}

.onUnload <- function(libpath) {
    library.dynam.unload("Ratings", libpath)
}

