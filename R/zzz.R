.onAttach <- function(...) {
 
   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   cat("##\n## Model-based Ratings Figures\n")
   cat("## Copyright (C) 2007-", this.year,
      " Kevin M. Quinn and Daniel E. Ho \n", sep="")
   require(coda, quietly=TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("Ratings", libpath)
}

