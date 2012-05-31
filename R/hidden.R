########## hidden functions to help in model implementation ##########


# return name of the calling function
"calling.function" <-
   function(parentheses=TRUE) {
     calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
     if (parentheses){
       calling.function <- paste(calling.function, "()", sep="")
     }
     return(calling.function)
   }


# parse the passed seeds
# 1] if a scalar is passed, it is used by Mersennse twister
# 2] if a list of length two is passed, a parallel-friendly stream is
#    created using L'Ecuyer
"form.seeds" <-
   function(seed) {
      if(length(seed)==1) {
         if(is.na(seed)) seed <- 12345
         seed <- as.integer(seed)
         if(seed < 0) {
            cat("Error: Mersenne seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)                       
         }
         seeds <- list(0, rep(seed,6), 0)
      }
      if(length(seed)==2) {
         if(!is.list(seed)) {
            cat("Error: List must be passed to use L'Ecuyer.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         lec.seed <- seed[[1]]
         lec.substream <- as.integer(seed[[2]])
         if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
         if(length(lec.seed) != 6) {
            cat("Error: L'Ecuyer seed not of length six.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if(!all(lec.seed >= 0))  {
             cat("Error: At least one L'Ecuyer seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if( max(lec.seed[1:3]) >= 4294967087){
           cat("Error: At least one of first three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294967087\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( all(lec.seed[1:3] == 0 )){
           cat("Error: first three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( max(lec.seed[4:6]) >= 4294944443){
           cat("Error: At least one of last three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294944443\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }         
         if( all(lec.seed[4:6] == 0 )){
           cat("Error: last three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if(lec.substream < 1) {
            cat("Error: L'Ecuyer substream number not positive.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)               
         }
         seeds <- list(1, lec.seed, lec.substream) 
      }
      if(length(seed)>2) {
            cat("Error: Seed passed as length greater than two.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)        
      }
      return(seeds)
   }

