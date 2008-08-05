

rgbinterp <- function(colvec, subn){

  oldn <- length(colvec)

  newred <- col2rgb(colvec[1])[1]
  newgreen <- col2rgb(colvec[1])[2]
  newblue <- col2rgb(colvec[1])[3]

  for (i in 2:oldn){
    newred <- c(newred, seq(from=col2rgb(colvec[i-1])[1],
                            to=col2rgb(colvec[i])[1],
                            length.out=subn))
    newgreen <- c(newgreen, seq(from=col2rgb(colvec[i-1])[2],
                            to=col2rgb(colvec[i])[2],
                            length.out=subn))
    newblue <- c(newblue, seq(from=col2rgb(colvec[i-1])[3],
                            to=col2rgb(colvec[i])[3],
                            length.out=subn))

  }

  return(rgb(red=newred, green=newgreen, blue=newblue, maxColorValue=255))
  
}



## barplotModelBased
##
## Input: 
##
## Output:
##
barplotModelBased <- function(tau.mat, barcol="darkgray", top.limit=0.5,
                              scale.factor=1.0, ...){

  ## some constants
  ncat <- ncol(tau.mat)
  npaper <- nrow(tau.mat)

  par(mar=c(1.5, 8, 3, 6) + 0.1, mgp=c(3,1,0))
  plot(1:ncat, 1:ncat, type="n",
       xlim=c(0, ncat+2), ylim=c(-0.5, npaper+0.8), axes=FALSE,
       xlab="", ylab="", cex.main=1.75, ...)
  
  mtext(side=1, line=-4, at=1:ncat, text=colnames(tau.mat),
        las=2, cex=1.25)
  mtext(side=2, line=0, at=1:npaper+.25, text=rownames(tau.mat),
        las=2, cex=1.25)
  for(i in 1:npaper){
    for(j in 1:ncat){
      rect(j-0.5, i,
           j+0.5, i+tau.mat[i,j]*scale.factor, col=barcol,
           border="white")
    }
    mtext(4,text="Prob.     ",line=-1,at=i+0.5, cex=1)
    segments(ncat+0.8, i,
             ncat+0.8, i+top.limit*scale.factor)
    segments(ncat+0.8, i,
             ncat+0.9, i)
    segments(ncat+0.8, i+(0.5*top.limit)*scale.factor,
             ncat+0.9, i+(0.5*top.limit)*scale.factor)
    segments(ncat+0.8, i+top.limit*scale.factor,
             ncat+0.9, i+top.limit*scale.factor)
    text(ncat+0.9, i, "0", pos=4, cex=.9)
    text(ncat+0.9, i+(0.5*top.limit)*scale.factor,
         as.character(round(0.5*top.limit, 2)), pos=4, cex=.9)
    text(ncat+0.9, i+top.limit*scale.factor,
         as.character(round(top.limit, 2)), pos=4, cex=.9)    
  }

  segments(0,1:npaper,5.8,1:npaper,
           col="lightgray")

  
}




## barplotObserved
##
## Input: Y
##
## Output:
##
barplotObserved <- function(Y, catnames, barcol="darkgray",  top.limit=0.5,
                            scale.factor=1.0, ...){

  ## convert Y to tau.mat matrix
  npaper <- ncol(Y)
  ncat <- length(catnames)
  tau.mat <- matrix(NA, npaper, ncat)
  rownames(tau.mat) <- colnames(Y)
  colnames(tau.mat) <- catnames

  nratings <- matrix(NA, npaper, 1)
  for (p in 1:npaper){
    nratings[p] <- sum(!is.na(Y[,p]))
  }
  
  for (p in 1:npaper){
    denom <- sum(!is.na(Y[,p]))
    for (c in 1:ncat){
      tau.mat[p,c] <- sum(Y[,p] == c, na.rm=TRUE) / denom      
    }
  }
  
  
  ## plot
  par(mar=c(1.5, 8, 3, 6) + 0.1, mgp=c(3,1,0))
  plot(1:ncat, 1:ncat, type="n",
       xlim=c(0, ncat+2), ylim=c(-0.5, npaper+0.8), axes=FALSE,
       xlab="", ylab="", cex.main=1.75, ...)
  
  mtext(side=1, line=-4, at=1:ncat, text=colnames(tau.mat),
        las=2, cex=1.25)
  mtext(side=2, line=0, at=1:npaper+.25, text=rownames(tau.mat),
        las=2, cex=1.25)
  for(i in 1:npaper){
    for(j in 1:ncat){
      rect(j-0.5, i,
           j+0.5, i+tau.mat[i,j]*scale.factor, col=barcol,
           border="white")
    }
    mtext(4, text="Frac.      ", line=-1, at=i+0.5, cex=1)

    if (nratings[i] > 1){
      ratingstring <- paste("(", nratings[i], " Ratings)", sep="")
    }
    else{
      ratingstring <- paste("(", nratings[i], " Rating)", sep="")
    }
    
    mtext(4, text=ratingstring, line=1, at=i+0.5, cex=1, las=2)


    segments(ncat+0.8, i,
             ncat+0.8, i+top.limit*scale.factor)
    segments(ncat+0.8, i,
             ncat+0.9, i)
    segments(ncat+0.8, i+(0.5*top.limit)*scale.factor,
             ncat+0.9, i+(0.5*top.limit)*scale.factor)
    segments(ncat+0.8, i+top.limit*scale.factor,
             ncat+0.9, i+top.limit*scale.factor)
    text(ncat+0.9, i, "0", pos=4, cex=.9)
    text(ncat+0.9, i+(0.5*top.limit)*scale.factor,
         as.character(round(0.5*top.limit, 2)), pos=4, cex=.9)
    text(ncat+0.9, i+top.limit*scale.factor,
         as.character(round(top.limit, 2)), pos=4, cex=.9)    

  }

  segments(0,1:npaper,5.8,1:npaper,
           col="lightgray")



}








## starplotObserved
##
## Input:
##
## Output:
##
starplotObserved <- function(Y, catnames, starcol="darkgray",
                             starsize=0.2, ...){

  npaper <- ncol(Y)
  ncat <- length(catnames)
  nstar <- matrix(NA, npaper, 1)
  nratings <- matrix(NA, npaper, 1)
  for (p in 1:npaper){
    nstar[p] <- mean(Y[,p], na.rm=TRUE)
    nratings[p] <- sum(!is.na(Y[,p]))
  }

  ## plot
  mymat <- matrix(c(.5, 1, .5, .5, 1, .5, .5, 1, .5, .5, 1, .5, .5, 1, .5),
                  1, 15, byrow=TRUE)

  
  par(mar=c(7, 10, 4, 0) + 0.1)
  plot(1:ncat, 1:ncat, type="n",
       xlim=c(0, ncat+4), ylim=c(0.5, npaper+1), axes=FALSE,
       xlab="", ylab="", cex.main=1.75, ...)
  mtext(side=2, line=0.5, at=1:npaper, 
        text=colnames(Y),
        las=2,
        cex=1.25)
  mtext(side=1, line=0.5, at=1:ncat, text=catnames,
        las=2, cex=1.25)
  segments(0.5,1:npaper, ncat+0.5, 1:npaper,
           col="lightgray")
  for(i in 1:npaper){
    for(j in 1:ncat){
      if(nstar[i]>=j){
        symbols(j, i, stars=mymat,
                bg=starcol,
                add=T, fg="lightgray", inches=starsize)
      } else {
        symbols(j, i, stars=mymat,
                bg="white",
                add=T, fg="lightgray", inches=starsize)
      }
    }
  }

  ratingstringvec <- rep("", npaper)
  for (p in 1:npaper){
    if (nratings[p] > 1){
      ratingstringvec[p] <- paste("(", nratings[p], " Ratings)", sep="")
    }
    else{
      ratingstringvec[p] <- paste("(", nratings[p], " Rating)", sep="")
    }
  }

  text(ncat+2, 1:npaper, ratingstringvec, cex=1.25)


}





## starplot
##
## Input:
##
## Output:
##
starplotModelBased <- function(tau.mat, colvec=NULL, starsize=0.2,
                               interpolation.level=200, ...){

  ## some constants
  ncat <- ncol(tau.mat)
  npaper <- nrow(tau.mat)
  

  if (is.null(colvec)){
    colvec <- rgb(red=c(255, 243, 231, 219, 207, 159, 82.5, 30, 0, 0, 0),
                  green=c(255.0, 250.1, 245.2, 240.3, 235.4, 215.8, 175, 144,
                    0, 0, 0),
                  blue=c(255, 255, 255, 255, 255, 255, 255, 255, 238, 183, 0),
                  maxColorValue=255
                  )
  }
  colvec <- rgbinterp(colvec, interpolation.level)
  
  ncolors <- length(colvec)

  cutpoints <- seq(from=0, to=0.5, length.out=ncolors+1)
  cutpoints[1] <- -1
  cutpoints[ncolors+1] <- 2
  
  ## plotting
  mymat <- matrix(c(.5, 1, .5, .5, 1, .5, .5, 1, .5, .5, 1, .5, .5, 1, .5),
                  1, 15, byrow=TRUE)

  layout(matrix(c(1,1,2,2), 2, 2, byrow=TRUE), heights=c(1, .2))
  
  ## plot
  par(mar=c(8, 13, 2, 12) + 0.1)
  plot(1:ncat, 1:ncat, type="n",
       xlim=c(0.5, ncat+0.5), ylim=c(0.5, npaper+1), axes=FALSE,
       xlab="", ylab="", cex.main=1.75, ...)

  mtext(side=2, line=0.5, at=1:npaper, #text=rownames(starprob2),
        text=rownames(tau.mat),
        las=2,
        cex=1.25)
  mtext(side=1, line=0.5, at=1:ncat, text=colnames(tau.mat),
        las=2, cex=1.25)

  segments(0.5, 1:npaper, ncat+.5, 1:npaper,
           col="lightgray")
  for(i in 1:npaper){
    for(j in 1:ncat){
      local.col <- colvec[findInterval(tau.mat[i,j], cutpoints)]
      symbols(j, i, stars=mymat, bg=local.col, fg="lightgray",
            inches=starsize, add=TRUE)
            
    }
  }

  ## legend
  par(mar=c(4.5, 12, 1.5, 13) + 0.1, mgp=c(1,0.2,0), yaxt="n")
  newmat <- matrix(seq(from=0, to=0.5, length.out=ncolors), ncolors, 1)
  xx <- seq(from=0, to=0.5, length.out=ncolors)
  image(x=xx, y=1, z=newmat, col=colvec, xlab="", 
        ylab="", breaks=cutpoints, xlim=c(0,0.5), axes=F)
 # axis(1, at=seq(0,0.5,0.1), labels=c(seq(0,0.4,0.1),"0.5+"), cex.axis=1.35,
 #      line=0)
  axis(1, at=seq(0,0.5,0.1), labels=rep("", 6), line=0)
  mtext(side=1, line=1, at=seq(0,0.5,0.1), text=c(seq(0,0.4,0.1),"0.5+"),
        cex=1.15) 
  mtext(side=3, line=1,"Probability", cex=1.15)

  
}





  
