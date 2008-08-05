library(Ratings)

## get full Mondo Times data
data(Mondo)

## subsetting the Mondo data to include only raters who rated 5 or more 
## outlets
Mondo.sub <- Mondo[apply(!is.na(Mondo), 1, sum) >= 5, ]

## also getting rid of outlets that are not rated now
Mondo.sub <- Mondo.sub[,apply(is.na(Mondo.sub), 2, mean) != 1] 



## make sure a device is open
if(dev.cur() <= 1) get(getOption("device"))()
## if(dev.interactive())
## is FALSE in demo() {using source()}:
## Fix it later, copy its code for now
is.dev.interactive <- .Device %in% c("X11", "GTK", "gnome",
				"quartz", "windows", "JavaGD")
op <- par(ask = is.dev.interactive)




## some simple non-model-based plots of selected outlets
starplotObserved(Mondo.sub[,1:2], starcol="green",
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))

starplotObserved(Mondo.sub[,5,drop=FALSE], starcol="green",
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))


starplotObserved(Mondo.sub[,1:12], starcol="green",
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))

barplotObserved(Mondo.sub[,1:2], barcol="blue",
                top.limit=.8, scale.factor=.9,
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))

barplotObserved(Mondo.sub[,5,drop=FALSE], barcol="blue",
                top.limit=1, scale.factor=.7,
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))


barplotObserved(Mondo.sub[,1:12], barcol="blue",
                top.limit=1, scale.factor=.6,
                 catnames = c("Awful", "Poor", "Average", 
                   "Very Good", "Great"))






## fit an ordinal IRT model to the data
ord.out <- ordrating(Mondo.sub, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=5, burnin=5000, mcmc=20000, verbose=1000)

## get posterior predictive rating probabilities
tau <- tauCalculate(ord.out, 500)
colnames(tau) <- c("Awful", "Poor", "Average", 
                   "Very Good", "Great")

## just the labeled outlets
tau.sub <- tau[-grep("thetaOutlet", rownames(tau)),]

## clean up names
rownames(tau.sub) <-  gsub("theta", "", rownames(tau.sub))



## some model-based barplots
barplotModelBased(tau.sub, barcol="blue")

barplotModelBased(tau.sub[c(3,7,9),], barcol="blue")


## model-based starplot
starplotModelBased(tau.sub)

## different color scheme
mycol <- rgb(red=c(255, 243, 231, 219, 207, 159, 82.5, 30, 0, 0, 255),
             green=c(255.0, 250.1, 245.2, 240.3, 235.4, 215.8, 175, 144,
               0, 0, 69),
             blue=c(255, 255, 255, 255, 255, 255, 255, 255, 238, 183, 0),
             maxColorValue=255
             )

starplotModelBased(tau.sub, colvec=mycol)

