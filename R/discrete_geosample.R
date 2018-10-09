#' A function to geosample a population of spatial locations.
#'
#' @param population A data frame containing all potential sampling locations and covariates (if any). If there are no covariates, this will be a matrix of x-y coordinates (in projected coordinate system) for all potential sampling locations.
#' @param samplesize The required total sample size n.
#' @param xcolumn An integer indicating the column number in the population dataframe for the projected coordinates (X) in metres.
#' @param ycolumn An integer indicating the column number in the population dataframe for the projected coordinates (Y) in metres.
#' @param minimumdistance Inhibition distance or minimum distance between any two locations in the population sample.
#' @param closepairs Number of close pairs locations (must be between 0 and n/2).
#' @param circleradius Radius of a circle with centre x being one of the primary (pupulation - closepairs) points within which close pairs are placed.
#' @keywords geosampling
#' @export
#' @return The function returns a list of two items consisting of a data frame for sampled locations and their covariates (if any). Otherwise, this will be an n x 2 matrix of X - Y coordinates for sampled locations.
#' @examples
#' ##Step 1 Download OpenStreetMap.org buildings data (e.g. Sasa in Nigeria)
#' ##Step 2 Transform the WGS84 coordinate system into a projected coordinate system (e.g. in the case of Sasa using the EPSG:26391 in QGIS software)
#' ##Step 3 Generate the centroids of the buildings and convert file into a comma separated file (a .csv file)
#' ##Step 4 Read the *.csv file or sample data "SasaSampleData.csv"
#' csvdata<-read.csv('./sampledata/SasaSampleData.csv', header=TRUE, stringsAsFactors=FALSE)
#' ##Step 5 Set the seed number to avoid always generating different sets of random samples everything you execute the function
#' set.seed(16713)
#' ##Step 6 Set parameters for the function and stored the return value: my.discrete.geosample <- discrete.geosample(population, xcolmun, ycolumn, minimumdistance, closepairs, samplesize, circleradius)
#' my.discrete.geosample <- discrete_geosample(csvdata, 385,1, 2, 18, 10,  5)
#' ##Step 7 Extract the data frame data into newsample variable
#' newsample<-data.frame(my.discrete.geosample[1])
#' ##Step 8 Plot original population data
#' par(pty="s",mfrow=c(1,2))
#' plot(csvdata[,1],csvdata[,2],pch=19,cex=0.25,xlab="X",ylab="Y",cex.lab=1,cex.axis=1,cex.main=1, main = "geosampled population in blue")
#' ##Step 9 Plot new sampled population data
#' points(newsample[,1],newsample[,2],pch=19,col="blue")


discrete_geosample <- function (population, samplesize, xcolumn, ycolumn, minimumdistance, closepairs, circleradius) {

  # set coordinates columns
  xycolmuns<-xcolumn:ycolumn

  if (!is.matrix(population) & !is.data.frame(population))
    stop("object must be a matrix or data.frame.")
  if (any(is.na(population[, xycolmuns]))) {
    warning("NA's not allowed in the coordinates.")
    population <- population[complete.cases(population), drop = FALSE]
    warning("eliminating rows with NA's.")
  }
  if(any(closepairs>samplesize/2)){
    stop("Close pairs must be between 0 and samplesize/2.")
  }
  # Inhibition distance varying with closepairs
  minimumdistance <- minimumdistance * sqrt(samplesize/(samplesize - closepairs))
  dsq <- minimumdistance*minimumdistance
  dif <- samplesize-closepairs
  if(any(circleradius>minimumdistance/2)){
    circleradius = minimumdistance/2
    warning("circleradius > minimumdistance/2, circleradius=minimumdistance/2 will be used.")
  }
  # Random sample without replacement.
  xy.all <- population[, xycolmuns]
  N <- dim(xy.all)[1]
  index <- 1:N
  index.sample <- sample(index, dif, replace = FALSE)
  xy.sample <- xy.all[index.sample,]
  # Inhibition process for the "pupulaton - closepairs" design points.
  for (i in 2:dif){
    dmin <- 0
    while (dmin < dsq){
      take <- sample(index, 1)
      dvec <- (xy.all[take, 1] - xy.sample[, 1])^2 +
        (xy.all[take, 2] - xy.sample[,2])^2;dvec
      dmin <- min(dvec);dmin
    }
    xy.sample[i,] <- xy.all[take,]
  }
  colnames(xy.sample) <- c("x", "y")
  # Close pairs sampling.
  if (closepairs>0) {
    xy.cp <- matrix(NA, nrow = closepairs, ncol = 2)
    cp.mat<-matrix(sample(1:dif,closepairs,replace=FALSE),closepairs,2)
    for (j in 1:closepairs){
      take1<-cp.mat[j,1]; take2<-cp.mat[j,2]
      xy1<-c(xy.sample[take1,]); xy1 <- as.numeric(unlist(xy1))
      angle<-2*pi*runif(1, min = 0, max = 1)
      radius<-circleradius*sqrt(runif(1, min = 0, max = 1))
      if(any(radius<minimumdistance/4)){
        radius = minimumdistance/4
      }
      xy.cp[j,] <-xy1+radius*c(cos(angle),sin(angle))
    }
    colnames(xy.cp) <- c("x", "y")
    xy.sample <- rbind(xy.sample, xy.cp)
  }
  # Subset population for sampled locations.
  ind.coords <- NULL
  for(i in 1:nrow(xy.sample)) {
    ind.sel <- which(xy.sample[i,1]==
                       population[,xycolmuns[1]] &
                       xy.sample[i,2]==
                       population[,xycolmuns[2]])
    ind.coords <- c(ind.coords,ind.sel)
  }
  sampledpopulation <- population[ind.coords,]
  # Return results.
  return(list(sampledpopulation = sampledpopulation, minimumdistance = minimumdistance))
}
