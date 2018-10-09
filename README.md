GEOSAMPLE - r package for generating a sample from a population of spatial locations.

TODO:
a) Add functionality to allow WGS84 longitude and latitude coordinates. Function currently uses centroids coordinates generated from a projected coordinates system.

b) Add functinality to allow direct extraction of buildings (or using its centroids) directly from OpenStreetMap.org


Stepwise example for current working function
##Step 1 Download OpenStreetMap.org buildings data (e.g. Sasa in Nigeria)
##Step 2 Transform the WGS84 coordinate system into a projected coordinate system (e.g. in the case of Sasa using the EPSG:26391 in QGIS software)
##Step 3 Generate the centroids of the buildings and convert file into a comma separated file (a .csv file)
##Step 4 Read the *.csv file or sample data "SasaSampleData.csv"

csvdata<-read.csv('./sampledata/SasaSampleData.csv', header=TRUE, stringsAsFactors=FALSE)

##Step 5 Set the seed number to avoid always generating different sets of random samples everything you execute the function
set.seed(16713)

##Step 6 Set parameters for the function and stored the return value: my.discrete.geosample <- discrete.geosample(population, xcolmun, ycolumn, minimumdistance, closepairs, samplesize, circleradius)

my.discrete.geosample <- discrete_geosample(csvdata, 385,1, 2, 18, 10,  5)

##Step 7 Extract the data frame data into newsample variable

newsample<-data.frame(my.discrete.geosample[1])

##Step 8 Plot original population data

par(pty="s",mfrow=c(1,2))

plot(csvdata[,1],csvdata[,2],pch=19,cex=0.25,xlab="X",ylab="Y",cex.lab=1,cex.axis=1,cex.main=1, main = "geosampled population in blue")

##Step 9 Plot new sampled population data

points(newsample[,1],newsample[,2],pch=19,col="blue")


Acknowledgement

[1] Chipeta, M.G., Terlouw, D.J., Phiri, K.S. and Diggle, P.J. (2017). Inhibitory
geostatistical designs for spatial prediction taking account of uncertain
covariance structure. Environmetrics, 28, DOI: 10.1002/env.2425

[2] Chipeta, M. G. (2016). Geostatistical design and analysis for estimating local variations in malaria disease burden Lancaster University 
