##### Read in libraries and data ---------------

library(adegenet)
library(spdep)
library(adehabitat)
library(vcftools)

#Import and convert vcf into genind object asdd raster and spatial data 

# Read the VCF file
vcf <- read.vcfR("path_to_your_file.vcf")

# Convert to genind object
genind_obj <- vcfR2genind(vcf)

# Check the genind object
print(genind_obj)

# Pull in lat and long coordinates ensuring the order of individuals in spatial_data matches the order in genind_obj
genind_obj@other$xy <- spatial_data[, c("X", "Y")]

# Check the genind object with spatial data
print(genind_obj@other$xy)

# Load your raster data
raster_data <- raster("path_to_your_raster.tif")

# Assuming you have a genind object named genind_obj
# And a raster object named raster_data

# Extract raster values at the locations of the individuals
coordinates <- genind_obj@other$xy
raster_values <- extract(raster_data, coordinates)

# Attach the raster values to the genind object
genind_obj@other$raster_values <- raster_values

# Check the genind object with raster data
print(genind_obj@other$raster_values)

##### Create Delauney network and Triangulation matrix ---------

myCn <- chooseCN(genind_obj$other$xy, type=6, k=10, plot=FALSE)
myCn

mySpca <- spca(genind_obj,type=5,d1=0,d2=2,scannf=FALSE)

#Plot Eigen Values 

barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

#Graphical display of model results and Monte Carlo test of global vs local scores 

screeplot(mySpca)
myGtest <- global.rtest(obj$tab,mySpca$lw,nperm=99)
myGtest


myLtest <- local.rtest(obj$tab,mySpca$lw,nperm=99)
myLtest

plot(myLtest)

plot(mySpca)

#Plot first global score dependent oon Principal wieghtings 

colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(obj)$xy[,1]
y <- other(obj)$xy[,2]
temp <- interp(x, y, mySpca$li[,1])
image(temp, col=azur(100))
points(x,y)

interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)

myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
  title("sPCA - interpolated map of individual scores")
  points(x,y)
}
filled.contour(temp, color.pal=myPal, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())

# Establish contribution of alleles to clustering 

temp <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    fac=obj$loc.fac, cex.fac=0.6)
#Mapping and clustering of scores --------------

s.value(cbind(1:11,rep(1,11)), -5:5, cleg=0)
text(1:11,rep(1,11), -5:5, col="red",cex=1.5)

#Establish significance with Morans I

pc1.mctest <- moran.mc(rupica.pca1$li[,1], rupica.graph, 999)
plot(pc1.mctest)
moran.plot(rupica.pca1$li[,1], rupica.graph)
#Second score 
pc2.mctest <- moran.mc(rupica.pca1$li[,2], rupica.graph, 999)
plot(pc2.mctest)

#####Multi Variate spatial score check ---------------

#Initial mySPCA
Korke.spca1 <- spca(genind_obj, cn=rupica.graph,scannf=FALSE,
                     nfposi=2,nfnega=0)
barplot(Korke.spca1$eig, col=rep(c("red","grey"), c(2,1000)),
        main="rupica dataset - sPCA eigenvalues")


screeplot(Korke.spca1)

#PLot both first and second eigenvalue scores to establish genetic variance and clustering spatial autocorrelation 

showKorke()
s.value(genind_obj$other$xy, rupica.spca1$ls[,1], add.p=TRUE, csize=0.7)
title("sPCA - first PC",col.main="yellow" ,line=-2, cex.main=2)

showKorke()
s.value(genind_obj$other$xy, rupica.spca1$ls[,2], add.p=TRUE, csize=0.7)
title("sPCA - second PC",col.main="yellow" ,line=-2, cex.main=2)

showKorke()
colorplot(genind_obj$other$xy, rupica.spca1$ls, axes=1:2, transp=TRUE, add=TRUE,
          cex=3)
title("sPCA - colorplot of PC 1 and 2\n(lagged scores)", col.main="yellow",
      line=-2, cex=2)


