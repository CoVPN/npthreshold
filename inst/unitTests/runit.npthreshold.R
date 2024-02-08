library("RUnit")
library("npthreshold")

test.random <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-6
#if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance=1e-6

set.seed(1)
#dat=rbigamma(n=500, shape.1=2, shape.2=2, rate.1=1, rate.2=2, rho=0.5) 
#checkEqualsNumeric(cor(dat[,1], dat[,2]), 0.4805187, tol=tolerance)





}
