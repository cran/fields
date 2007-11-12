library( fields)
# tests of predict.se
# using approximations for conditional simulation on a grid. 
#
options( echo=FALSE)

data( ozone2)
as.image(ozone2$y[16,], x= ozone2$lon.lat, ncol=24, nrow=20, 
          na.rm=TRUE)-> dtemp
#
# A useful discretized version of ozone2 data
 
x<- cbind(dtemp$x[dtemp$ind[,1]], dtemp$y[dtemp$ind[,2]])
y<- dtemp$z[ dtemp$ind]
weights<- dtemp$weights[ dtemp$ind]

Krig( x, y, Covariance="Matern", 
   theta=1.0, smoothness=1.0, weights=weights) -> out



# the grid ...

glist<- list( x= dtemp$x, y=dtemp$y)

set.seed( 233)
sim.Krig.grid( out, grid= glist, M=200, extrap=TRUE)-> look

predict.surface.se( out, grid=glist, extrap=TRUE)-> test

look2<- matrix( NA, 20,24)

for(  k in 1:24){
for ( j in 1:20){
look2[j,k] <- sqrt(var( look$z[j,k,], na.rm=TRUE))
}
}


test.for.zero(  mean( abs(look2- test$z)/test$z), 0, relative=FALSE,
tol=.05, tag="Conditional simulation marginal se for grid")

#
# test for covariances


ind0<- expand.grid( c(1,4,5,10), c(3,4,5,10, 15))

x0<- cbind( glist$x[ind0[,1]], glist$y[ind0[,2]])
look2<- matrix( NA, 200,20)
for(  k in 1:20){
look2[,k] <- look$z[ ind0[k,1], ind0[k,2],]}

predict.se( out, x0, cov=TRUE)-> test2
ds<- 1/sqrt(diag(test2))
test3<-  diag(ds)%*% test2 %*% diag(ds)

#check plot( diag( test2), diag( var( look2)))

# Another plot to look at plot( c(test3), c(cor(look2)))

hold<-cor(look2)
upper<- col(hold)> row( hold)
dd<- (c(hold)- c(test3))[upper]

test.for.zero(   mean( abs(dd)) ,0, relative=FALSE,
tol=.05, tag="Conditional simulation correlations for grid (RMSE) ")

cat( "all done with grid based se tests", fill=TRUE)
options( echo=TRUE)