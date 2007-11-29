library( fields)
options( echo=FALSE)

# test data
data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]

# turning spam on and off
Krig(x,y, cov.function = "stationary.taper.cov", theta=1.5,
      cov.args= list( spam.format=FALSE,
        Taper.args= list( theta=2.0,k=2) )
           ) -> out1

Krig(x,y, cov.function = "stationary.taper.cov", lambda=2.0, theta=1.5,
      cov.args= list( spam.format=TRUE,
        Taper.args= list( theta=2.0,k=2) )
           ) -> out2

temp1<- predict( out1,lambda=2.0)
temp2<- predict( out2)
test.for.zero( temp1, temp2)

#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]

# now look at mKrig w/o sparse matrix 
mKrig( x,y, cov.function="stationary.cov", theta=100, lambda=.3)-> look

Krig( x,y, cov.function="stationary.cov", theta=100, lambda=.3)-> look2

test.for.zero( look$d, look2$d)
test.for.zero( look$c, look2$c)


set.seed(123)
xnew<- cbind( (runif(20)-.5)*5, (runif(20)-.5)*5)
 temp<- predict( look, xnew)
 temp2<- predict( look2, xnew)
test.for.zero( temp, temp2)

# test of matrix of obs
N<- length( y)
Y<- cbind( runif(N), y,runif(N), y)

mKrig( x,Y, cov.function="stationary.cov", 
        theta=100, lambda=.3)-> lookY
temp3<-  predict( lookY, xnew)[,4]

test.for.zero( temp, temp3)
predict.surface( look)-> temp
predict.surface( look2)-> temp2

good<- !is.na( temp2$z)
test.for.zero( temp$z[good], temp2$z[good])

# testing stationary taper covariance with max.points option
# and also surface prediction

N<- length( y)
mKrig( x,y, cov.function="stationary.taper.cov", theta=2, 
spam.format=FALSE, lambda=.35, mean.neighbor=100 )-> look

Krig( x,y, cov.function="stationary.taper.cov", theta=2, 
spam.format=FALSE, lambda=.35)-> look2

predict.surface( look, nx=50, ny=45)-> temp
predict.surface( look2, nx=50, ny=45)-> temp2

good<- !is.na( temp2$z)
test.for.zero( temp$z[good], temp2$z[good])

# 
# Use Wendland with sparse off and on.
Krig( x,y, cov.function="wendland.cov", 
       cov.args=list( k=2, theta=2.8), 
       lambda=.3, spam.format=FALSE)-> look

mKrig( x,y, cov.function="wendland.cov",k=2, theta=2.8,
      spam.format=FALSE, lambda=.3)-> look2

mKrig( x,y, cov.function="wendland.cov",k=2, theta=2.8,
      spam.format=TRUE, lambda=.3)-> look3

# final tests is predict.
set.seed(223)
xnew<- cbind(runif( N)*.5 + x[,1], runif(N)*.5 + x[,2])
 temp<- predict( look, xnew)
 temp2<- predict( look2, xnew)
 temp3<- predict( look3, xnew)
test.for.zero( temp, temp2)
test.for.zero( temp2, temp3)

### bigger sample size
set.seed( 334)
N<- 1000
x<- matrix( runif(2*N),ncol=2)
y<- rnorm( N)

mKrig( x,y, cov.function="wendland.cov",k=2, theta=.1, lambda=.3,
max.points= 1e5)-> look2


test.for.zero( look2$non.zero.entires, 30146)

###### 
### test out passing to chol

data( ozone2)
     y<- ozone2$y[16,]
     good<- !is.na( y)
     y<-y[good]
     x<- ozone2$lon.lat[good,]

     # interpolate using defaults (Exponential)
     # stationary covariance
     mKrig( x,y, theta = 1.5, lambda=.2)-> out
     #
     # NOTE this should be identical to 
     Krig( x,y, theta=1.5, lambda=.2) -> out2
      temp<- predict( out)
      temp2<- predict( out2)

      test.for.zero( temp, temp2)

# test passing arguments for chol and max.points

set.seed( 334)
N<- 300
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 3*pi*x[,1])*sin( 3.5*pi*x[,2]) + rnorm( N)*.01


 Krig( x,y, Covariance="Wendland",
      cov.args= list(k=2, theta=.8, dimension=2),                   , 
       give.warnings=FALSE,
       lambda=1e2) -> out

 mKrig( x,y, 
            cov.function="wendland.cov",k=2, theta=.8, 
            lambda=1e2, 
            mean.neighbor=150, 
            chol.args=list( memory=list( tmpmax=1e5)), 
             )-> out2

 temp<- predict( out)
      temp2<- predict( out2)

      test.for.zero( temp, temp2)

cat("all done with mKrig tests", fill=TRUE)
options( echo=TRUE)

