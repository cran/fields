library( fields)
options(echo=FALSE)
test.for.zero.flag<- 1
#source("fields.diagonalize2.R")

data( ozone2)
x<- ozone2$lon.lat[1:50,]
knots<- cover.design( x, 10)$design

X<- cbind( fields.mkpoly(x,2), stationary.cov(x,knots, theta=4))
K<-  stationary.cov(knots,knots, theta=4)

M<- ncol( X)
M2<- nrow(K)
B<- matrix( 0, M,M)
B[ ((M-M2 +1): M),  ((M-M2 +1): M)]<- K
A<- t(X)%*%X

fields.diagonalize(A,B)-> look
fields.diagonalize2(A,B, verbose=FALSE)-> look2

test.for.zero( look$D, look2$D, tol=6E-8,tag="eigenvalues of both versions")

G1<- look$G
G2<-look2$G
a1<- sign( G1[1,])
a2<- sign(G2[1,])
a<- a1*a2


lambda<- .8
test.for.zero( solve( A + lambda* B), G2%*%diag( 1/(1+ lambda*look2$D))%*%t(G2), tag="inverse A+lambda*B" )
test.for.zero( solve( A + lambda* B), G1%*%diag( 1/(1+ lambda*look$D))%*%t(G1), tag="inverse A+lambda*B" )
test.for.zero( G2%*%diag( 1/(1+ lambda*look2$D))%*%t(G2) ,
                    G1%*%diag( 1/(1+ lambda*look$D))%*%t(G1), tag="inverse A+lambda*B" )


options( echo=TRUE)
cat("all done testing both versions of simultaneous diagonalization ", fill=TRUE)

