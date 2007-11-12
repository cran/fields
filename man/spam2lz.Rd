\name{spam2lz}
\alias{spind2spam}
\alias{spam2spind}
\alias{spind2full}
\alias{spam2full}

\title{Conversion of formats for sparse matrices}
\description{
Some supporting functions that are internal to fields top level
methods. These are used to convert between the efficient but 
opaque format used by spam and more easily checked format based directly
on the row and column indices of non zero elements.
}
\usage{
spind2full(obj)

spam2full(obj)

spind2spam(obj)

spam2spind(obj)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
Either a list with the sparse index components (spind) or an 
obj of class spam.}
}


\details{
The differencee in formats is best illustarted by an example:

A 4X5 sparse matrix:
\preformatted{
     [,1] [,2] [,3] [,4] [,5]
[1,]    1    9    0    0   33
[2,]    0    0    0   26   34
[3,]    3   11    0   27   35
[4,]    0   12   20    0   36
}

spind format is a list with components "ind", "ra" and "da"
here is how the matrix above would be encoded:

\preformatted{
ind
      I  
 [1,] 1 1
 [2,] 1 2
 [3,] 1 5
 [4,] 2 4
 [5,] 2 5
 [6,] 3 1
 [7,] 3 2
 [8,] 3 4
 [9,] 3 5
[10,] 4 2
[11,] 4 3
[12,] 4 5

da
[1] 4 5

ra
 [1]  1  9 33 26 34  3 11 27 35 12 20 36
}

spam format is an S4 class with slot names 
"entries",  "colindices",  "rowpointers" and "dimension".

entries

 [1]  1  9 33 26 34  3 11 27 35 12 20 36

colindices

 [1] 1 2 5 4 5 1 2 4 5 2 3 5

rowpointers

 [1]  1  4  6 10 13

dimension

 [1] 4 5

The row pointers are the position in the array of entries where the next row
starts. 

NOTE: It is possible for the spind format to have a missing row of all
zeroes but this not allowed in spam format and produces an error message. 


}
\author{Doug Nychka}
\seealso{as.spam}
\keyword{spatial}
% at least one, from doc/KEYWORDS