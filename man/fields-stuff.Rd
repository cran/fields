%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2024 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2
%##END HEADER
%##END HEADER

\name{fields-stuff}
\alias{fields.duplicated.matrix}
\alias{fields.mkpoly}
\alias{fields.derivative.poly}
\alias{fields.evlpoly}
\alias{fields.evlpoly2}

\title{Fields supporting functions}
\description{
Some supporting functions that are internal to fields top level
methods. Variants of these might be found in the R base but these
have been written for cleaner code or efficiency.
}
\usage{

fields.duplicated.matrix(mat, digits = 8) 

fields.mkpoly(x, m = 2, tag = "term")

fields.derivative.poly(x, m,dcoef)

fields.evlpoly( x, coef)

fields.evlpoly2( x, coef, ptab)
}

\arguments{
\item{mat}{Arbitrary matrix for examining rows}

\item{digits}{Number of significant digits to use for comparing
elements to determine duplciate values. }

\item{x}{Arbitrary matrix where rows are components of a
multidimensional vector}

\item{m}{ The null space degree  -- results in a polynomial of degree
(m-1) }

\item{dcoef}{ Coefficients of a multidimensional polynomial}

\item{coef}{Polynomial coefficients.}

\item{ptab}{Table of powers of different polnomial terms.}

\item{tag}{ \code{mkpoly} fills in as columns names the higher order
terms of the polynomial terms. The tag is the  text string prefix for
thes column names and the powers of the individual variables are
appended. Default is just "terms".}

}

\details{

\code{fields.duplicated} finds duplicate rows in a matrix.  The digits
arguments is the number of digits that are considered in the
comparison. 
The returned value is an array of integers from 1:M where M is the
number of unique rows and duplicate rows are referenced in the same
order that they appear as the rows of \code{mat}. 

\code{fields.mkpoly}  computes the complete matrix of all monomial
terms up to degree (m-1). Each row of \code{x} is are the componets of
a vector.  (The fields function mkpoly returns the number of these
terms.) In 2 dimensions with m=3 there   6 polynomial terms  up to
quadratic ( 3-1 =2) order and will be returned as the matrix:

cbind(  1 , x[,1], x[,2], x[,1]**2, x[,1]*x[,2], x[,2]**2  )

This function is used for the fixed effects polynomial  or spatial
drift used in spatial estimating functions Krig, Tps and mKrig. 
The matrix ptab is a table of the powers in each term for each
variable and is included as an attribute to the matrix returned by
this function. 
See the \code{attr} function for extracting an attribute from an
object. 

\code{ptab} for the example above is 
\preformatted{
    [,1] [,2]
[1,]    0    0
[2,]    1    0
[3,]    0    1
[4,]    2    0
[5,]    1    1
[6,]    0    2
}

This information is used in finding derivatives of the polynomial is also used
to create column names  for the terms that are of higher order than linear. 

\code{fields.deriviative.poly} finds the partial derivative matrix of
a multidimensional polynomial of degree (m-1) at different vector
values and with coefficients \code{dcoef}. 
This function has been orgainzed to be a clean utility for the
predicting the derivative of the estimated function from Krig or
mKrig. 
Within the fields context
the polynomial itself would be evaluated as 
     fields.mkpoly( x,m)\%*\%dcoef.
If x has  d columns ( also the dimension of the polynomial) and n rows
the partial derivatives of this polynomial at the locations x can be
organized in a nXd matrix.  This is the object returned by ths
function.   

\code{evlpoly} and \code{evlpoly2} are FORTRAN based functions for 
evaluating univariate polynomials and multivariate polynomials.  The
table of powers (ptab) needed for evlpoly2 is the same format as that
returned my the fields.mkpoly function. 

}
\author{Doug Nychka}
\keyword{spatial}

