# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"qr.yq2" <-
function (qr, y) 
{
    t(qr.q2ty(qr, t(y)))
}
