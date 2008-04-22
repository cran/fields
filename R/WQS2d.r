# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"WQS2d" <-
function (x, transpose = FALSE) 
{
    if (!transpose) 
        t(WQS(t(WQS(x))))
    else {
        WQS.T(t(WQS.T(t(x))))
    }
}
