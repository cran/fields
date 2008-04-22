# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

".First.lib" <-
function (lib, pkg) 
{library.dynam("fields",pkg, lib)
cat(" Try help(fields) for an overview of this library ", 
           fill=TRUE)
#"  Some name changes have been made to several common functions:",
#"  exp.-> Exp. and rad. -> Rad. See help( Exp.cov)", fill=TRUE)
}
