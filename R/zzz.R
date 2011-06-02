# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
".First.lib" <- function(lib, pkg) {
    library.dynam("fields", pkg, lib)
     packageStartupMessage(" Use help(fields) for an overview of this library\n
library( fields, keep.source=TRUE) retains comments in the source code.\n
Copyright 2004-2011, Licensed under GPL, www.gpl.org/licenses/gpl.html ")
}
