# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
".First.lib" <- function(lib, pkg) {
    library.dynam("fields", pkg, lib)
    packageStartupMessage(" Try help(fields) for an overview of this library\nfields web: http://www.image.ucar.edu/Software/Fields ")
}
