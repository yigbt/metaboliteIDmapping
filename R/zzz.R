##
## Load  metabolitesMapping into namespace
##

.onLoad <- function(libname, pkgname) {
    
    ns <- asNamespace( pkgname)
    ah <- AnnotationHub::AnnotationHub()
    metabolitesMapping <- ah[["AH91792"]]
        
    assign( "metabolitesMapping", metabolitesMapping, envir = ns)
    
    namespaceExport( ns, "metabolitesMapping")
    
}