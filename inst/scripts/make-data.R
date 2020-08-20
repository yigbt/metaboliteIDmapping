# This pipeline merges four different source of metabolite ID formats
# into one large data structure.
# Sources are: HMDB, Comptox, ChEBI, and graphite R package
#
# Final tibble contains more than 1.000.000 rows with 9 different
# ID formats.


# This pipeline utilizes several publicly accessible databases to
# retrieve different metabolite ID formats.
# In the following we will briefly describe the original sources:
#
#
# 1) Human Metabolome Database ( HMDB)
#
# Website: https://hmdb.ca/
# Current version: 4.0
# Download link: https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
# File format: XML
# ID formats: HMDB, CAS, Pubchem CID, KEGG, ChEBi, Drugbank
# Number of metabolites: 114100
#
#
#
# 2) Chemical Entities of Biological Interest (ChEBI)
#
# Website: https://www.ebi.ac.uk/chebi/
# Current version: Jan 5, 2020
# Download link: ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv
# File format: TSV
# ID formats: ChEBI, CAS, KEGG
# Number of metabolites: 17227
#
#
#
# 3) Comptox Chemical Dashboard
#
# Website: https://comptox.epa.gov/dashboard
# From comptox we retrieved two separated files, one linking to Pubchem and
# one linking to CAS numbers.
#
# I - Linking to Pubchem
# Download link: ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/PubChem_DTXSID_mapping_file.txt
# Current version: Nov 14, 2016
# File format: TSV
# ID formats: DTXSID, CID, SID
# Number of metabolites: 735553
#
# II - Linking to CAS registry numbers
# Download link: ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/2019/April/DSSTox_Identifiers_and_CASRN.xlsx
# Current version: Apr, 2019
# File format: XLSX
# ID formats: DTXCID, DTXSID, CAS
# Number of metabolites: 875755
#
# III - Full-join on both tables based on DTXSID
# ID formats: DTXCID, DTXSID, CAS, CID, SID
# Number of metabolites: 875796
#
#
#
# 4) Graphite R package
#
# Website: https://www.bioconductor.org/packages/release/bioc/html/graphite.html
# Current Version: Bioconductor release 3.11
# Data structure: date.frame
# Access from R package: graphite:::loadMetaboliteDb()@table
# ID formats: KEGG, ChEBI, CAS, Pubchem CID
# Number of metabolites: 155651




#####################################################################
#####################################################################
##
## A - Functions
##
#####################################################################
#####################################################################

#' Load a list of libraries
#'
#' @param packages A list of librarie to loaded
#' @param verbose Print the name of the library that is loaded. Default:FALSE
#'
load_libraries <- function(libraries, verbose=FALSE) {
    for(i in 1:length(libraries)) {
        if(verbose) {
            print(paste("Loading library :",libraries[i]," i:", i))
        }
        suppressPackageStartupMessages(library(libraries[i], character.only = TRUE))
        #library(libraries[i], character.only = T)
    }
}



#' Function for parsing huge xml files
#'
#' The XML file is retrieved from HMDB and is nearly 4GB big.
#' This function helps to extract solely the IDs of metabolites.
#'
ourBranches <- function() {

    store <- new.env()

    xmlGetValue <- function(x, node){
        val <- xpathSApply(x, node, xmlValue)
        ifelse( length( val) == 0, NA, val)
    }

    accession <- function(x, ...) {
        key <- xmlValue( getNodeSet(x, "//accession")[[1]])
        ids <- c( HMDB = key,
                  CAS = xmlGetValue(x, "//cas_registry_number"),
                  ChEBI = xmlGetValue(x, "//chebi_id"),
                  KEGG = xmlGetValue(x, "//kegg_id"),
                  CID = xmlGetValue(x, "//pubchem_compound_id"),
                  Drugbank = xmlGetValue(x, "//drugbank_id"),
                  Name = xmlGetValue(x, "//name"))
        store[[ key]] <- ids
    }

    getStore <- function() as.list(store)

    list(metabolite = accession, getStore=getStore)
}


#' Create tibble structure based on rows from HMDB DB
#'
#' Create an tibble structure for the intermediate joining of
#' comptox/graphite with HMDB database
#'
#' @param row Row entry from the HMBD DB, that could not be mapped
#'            Starting from these values, a tibble for the joined
#'            data frame will be created.
#'
#' @return tibble to be concatenated to the combined data structure
#'
create_tibble <- function( row = data.frame( HMDB = NA, CAS = NA, CID = NA,
                                             KEGG = NA, ChEBI = NA, Drugbank = NA)){

    df <- as_tibble( data.frame( CAS = ifelse( !is.na( row[['CAS']]), row[['CAS']], NA),
                                 DTXSID = NA, DTXCID = NA, SID = NA,
                                 CID = ifelse( !is.na( row[['CID']]), row[['CID']], NA),
                                 KEGG = ifelse( !is.na( row[['KEGG']]), row[['KEGG']], NA),
                                 ChEBI = ifelse( !is.na( row[['ChEBI']]), row[['ChEBI']], NA),
                                 HMDB = ifelse( !is.na( row[['HMDB']]), row[['HMDB']], NA),
                                 Drugbank = ifelse( !is.na( row[['Drugbank']]), row[['Drugbank']], NA))) %>%
        mutate_all( as.character)

    return( df)

}

#' Create tibble structure based on rows from ChEBI DB
#'
#' Create an tibble structure for the intermediate joining of
#' comptox/graphite/HMDB with ChEBI database
#'
#' @param row Row entry from the ChEBI DB, that could not be mapped
#'            Starting from these values, a tibble for the joined
#'            data frame will be created.
#'
#' @return tibble to be concatenated to the combined data structure
#'
create_tibble_chebi <- function( row = data.frame( CAS = NA, KEGG = NA, ChEBI = NA)){

    df <- as_tibble( data.frame( CAS = ifelse( !is.na( row[['CAS']]), row[['CAS']], NA),
                                 DTXSID = NA, DTXCID = NA, SID = NA, CID = NA,
                                 KEGG = ifelse( !is.na( row[['KEGG']]), row[['KEGG']], NA),
                                 ChEBI = ifelse( !is.na( row[['ChEBI']]), row[['ChEBI']], NA),
                                 HMDB = NA, Drugbank = NA)) %>%
        mutate_all( as.character)

    return( df)

}



#####################################################################
#####################################################################
##
## B - Prerequisites
##
#####################################################################
#####################################################################


Bioconductor_package_list <- c( "graphite")
R_package_list <-c("dplyr",  "magrittr", "XML", "readxl", "rappdirs", "tidyr", "stringr")

load_libraries( c( R_package_list, Bioconductor_package_list), verbose = TRUE)


recompute <- 0




#####################################################################
#####################################################################
##
## 1. Preprocessing
##
#####################################################################
#####################################################################


## The proprocessing step contais the downloading and re-formating
## of online metabolite databases
## The following online resources are used:
## 1) HMDB
## 2) ChEBI
## 3) Comptox
## 4) graphite R package


cat( "\nSection 1 - Preprocessing \n\n")


#####################################################################
## STEP 1 - HMDB IDs
#####################################################################

cat( "Step 1 - HMDB\n")

# 1) Download metabolites file from hmdb.ca
# 2) Unzip the file to get xml fileU
# 3) Extract ID-specific information of metabolites from xml file
# 4) create HMDB-specific tibble

hmdbfile <- file.path( rappdirs::user_cache_dir(),
                       "hmdb_mapping.rda",
                       fsep = .Platform$file.sep)

if( !file.exists( hmdbfile) || recompute){

    cat( "Download zipped HMDB database...\n")
    zipfile <- file.path( rappdirs::user_cache_dir(),
                          "hmdb.zip",
                          fsep = .Platform$file.sep)

    url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
    if( !file.exists( zipfile)) {
        download.file( url = url, destfile = zipfile,
                       method = "wget", quiet = TRUE)
    }

    cat( "Unzip Database...\n")
    unzippedfile <- file.path( rappdirs::user_cache_dir(),
                               "hmdb_metabolites.xml",
                               fsep = .Platform$file.sep)

    if( !file.exists( unzippedfile)){
        unzip( zipfile = zipfile, exdir = rappdirs::user_cache_dir())
    }

    cat( "Parse XML file...\n")
    branches <- ourBranches()
    xmlEventParse( unzippedfile, list(), branches = branches)

    cat( "Create tibble structure...\n")
    df <- lapply( branches$getStore(), function(x) x)

    cat( "Remove empty strings ...\n")
    hmdb <- as_tibble( do.call( rbind, df)) %>%
        dplyr::arrange( HMDB) %>%
        dplyr::mutate_if(is.character, list(~ dplyr::na_if(.,"")))

    cat( "Save extracted data...\n")
    save( hmdb, file = hmdbfile)

}else{

    load( file = hmdbfile)

}


## keep names of metabolites in separate tibble
hmdb_names <- hmdb %>% select( HMDB, Name)
hmdb <- hmdb %>% select( -Name)



#####################################################################
## Step 2 - ChEBI IDs
#####################################################################

cat( "Step 2 - ChEBI\n")


# 1) Download ID mapping table from ebi.ac.uk
# 2) Remove everything except for KEGG and CAS IDs
# 3) Create a tibble from overlapping IDs

chebifile <- file.path( rappdirs::user_cache_dir(),
                        "chebi_mapping.rda",
                        fsep = .Platform$file.sep)

if( !file.exists( chebifile) || recompute){

    destfile <- file.path( rappdirs::user_cache_dir(),
                           "chebi.tsv",
                           fsep = .Platform$file.sep)

    url <- "ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv"
    if( !file.exists( destfile)){

        cat( "Download ChEBI file ...\n")

        download.file( url = url,
                    destfile = destfile,
                    method = "wget",
                    quiet = TRUE)
    }

    cat( "Filter and rearrange ChEBI information ...\n")
    filecontent <- read.csv( file = destfile, sep = "\t", header= TRUE)

    chebi <- as_tibble( filecontent) %>%
        dplyr::filter( SOURCE == "KEGG COMPOUND") %>%
        dplyr::filter( TYPE == "KEGG COMPOUND accession" | TYPE == "CAS Registry Number" ) %>%
        dplyr::select( c(ID, COMPOUND_ID, TYPE, ACCESSION_NUMBER)) %>%
        dplyr::distinct( COMPOUND_ID, TYPE, ACCESSION_NUMBER, .keep_all = TRUE) %>%
        tidyr::spread( TYPE, ACCESSION_NUMBER) %>%
        dplyr::select( -ID) %>%
        dplyr::mutate_all( as.character)


    colnames( chebi) <- c( "ChEBI", "CAS", "KEGG")

    chebi <- chebi %>%
        dplyr::group_by( ChEBI) %>%
        dplyr::arrange( ChEBI) %>%
        tidyr::fill( KEGG, CAS, .direction = "downup") %>%
        dplyr::distinct() %>%
        ungroup()

    ## save rda file
    save( chebi, file = chebifile)

}else{

    load( file = chebifile)

}





#####################################################################
## Step 3 - Comptox IDs
#####################################################################

cat( "Step 3 - Comptox\n")

# 1) Download IDs from Comptox Dashboard in form of a html file
# 2) Get Pubchem-specific file
# 3) Get CAS-specific file
# 4) Merge them based on DTXSIDs

comptoxfile <-  file.path( rappdirs::user_cache_dir(),
                           "comptox_mapping.rda",
                           fsep = .Platform$file.sep)

if( !file.exists( comptoxfile) || recompute){

    masterfile <-  file.path( rappdirs::user_cache_dir(),
                              "comptox_master_download.html",
                              fsep = .Platform$file.sep)

    if( !file.exists( masterfile) || recompute){
        download.file( "https://comptox.epa.gov/dashboard/downloads",
                       destfile = masterfile,
                       method = "wget",
                       quiet = TRUE)
    }

    linenr <- grep( "DSSTox Identifier to PubChem Identifier Mapping File",
                    readLines( masterfile))
    url <- gsub( "^.*href=\"(.*txt)\".*", "\\1", readLines( masterfile)[ linenr])

    cat( "Download Pubchem IDs ...\n")
    destfile <- file.path( rappdirs::user_cache_dir(),
                           "comptox_pubchemd.tsv",
                           fsep = .Platform$file.sep)

    if( !file.exists( destfile)){
        download.file( url = url, destfile = destfile,
                       method =  "wget", quiet = TRUE)
    }

    comptox_pubchem <- read.csv( file = destfile, header = TRUE, sep = "\t")
    comptox_pubchem <- as_tibble( comptox_pubchem) %>% mutate_all( as.character)

    cat( "Download Comptox CAS IDs ...\n")
    linenr <- grep( "DSSTox identifiers mapped to CAS Numbers and Names File",
                    readLines( masterfile))
    url <- gsub( "^.*href=\"(.*xlsx)\".*", "\\1", readLines( masterfile)[ linenr])

    destfile <- file.path( rappdirs::user_cache_dir(),
                           "comptox_cas.xlsx",
                           fsep = .Platform$file.sep)

    if( !file.exists( destfile)){
        download.file( url = url, destfile = destfile,
                       method =  "wget", quiet = TRUE)
    }

    comptox_cas <- readxl::read_xlsx( path = destfile, sheet = 1) %>%
        dplyr::select( casrn, dsstox_substance_id, dsstox_structure_id, preferred_name)
    colnames( comptox_cas) <-  c( "CAS", "DTXSID", "DTXCID", "Name")

    ## combine both comptox data sets
    ## and remove those rows where DTXSIDs do not start with "DTXSID"
    comptox <- dplyr::full_join( comptox_cas, comptox_pubchem, by = "DTXSID") %>%
        filter( str_detect( DTXSID, "DTXSID"))

    ## save the tibble
    save( comptox, file = comptoxfile)

}else{

    load( file = comptoxfile)

}

## keep names of metabolites in separate tibble
comptox_names <- comptox %>% select( DTXSID, Name)
comptox <- comptox %>% select( -Name)




#####################################################################
## Step 4 - graphite IDs
#####################################################################

cat( "Step 4 - graphite\n")

# 1) use the in-house `graphite` IDS

graphite <- as_tibble( graphite:::loadMetaboliteDb()@table) %>%
    dplyr::select( PUBCHEM, KEGGCOMP, CHEBI, CAS) %>%
    dplyr::distinct()
colnames( graphite) <- c( "SID", "KEGG", "ChEBI", "CAS")

# subset the tibble to those which have at least a CAS OR SID
graphite_noNA <- graphite %>% filter( !is.na( CAS) | !is.na( SID))
graphite_NA <- graphite %>% filter( is.na( CAS) & is.na( SID))



#####################################################################
#####################################################################
##
## Merging
##
#####################################################################
#####################################################################

## Merge one database at a time using the following order:
## 1) Comptox and graphite
## 2) Add HMDB DB
## 3) Add ChEBI DB
## Since most of the steps are fairly time consuming, use cached
## versions if present.

cat( "\n\nSection 2 - Mergin of databases\n\n")



#####################################################################
## Step 1 - Merge comptox and graphite
#####################################################################

cat( "Step 1 - Comptox and graphite\n")

joined_comptox_graphite <-  file.path( rappdirs::user_cache_dir(),
                                       "joined_comptox_graphite.rda",
                                       fsep = .Platform$file.sep)

if( !file.exists( joined_comptox_graphite) || recompute){

    cat( "Join Comptox and graphite ...\n")
    j <- dplyr::full_join( comptox, graphite_noNA, by = c( "CAS", "SID"))

    cat( "Fill duplicated rows ... \n")
    nona <- j %>%
        dplyr::filter( !is.na( CAS)) %>%
        dplyr::group_by( CAS) %>%
        tidyr::fill( DTXCID, DTXSID, CID, KEGG, ChEBI, .direction = "downup") %>%
        dplyr::ungroup()
    isna <- j %>% dplyr::filter( is.na( CAS))


    comptox_graphite <- bind_rows( nona, isna)
    save( comptox_graphite, file = joined_comptox_graphite)

}else{

    load( file = joined_comptox_graphite)

}


#####################################################################
## Step 2 - Add HMDB database
#####################################################################


cat( "Step 2 - Add HMDB DB\n")

joined_comptox_graphite_hmdb <-  file.path( rappdirs::user_cache_dir(),
                                       "joined_comptox_graphite_hmdb.rda",
                                       fsep = .Platform$file.sep)

hmdb_unmappable <- hmdb %>%
    filter( is.na( CAS), is.na( ChEBI), is.na( KEGG), is.na( CID))
hmdb <- dplyr::setdiff( hmdb, hmdb_unmappable)

if( !file.exists( joined_comptox_graphite_hmdb) || recompute){


    comptox_graphite_hmdb <- comptox_graphite
    comptox_graphite_hmdb$HMDB <- as.character( NA)
    comptox_graphite_hmdb$Drugbank <- as.character( NA)

    noMatch <- create_tibble( )

    for( r in 1:nrow( hmdb)){

        if( r %% 1000 == 0){
            cat( "\t", r, " rows processed ... \n")
        }
        if( nrow( noMatch) %% 1000 == 0){
            cat( "\t", nrow( noMatch), " rows not matched ... \n")
        }

        row <- hmdb[ r, ]
        ids <- c()

        if( !is.na( row[['CAS']])){
            ids <- which( comptox_graphite_hmdb$CAS == row[['CAS']])
        }

        if( length( ids) == 0){

            if( !is.na( row[['ChEBI']])){
                ids <- which( comptox_graphite_hmdb$ChEBI == row[['ChEBI']])
            }
        }

        if( length( ids) == 0){

            if( !is.na( row[['KEGG']])){
                ids <- which( comptox_graphite_hmdb$KEGG == row[['KEGG']])
            }
        }

        if( length( ids) == 0){

            if( !is.na( row[['CID']])){
                cid_ids <- which( comptox_graphite_hmdb$CID == row[['CID']])
            }
        }


        if( length( ids) == 0){

            noMatch <- bind_rows( noMatch, create_tibble( row))

        }else{


            ## check if all retrieved rows from comptox_graphite have non-existing HMDB and DB IDs
            if( sum( is.na( comptox_graphite_hmdb[ ids, 8])) == length( ids) &&
                sum( is.na( comptox_graphite_hmdb[ ids, 9])) == length( ids) ){

                comptox_graphite_hmdb[ ids, 8] <- rep( row[['HMDB']], length( ids))
                if( !is.na( row[['Drugbank']])){
                    comptox_graphite_hmdb[ ids, 9] <- rep( row[['Drugbank']], length( ids))
                }

            }else{

                noMatch <- bind_rows( noMatch, create_tibble( row))

            }

        }

    }


    ## add the unmapped rows
    comptox_graphite_hmdb <- bind_rows( comptox_graphite_hmdb,
                                        noMatch)

    ## add those HMDB rows where only HMDB IDs were listed
    df <- apply( hmdb_unmappable, 1, create_tibble)
    comptox_graphite_hmdb <- bind_rows( comptox_graphite_hmdb,
                                        do.call( bind_rows, df))

    save( comptox_graphite_hmdb, file = joined_comptox_graphite_hmdb)

}else{

    load( file = joined_comptox_graphite_hmdb)

}

#####################################################################
## Step 3 - Add ChEBI database
#####################################################################

cat( "Step 3 - Add ChEBI DB\n")

joined_full <- file.path( rappdirs::user_cache_dir(),
                          "metabolitesMapping.rda",
                          fsep = .Platform$file.sep)

`%notin%` <- Negate( `%in%`)
metabolitesMapping <- comptox_graphite_hmdb

if( !file.exists( joined_full) || recompute){

    noMatch <- create_tibble_chebi( )

    ## Go through each row of the ChEBI data frame and check
    ## whether the ChEBI ID is already present in the mapping table.
    ## If yes, check if additional information on KEGG IDs or CAS numbers
    ## can be added by either filling NA gaps or adding complete new rows.
    ## Do the same in cases where the ChEBI IDs are not present.
    for( r in 1:nrow( chebi)){

        if( r %% 1000 == 0){
            cat("\t", r, " rows processed ... \n")
        }

        row <- chebi[r, c( 2,3,1)]
        ids <- which( metabolitesMapping$ChEBI == row[[ 'ChEBI']])

        ## ChEBI found
        if( length( ids) > 0){
            comp <- metabolitesMapping[ ids, c( 1,6,7)]
            check <- apply( comp, 1, function(x) all( x == row))

            if( TRUE %notin% check){

                if( !is.na( row[['KEGG']]) && is.na( row[[ 'CAS']])){


                    comp <- metabolitesMapping[ ids, c( 6,7)]
                    row_sub <- row[c(2,3)]
                    check <- apply( comp, 1, function(x) all( x == row_sub))

                    if( TRUE %notin% check){

                        boolean <- is.na( metabolitesMapping$KEGG[ ids])
                        if( TRUE %in% boolean){
                            ## in case NA is written in table
                            metabolitesMapping$KEGG[ ids][ boolean ] <- row[[ 'KEGG']]
                        }else{
                            noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                        }

                    }

                }

                if( !is.na( row[[ 'CAS']]) && is.na( row[[ 'KEGG']])){
                    comp <- metabolitesMapping[ ids, c( 1,7)]
                    row_sub <- row[ c(1,3)]
                    check  <- apply( comp, 1, function(x) all( x == row_sub))

                    if( TRUE %notin% check){

                        boolean <- is.na( metabolitesMapping$CAS[ ids])
                        if( TRUE %in% boolean){
                            ## in case NA is written in table
                            metabolitesMapping$CAS[ ids][ boolean] <- row[[ 'CAS']]
                        }else{
                            noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                        }

                    }
                }

                if( !is.na( row[[ 'CAS']]) && !is.na( row[[ 'KEGG']])){

                    boolean_cas <- is.na( metabolitesMapping$CAS[ ids])
                    boolean_kegg <- is.na( metabolitesMapping$KEGG[ ids])
                    if( TRUE %in% boolean_cas || TRUE %in% boolean_kegg){
                        ## in case NA is written in table
                        metabolitesMapping$CAS[ ids][ boolean_cas] <- row[[ 'CAS']]
                        metabolitesMapping$KEGG[ ids][ boolean_kegg] <- row[[ 'KEGG']]
                    }else{
                        noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                    }

                }

            }


        }else{

            ## ChEBI is not present in merged data set
            ## easiest case, both other IDs are NA -> add row
            if( is.na( row[[ 'KEGG']]) && is.na( row[[ 'CAS']])){
                noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
            }

            ## only CAS is present
            if( is.na( row[[ 'KEGG']]) && !is.na( row[[ 'CAS']])){
                ids <- which( metabolitesMapping$CAS == row[[ 'CAS']])
                if( length( ids) > 0){
                    boolean <- is.na( metabolitesMapping$ChEBI[ ids])
                    if( TRUE %in% boolean){
                        ## in case NA is written in table
                        metabolitesMapping$ChEBI[ ids][ boolean] <- row[[ 'ChEBI']]
                    }else{
                        noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                    }
                }else{
                    noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                }
            }

            ## only KEGG is present
            if( !is.na( row[[ 'KEGG']]) && is.na( row[[ 'CAS']])){

                ids <- which( metabolitesMapping$KEGG == row[[ 'KEGG']])
                if( length( ids) > 0){
                    boolean <- is.na( metabolitesMapping$ChEBI[ ids])
                    if( TRUE %in% boolean){
                        ## in case NA is written in table
                        metabolitesMapping$ChEBI[ ids][ boolean] <- row[[ 'ChEBI']]
                    }else{
                        noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                    }
                }else{
                    noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                }
            }

            ## both KEGG and CAS are present
            if( !is.na( row[[ 'KEGG']]) && !is.na( row[[ 'CAS']])){

                ids <- which( metabolitesMapping$KEGG == row[[ 'KEGG']])
                ids <- unique( ids, which( metabolitesMapping$CAS == row[[ 'CAS']]))

                if( length( ids) > 0){
                    boolean <- is.na( metabolitesMapping$ChEBI[ ids])
                    if( TRUE %in% boolean){
                        metabolitesMapping$ChEBI[ ids][ boolean] <- row[[ 'ChEBI']]
                    }else{
                        noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                    }
                }else{
                    noMatch <- bind_rows( noMatch, create_tibble_chebi( row))
                }

            }
        }
    }
    joined_full <- file.path( rappdirs::user_cache_dir(),
                              "metabolitesMapping.rda",
                              fsep = .Platform$file.sep)

    metabolitesMapping <- bind_rows( metabolitesMapping, noMatch)
    save( metabolitesMapping, file = joined_full, compress = "xz")

}else{

    load( file = joined_full)

}



#####################################################################
#####################################################################
##
## Adding common names
##
#####################################################################
#####################################################################


## From Comptox Dashboard and HMDB we were able to extract common names for
## the retrieved metabolites. These shall now be merged with the mapping table.
## 1) Join Comptox names based on their DTXSID
## 2) Fill the remaining HMDB gaps based on their HMDB identifier

## join names from Comptox
metabolitesMapping <- metabolitesMapping %>% full_join( comptox_names)

## joint remaining HMDB names
## These are 99560 entries.
## Get all HMDB entries without an associated Name
## Get their respective position in the metabolitesMapping table and the
## hmdb_names mapping table
## Insert a Position value in metabolitesMapping to make sure that duplicated
## HMDB entries both get a valid Name
## Remove the Position column in the end.
missing_values <- metabolitesMapping %>% filter( !is.na( HMDB) & is.na( Name)) %>% pull( HMDB)
missing_positions <- metabolitesMapping %>% filter( !is.na( HMDB) & is.na( Name)) %>% pull( Position)
hmdb_values_present <- metabolitesMapping %>% filter( !is.na( HMDB) & !is.na( Name)) %>% pull( HMDB)
tmp <- missing_values %notin% hmdb_values_present

missing_positions <- missing_positions[ tmp]
missing_values <- missing_values[ tmp]

hmdb_missing <- match( missing_values, hmdb_names$HMDB)
metabolitesMapping$Name[ missing_positions] <- hmdb_names$Name[ hmdb_missing]

metabolitesMapping <- metabolitesMapping %>% select( -Position)

save( metabolitesMapping, file = joined_full, compress = "xz")

