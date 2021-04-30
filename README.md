# The `metabolitIDmapping` R package

The `metaboliteIDmapping` AnnotationHub package provides a
comprehensive ID mapping for various metabolite ID formats. Within
this annotation package, nine different ID formats and metabolite
common names are merged in one large mapping table. ID formats include
[Comptox Chemical Dashboard](https://comptox.epa.gov/dashboard) IDs
(DTXCID, DTXSID), [Pubchem](https://pubchem.ncbi.nlm.nih.gov/) IDs
(CID, SID), [CAS Registry
numbers](https://www.cas.org/support/documentation/references)
(CAS-RN), [Human Metabolome Database](https://hmdb.ca/) (HMDB),
[Chemical Entities of Biological
Interest](https://www.ebi.ac.uk/chebi/) (ChEBI), [KEGG
Compounds](https://www.genome.jp/kegg/compound/) (KEGG), and
[Drugbank](https://www.drugbank.ca/) (Drugbank)

The metabolite IDs and names were retrieved from four different
publicly available sources and merged into one mapping table by means
of the R script that is distributed alongside the AnnotationHub
package. 

For detailed information about the data sources please have a look in
the vignette at our [Bioconductor page](https://bioconductor.org/packages/devel/data/annotation/vignettes/metaboliteIDmapping/inst/doc/metaboliteIDmapping.html)

# Installation

It is recommended to install the `metaboliteIDmapping` package via Bioconductor.
Therefore, start `R` (version 4.0) and enter:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("metaboliteIDmapping")
```

# Usage
There are two different ways to load the mapping ID table from this package.

First, simply load the `metaboliteIDmapping` package into your R session.
When the package is loaded, the data will be available as tibble:

```R
library( metaboliteIDmapping)

metabolitesMapping
```

Second, search for the mapping table in the AnnotationHub resource interface:

```R
library( AnnotationHub)

ah <- AnnotationHub()
datasets <- query( ah, "metaboliteIDmapping")
```

Currently, there are three versions of the mapping table. 

* AH79817 represents the original ID mapping containing 9 different ID formats
* AH83115 mapping table which also includes common names for each compound
* AH91792 current version of the table that accounts for tautomers 

For implanting this data in your code, it is recommended to use the
AHid for retrieval:
```R
data <- ah[["AH91792"]]
```

# LICENSE

Copyright (C) 2011 - 2020 Helmholtz Centre for Environmental Research
UFZ.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the UFZ
License document for more details:
<https://github.com/yigbt/metaboliteIDmapping/blob/master/LICENSE.md>
