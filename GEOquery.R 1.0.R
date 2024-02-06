# Test to run code to upload data from GEO
# The script is copied from the bioconductor manual:
# https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html

BiocManager::install("GEOquery", update = FALSE)
library (GEOquery)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Now, we are free to access any GEO accession. Note that in the following, I use a file packaged with 
# the GEOquery package. In general, you will use only the GEO accession, as noted in the code comments.

## If you have network access, the more typical way to do this would be to use this:
## gds <- getGEO("GDS507")
gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

# Now, gds contains the R data structure (of class GDS) that represents the GDS507 entry from GEO. You’ll note 
# that the filename used to store the download was output to the screen (but not saved anywhere) for later use to a call 
# to getGEO(filename=...).

# We can do the same with any other GEO accession, such as GSM11805, a GEO sample.


## If you have network access, the more typical way to do this would be to use this:
## gsm <- getGEO("GSM11805") 
gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))

# Each of these classes is comprised of a metadata header (taken nearly verbatim from the SOFT format header) 
# and a GEODataTable. The GEODataTable has two simple parts, a Columns part which describes the column headers on the Table part. 
# There is also a show method for each class. For example, using the gsm from above

# Look at gsm metadata:
head(Meta(gds))

# Look at data associated with the GSM:
# but restrict to only first 5 rows, for brevity
Table(gsm)[1:5,]

# Look at Column descriptions:
Columns(gsm)

