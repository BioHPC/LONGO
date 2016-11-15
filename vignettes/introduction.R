## ---- echo=FALSE---------------------------------------------------------
   biomaRt::listDatasets(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))[1]

## ---- echo=FALSE---------------------------------------------------------
   biomaRt::listAttributes(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org",
                         dataset = "hsapiens_gene_ensembl"))[1]

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

