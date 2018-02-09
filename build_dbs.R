
# Build database for GO annotations
# using the bioconductor annotation packages


if (0) {
  library(RMariaDB)
  library(DBI)
  con2 <- dbConnect(RMySQL::MySQL(), username="root", password=rstudioapi::askForPassword(), dbname="annotate")
}


if (0) {
  # central annotation tables -- TABLE 'anno_[species]'
  library(org.Ce.eg.db)
  library(org.Dm.eg.db)
  
  drome_list = as.list(org.Dm.egGO2ALLEGS)
  cel_list = as.list(org.Ce.egGO2ALLEGS)
  
  drome = data.frame("goid"= gsub("\\..*","",names( unlist(drome_list) )),
                     "entrez"=unlist(drome_list),
                     "evidence"=gsub(".+\\.","",names( unlist(drome_list) )),
                     stringsAsFactors = F )
  
  cel = data.frame("goid"= gsub("\\..*","",names( unlist(cel_list) )),
                   "entrez"=unlist(cel_list),
                   "evidence"=gsub(".+\\.","",names( unlist(cel_list) )),
                   stringsAsFactors = F )
  
  dbWriteTable(con2, "anno_dro", drome,
               field.types=c(goid="CHAR(10)", 
                             entrez="INT UNSIGNED", 
                             evidence="CHAR(3)"), 
               row.names=F, overwrite=T)
  
  dbWriteTable(con2, "anno_cel", cel,
               field.types=c(goid="CHAR(10)", 
                             entrez="INT UNSIGNED", 
                             evidence="CHAR(3)"), 
               row.names=F, overwrite=T)
  
  
  # --> GO meta-data -- TABLE 'gon_[species]'
  drome_gon = data.frame("goid"=names(drome_list),
                         "ontology"=Ontology( names(drome_list) ),
                         "n"=sapply(drome_list, function(x) length(unique(x))), row.names = NULL )
  
  cel_gon = data.frame("goid"=names(cel_list),
                       "ontology"=Ontology( names(cel_list) ),
                       "n"=sapply(cel_list, function(x) length(unique(x))), row.names = NULL )
  
  dbWriteTable(con2, name = "gon_dro", value = drome_gon, overwrite=T, 
                field.types=c(goid="CHAR(10)", ontology="CHAR(2)", n="SMALLINT UNSIGNED"), row.names=F )
  dbWriteTable(con2, name = "gon_cel", value = cel_gon, overwrite=T, 
               field.types=c(goid="CHAR(10)", ontology="CHAR(2)", n="SMALLINT UNSIGNED"), row.names=F )
  
  
  # --> genes meta-data -- TABLE 'genes_[species]'
  library(biomaRt)
  ENS = list()
  ENS$dro_mart = useEnsembl(biomart="ensembl", dataset= "dmelanogaster_gene_ensembl")
  ENS$dro_entrezSymbols = getBM(attributes=c("entrezgene", "flybase_gene_id" , "external_gene_name"), mart=ENS$dro_mart)
  ENS$dro_entrezSymbols = ENS$dro_entrezSymbols[ which(!is.na(ENS$dro_entrezSymbols$entrezgene)), ]
  colnames(ENS$dro_entrezSymbols) = c("entrez","geneid","genename")
  
  ENS$cel_mart = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl")
  ENS$cel_entrezSymbols = getBM(attributes=c("entrezgene", "wormbase_gene" , "external_gene_name"), mart=ENS$cel_mart)
  ENS$cel_entrezSymbols = ENS$cel_entrezSymbols[ which(!is.na(ENS$cel_entrezSymbols$entrezgene)), ]
  colnames(ENS$cel_entrezSymbols) = c("entrez","geneid","genename")
  
  dbWriteTable(con2, name = "genes_dro", value = ENS$dro_entrezSymbols, overwrite=T, 
               field.types=c(species="CHAR(3)", entrez="INT UNSIGNED", geneid="CHAR(14)", genename="VARCHAR(20)"), row.names=F )
  dbWriteTable(con2, name = "genes_cel", value = ENS$cel_entrezSymbols, overwrite=T, 
               field.types=c(species="CHAR(3)", entrez="INT UNSIGNED", geneid="CHAR(14)", genename="VARCHAR(20)"), row.names=F )
  
  
  # create INDICES
  
  dbGetQuery(con2, "CREATE INDEX `idx_goid` ON anno_dro (`goid`)")
  dbGetQuery(con2, "CREATE INDEX `idx_goid` ON anno_cel (`goid`)")
  
  dbGetQuery(con2, "CREATE INDEX `idx_evidence` ON anno_dro (`evidence`)")
  dbGetQuery(con2, "CREATE INDEX `idx_evidence` ON anno_cel (`evidence`)")
  
  dbGetQuery(con2, "CREATE INDEX `idx_entrez` ON anno_dro (`entrez`)")
  dbGetQuery(con2, "CREATE INDEX `idx_entrez` ON anno_cel (`entrez`)")
  
  dbGetQuery(con2, "CREATE INDEX `idx_gon` ON gon_dro (`goid`,`n`,`ontology`)")
  dbGetQuery(con2, "CREATE INDEX `idx_gon` ON gon_cel (`goid`,`n`,`ontology`)")
  
  dbGetQuery(con2, "CREATE INDEX `idx_genes` ON genes_dro (`entrez`)")
  dbGetQuery(con2, "CREATE INDEX `idx_genes` ON genes_cel (`entrez`)")
  
  
  # check if the JOIN types for intended query make sense
  dbGetQuery(con2, 
             "EXPLAIN SELECT anno_dro.goid, anno_dro.entrez, genename FROM anno_dro 
             JOIN gon_dro ON anno_dro.goid = gon_dro.goid 
             JOIN genes_dro ON anno_dro.entrez = genes_dro.entrez
             WHERE n BETWEEN 2 AND 500 
             AND ontology IN ('BP','CC','MF')")
}
















