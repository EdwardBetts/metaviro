suppressMessages(library(data.table))
gi_taxid_nucl=fread("gi_taxid_nucl.dmp") # Very long !
setnames(gi_taxid_nucl,c("gi","TaxID"))
save(gi_taxid_nucl,file="gi_mapping/gi_taxid_nucl.RData")

