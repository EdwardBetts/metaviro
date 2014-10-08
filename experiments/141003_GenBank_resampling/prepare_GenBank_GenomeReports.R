# Viruses 
load("gi_mapping/gi_taxid_nucl.RData")

viruses_annotations=fread("../../data/GENOME_REPORTS/viruses.txt",na.strings=c("-","NA"))
setkey(viruses_annotations,"TaxID")
# Remove duplicate TaxID 
setnames(viruses_annotations,make.names(colnames(viruses_annotations)))

# flag duplicates 
viruses_annotations[TaxID %in% viruses_annotations[,.N,by=TaxID][N>1][,TaxID],is_duplicate:=T]
viruses_annotations[is_duplicate==T]
viruses_annotations[is.na(GC.) & (is_duplicate==T),to_remove:=T] 
# flag for removal incomplete entries from the duplicates 
# viruses_annotations[to_remove==T]
viruses_annotations=viruses_annotations[is.na(to_remove)]
# Reflag remaining duplicates 
viruses_annotations[,is_duplicate:=F]
viruses_annotations[TaxID %in% viruses_annotations[,.N,by=TaxID][N>1][,TaxID],is_duplicate:=T]
viruses_annotations[is_duplicate==T][order(TaxID)]
# Take randomly one of the entries per each TaxID
viruses_annotations[is_duplicate==T,rank:=rank(BioProject.ID),by=TaxID]
viruses_annotations[is_duplicate==T & rank>1,to_remove:=T]
viruses_annotations=viruses_annotations[is.na(to_remove)]

# Last flight, should be unique now 
viruses_annotations[,is_duplicate:=F]
viruses_annotations[TaxID %in% viruses_annotations[,.N,by=TaxID][N>1][,TaxID],is_duplicate:=T]
viruses_annotations[is_duplicate==T][order(TaxID)]

# Clear accessory columns 
viruses_annotations$is_duplicate<-NULL
viruses_annotations$to_remove<-NULL
viruses_annotations$rank<-NULL

viruses_taxId2gi=gi_taxid_nucl[ TaxID %in% viruses_annotations$TaxID]

setnames(viruses_annotations,c("X.Organism.Name"),c("Organism.Name"))
setnames(viruses_annotations,c("Size..Kb."),c("Size"))

GenBank_annotations=viruses_annotations
GenBank_taxID2gi=viruses_taxId2gi
# viruses_annotations_gi=merge(viruses_annotations,viruses_taxId2gi)
save(GenBank_annotations,GenBank_taxID2gi,file="gi_mapping/viruses_annotations_gi.RData")


# Bacteria 

bact_annotations=fread("../../data/GENOME_REPORTS/prokaryotes.txt",na.strings=c("-","NA"))
bact_annotations[is.na(Group)]
bact_annotations[is.na(SubGroup)]
# We provide an order of preferences over the possible status 

bact_annotations$Status=factor(bact_annotations$Status,levels=c("Complete Genome","Complete","Chromosome","Chromosome with gaps","Scaffold","Contig"))
levels(bact_annotations$Status)
setkey(bact_annotations,"TaxID")
setnames(bact_annotations,make.names(colnames(bact_annotations)))

bact_annotations[,is_duplicate:=F]
bact_annotations[TaxID %in% bact_annotations[,.N,by=TaxID][N>1][,TaxID],is_duplicate:=T]
# bact_annotations[is_duplicate==T]
bact_annotations[is_duplicate==T,rank:=rank(Status,ties.method="random"),by=TaxID]
bact_annotations[is_duplicate==T & rank==1]
bact_annotations[TaxID==210][order(rank)]

# check that group and subgroups agrees over duplicate BioProject 
compared_dupl=dcast.data.table(bact_annotations[is_duplicate==T],TaxID~rank,value.var="Status")
setnames(compared_dupl,make.names(colnames(compared_dupl)))
compared_dupl[,list(TaxID,X1,X2)]
compared_dupl[X1!=X2] #all agrees on groups, SubGroups
# Disagree on  Status...

bact_annotations[rank>1,to_remove:=T]
bact_annotations=bact_annotations[is.na(to_remove)]
bact_annotations[,is_duplicate:=F]
bact_annotations[TaxID %in% bact_annotations[,.N,by=TaxID][N>1][,TaxID],is_duplicate:=T]
bact_annotations[is_duplicate==T]

bact_annotations$is_duplicate<-NULL
bact_annotations$to_remove<-NULL
bact_annotations$rank<-NULL

setnames(bact_annotations,c("Size..Mb."),c("Size"))

bact_taxId2gi=gi_taxid_nucl[ TaxID %in% bact_annotations$TaxID]

GenBank_annotations=bact_annotations
GenBank_taxID2gi=bact_taxId2gi


# bact_annotations_gi=merge(bact_annotations,gi_taxid_nucl)
save(GenBank_annotations,GenBank_taxID2gi,file="gi_mapping/bact_annotations_gi.RData")
# load("gi_mapping/bact_annotations_gi.RData")
