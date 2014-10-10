library(data.table)
this_input="blast_dbs/viruses_blastdb_content.txt blast_dbs/bact_blastdb_content.txt blast_dbs/euk_blastdb_content.txt blast_dbs/archaea_blastdb_content.txt"
this_input=unlist(strsplit(this_input,split=" "))

str(this_input)
all_content=data.table()
for(a_file in this_input){
	this_content=fread(a_file,sep='\t')
	this_content$file=a_file
	all_content=rbind(all_content,this_content)
}

# known bads 

bad_gis=fread("blast_dbs/bad_bact_headers.txt",sep="|")$V2

#detect length outliers 

setnames(all_content,c("gi","acc",'acc2',"description","length","file"))
all_content$file=factor(all_content$file)

save(all_content,file="all_db_content.RData")

# load("all_db_content.RData")
all_content[,domain:=""]
all_content[,to_remove:=F]
all_content[file=="blast_dbs/archaea_blastdb_content.txt",domain:="archaea"]
all_content[file=="blast_dbs/bact_blastdb_content.txt",domain:="bact"]
all_content[file=="blast_dbs/euk_blastdb_content.txt",domain:="euk"]
all_content[file=="blast_dbs/viruses_blastdb_content.txt",domain:="viruses"]



all_content[order(length)]


by(all_content,all_content$file,function(df)summary(df$length))
by(all_content,all_content$file,function(df)quantile(df$length,probs=c(0.95)))

# check all megabases genomes 

all_content[length>=10e6][order(file)]

# check largest 95% 

all_content[file=="blast_dbs/archaea_blastdb_content.txt" & length>=154599]


# check multiple appearance 

bad_gis=unique(c(all_content[,.N,by=list(file,gi)][,.N,by=gi][N>1,gi],bad_gis))


# where from ? 

all_content[gi %in% bad_gis,.N,by=file]


# filter first 

all_content=all_content[!(gi %in% bad_gis)]

by(all_content,all_content$file,function(df)summary(df$length))
by(all_content,all_content$file,function(df)quantile(df$length,probs=c(0.95)))


# check the largest per file 

all_content[,size_rank:=as.numeric(rank(-1*length,ties.method="first")),by=file]
all_content[size_rank<=4,list(file,gi,description,length,size_rank)]

# mosquitos!

mosquitos_gi=all_content[file=="blast_dbs/bact_blastdb_content.txt" & grepl("Anopheles gambiae",description),gi]
bad_gis=unique(c(bad_gis,mosquitos_gi))
length(bad_gis)

all_content=all_content[!(gi %in% bad_gis)]




# check the largest per file 

all_content[,size_rank:=as.numeric(rank(-1*length,ties.method="first")),by=file]
all_content[size_rank<=4,list(file,gi,description,length,size_rank)]

# seems ok 


# viruses in bact ? 


all_content[file=="blast_dbs/bact_blastdb_content.txt"][grep('virus',description)] #nothing

all_content[domain=="bact" & grepl("saccharomyces",description,ignore.case=T),to_remove:=T,by=.I]

# all_content[domain=="viruses" & grepl("saccharomyces",description,ignore.case=T)]
# all_content[domain=="viruses" & grepl("saccharomyces",description,ignore.case=T),to_remove:=F,by=.I]

all_content[domain=="archaea" & grepl("saccharomyces",description,ignore.case=T)]

all_content[domain=="archaea" & grepl("saccharomyces",description,ignore.case=T),to_remove:=T,by=.I]

all_content[to_remove==T]


with_chr=all_content[grep("chromosome",description)]

write.table(with_chr,file="with_chr.txt")


# add the SACE 

bad_gis=unique(c(bad_gis,all_content[to_remove==T,gi]))

write.table(bad_gis,file="blast_dbs/bad_gis.txt")