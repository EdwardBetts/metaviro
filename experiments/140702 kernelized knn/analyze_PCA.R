library(ggplot2)
library(reshape2)
library(gridExtra)
library(data.table)

pca_coord=fread("pca_coord_3mers.csv")

ggplot(pca_coord,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)
ggsave("pca_coord_3mers.pdf")


ggplot(pca_coord,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)+facet_wrap(~class)
ggsave("pca_coord_3mers_by_class.pdf")


ggplot(pca_coord,aes(x=coord_1,y=coord_2,colour=GC))+geom_point()+facet_wrap(~class)
ggsave("pca_coord_3mers_gc.pdf")

ggplot(pca_coord,aes(x=coord_1,y=coord_2,colour=sequence_length))+geom_point()+facet_wrap(~class)
ggsave("pca_coord_3mers_length.pdf")


ggplot(pca_coord,aes(x=sequence_length,y=GC,colour=coord_1))+geom_point()+facet_wrap(~class)
ggplot(pca_coord,aes(x=sequence_length,y=GC,colour=coord_2))+geom_point()+facet_wrap(~class)



# on the length normalized data 

pca_coord_norm= fread("pca_coord_length_normalized_3mers.csv")

ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)+geom_density2d()

ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)+facet_wrap(~class)
ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=GC))+geom_point()


# We try a GC content normalization of the counts 

ggplot(pca_coord_norm,aes(x=GC,y=GCG,colour=class))+geom_point(alpha=0.2)+geom_density2d()
ggplot(pca_coord_norm,aes(x=GC,y=GCG,colour=class))+geom_point(alpha=0.2)+geom_smooth()

ggplot(pca_coord_norm,aes(x=GC,y=GGG_x,colour=class))+geom_point(alpha=0.2)+geom_smooth()
ggplot(pca_coord_norm,aes(x=GC,y=GGG_y,colour=class))+geom_point(alpha=0.2)+geom_smooth()

g1=ggplot(pca_coord_norm,aes(x=GC,y=TTT_x,colour=class))+geom_point(alpha=0.2)+geom_smooth()
g2=ggplot(pca_coord_norm,aes(x=GC,y=TTT_y,colour=class))+geom_point(alpha=0.2)+geom_smooth()
grid.arrange(g1,g2)

g1=ggplot(pca_coord_norm,aes(x=GC,y=ATT_x,colour=class))+geom_point(alpha=0.2)+geom_smooth()
g2=ggplot(pca_coord_norm,aes(x=GC,y=ATT_y,colour=class))+geom_point(alpha=0.2)+geom_smooth()
grid.arrange(g1,g2)

ggplot(pca_coord_norm,aes(x=GC,y=ATT,colour=class))+geom_point(alpha=0.2)+geom_smooth()
ggplot(pca_coord_norm,aes(x=GC,y=CCC,colour=class))+geom_point(alpha=0.2)+geom_smooth()
ggplot(pca_coord_norm,aes(x=GC,y=CGC,colour=class))+geom_point(alpha=0.2)+geom_smooth()


pca_coord_norm= fread("pca_coord_length_gc_normalized_3mers.csv")

ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)+geom_density2d()

ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=class))+geom_point(alpha=0.5)+facet_wrap(~class,scale="free")
ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=GC))+geom_point()+geom_density2d()+facet_wrap(~class)
ggplot(pca_coord_norm,aes(x=coord_2,y=coord_3,colour=GC))+geom_point()+geom_density2d()+facet_wrap(~class)
ggplot(pca_coord_norm,aes(x=coord_2,y=coord_4,colour=GC))+geom_point()+geom_density2d()+facet_wrap(~class)
ggplot(pca_coord_norm,aes(x=coord_3,y=coord_4,colour=GC))+geom_point()+geom_density2d()+facet_wrap(~class)

ggplot(pca_coord_norm,aes(x=coord_1,y=coord_2,colour=sequence_length))+geom_point()+geom_density2d()



# if we use the other dimensions, can we correct for the GC content variation ? 