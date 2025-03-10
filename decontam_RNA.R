# library(sva)

# filter to only frozen (no FFPE) samples sent to Nationwide
rna.adj.vars<-
  rna_nationwide_frozen_Samples%>%
  pull(Site_REF,RNAseq_SampleID)
# rna.adj.vars2<-
#   rna_nationwide_frozen_Samples%>%
#   pull(batch,RNAseq_SampleID)
rna.bio.vars<-
  rna_nationwide_frozen_Samples%>%
  pull(Type,RNAseq_SampleID)

rna_combat_in<-
  rna_NW_FF%>%
  select(tax_id,any_of(names(rna.adj.vars)))%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)
rna_combat_in[rna_combat_in<2]<-0
rna_combat_in<-
  rna_combat_in%>%
  filter(rowMeans(.>0)>0.01)

# rna_PLSDA_out<-
#   PLSDA_batch(rna_plsda_in,
#               Y.trt = rna_bio.var,
#               Y.bat = rna.adj.vars,
#               ncomp.trt = nlevels(as.factor(rna.bio.vars))-1,
#               ncomp.bat = nlevels(as.factor(rna.adj.vars))-1,
#               balance=FALSE)

rna_combatd_decontamd<-
  rna_combat_in%>%
  filter(rownames(.) %in% (unified_taxonomy%>%
                             filter(str_detect(taxonomy,"Bacteria"))%>%
                             filter(type=="genus")%>%
                             pull(tax_id)))%>%
  as.matrix()%>%
  ComBat_seq(group=rna.bio.vars,
             batch = rna.adj.vars
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")%>%
  filter(!tax_id %in% salter_list_nonHuman$tax_id)

rna_combatd_decontamd_phylum<-
  rna_combat_in%>%
  filter(rownames(.) %in% (unified_taxonomy%>%
                             filter(str_detect(taxonomy,"Bacteria"))%>%
                             filter(type=="phylum")%>%
                             pull(tax_id)))%>%
  as.matrix()%>%
  ComBat_seq(group=rna.bio.vars,
             batch = rna.adj.vars
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")

rna_rawCounts_decontamd<-
  decontaminate_up(rna_combat_in%>%
                     rownames_to_column("tax_id"),
                   contaminants = salter_list_nonHuman$tax_id,
                   decontam_level = "genus",
                   taxa_table = tax_table)%>%
  .[[1]]%>%
  rownames_to_column("tax_id")

## Before batch correction
RNA_batchBefore<-
  rna%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>=5)>0.05)%>%
  select(any_of(colnames(rna_combatd_decontamd)))%>%
  t()%>%
  avgdist(sample = 500,iterations=50)
RNA_batchBefore_cmdscale<-cmdscale(RNA_batchBefore,list. = T,eig=T,k=10)
RNA_batchBefore_cmdscale$eig<-(RNA_batchBefore_cmdscale$eig/sum(RNA_batchBefore_cmdscale$eig))*100
RNA_batchBefore_plot<-
  RNA_batchBefore_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(rna.adj.vars_setup,by=c("SampleID"="RNAseq_SampleID"))%>%
  mutate(Site_REF=case_when(Site_REF=="Yale"~"Connecticut",
                            Site_REF=="Harvard"~"Massachusetts",
                            Site_REF=="INCAN"~"Mexico City",
                            Site_REF=="Roswell Park"~"New York",
                            Site_REF=="Mayo"~"Minnesota",
                            Site_REF=="Moffitt"~"Florida",
                            TRUE~Site_REF))%>%
  ggplot(aes(x=PC1,y=PC2,group=Site_REF,label=NA))+
  stat_ellipse(aes(color=Site_REF),show.legend = F)+
  geom_point(aes(fill=Site_REF),shape=21,color="gray80",show.legend = F)+
  labs(x=paste0("PC1 (",round(RNA_batchBefore_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(RNA_batchBefore_cmdscale$eig[2],1),"%)"))+
  scale_fill_manual(values=rev(SBScolor))+
  scale_color_manual(values=rev(SBScolor))
RNA_batchBefore_plot

# After batch correction
RNA_batchAfter<-
  rna_combatd_decontamd%>%
  column_to_rownames("tax_id")%>%
  floor()%>%
  t()%>%as.data.frame()%>%select(-any_of(salter_list_nonHuman$tax_id))%>%
  avgdist(sample = 500,iterations=50)
# RNA_batchAfter<-
#   rna_PLSDA_out$X.nobatch%>%
#   as.data.frame()%>%
#   vegdist(method="euclidean")
RNA_batchAfter_cmdscale<-
  cmdscale(RNA_batchAfter,list. = T,eig=T,k=10)
RNA_batchAfter_cmdscale$eig<-
  (RNA_batchAfter_cmdscale$eig/sum(RNA_batchAfter_cmdscale$eig))*100
RNA_batchAfter_plot<-
  RNA_batchAfter_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(rna_nationwide_frozen_Samples,by=c("SampleID"="RNAseq_SampleID"))%>%
  mutate_at("batch",str_replace,"Sherlock","Unavailable")%>%
  mutate(Site_REF=case_when(Site_REF=="Yale"~"Connecticut",
                            Site_REF=="Harvard"~"Massachusetts",
                            Site_REF=="INCAN"~"Mexico City",
                            Site_REF=="Roswell Park"~"New York",
                            Site_REF=="Mayo"~"Minnesota",
                            Site_REF=="Moffitt"~"Florida",
                            TRUE~Site_REF))%>%
  ggplot(aes(x=PC1,y=PC2,group=Site_REF,label=NA))+
  stat_ellipse(aes(color=Site_REF),show.legend = F)+
  geom_point(aes(fill=Site_REF),shape=21,color="gray80")+
  labs(x=paste0("PC1 (",round(RNA_batchAfter_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(RNA_batchAfter_cmdscale$eig[2],1),"%)"),
       fill="Study Site")+
  scale_fill_manual(values=rev(SBScolor))+
  scale_color_manual(values=rev(SBScolor))+
  theme(legend.text = element_text(size=8))+
  guides(fill=guide_legend(ncol=2))
RNA_batchAfter_plot
