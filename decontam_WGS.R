##################### BATCH CORRECTION PLOTS - BEFORE
#################### LOAD EXTRACTION DATA
WGS_extractSites<-
  read_delim("extra_data/DNA-extraction_batch.txt")%>%
  left_join(broad_big_list)%>%
  mutate(extraction_site1=
           if_else(`Created By`=="Nationwide"|`Created By`=="CGR",
                   `Created By`,
                   Study),
         extraction_site2=
           case_when(`Created By`=="Nationwide" ~ paste0("Nationwide-",`NW Batch #`),
                     `Created By`=="CGR" ~ "CGR",
                     T ~ Study))%>%
  select(Barcode=`TZ Barcode`,extraction_site1,extraction_site2)%>%
  filter(Barcode %in% wgs_nopublic$Barcode,
         Barcode %in% colnames(wgs_sherlockOnly))%>%
  select(Barcode,extraction_site1,extraction_site2)%>%
  unique()

WGS_extractColor<-SBScolor[1:nrow(count(WGS_extractSites,extraction_site2))]
names(WGS_extractColor)<-
  count(WGS_extractSites,extraction_site2)%>%
  pull(extraction_site2)

WGS_batchBefore<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>2)>0.05)%>%
  select(any_of(wgs_nopublic$Barcode))%>%
  #select(any_of(natgen_232_barcodes))%>%
  t()%>%
  avgdist(sample = 100,iterations=50)
WGS_batchBefore_cmdscale<-cmdscale(WGS_batchBefore,list. = T,eig=T,k=10)
WGS_batchBefore_cmdscale$eig<-(WGS_batchBefore_cmdscale$eig/sum(WGS_batchBefore_cmdscale$eig))*100
loadings<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  filter(rowMeans(.>10)>0.05)%>%
  column_to_rownames("tax_id")%>%
  .[,labels(WGS_batchBefore)]%>%
  apply(MARGIN = 2,FUN = function(x){x/sum(x)})%>%
  apply(MARGIN = 1,FUN = function(x){x/sum(x)})%>%
  t()%>%
  `%*%`(WGS_batchBefore_cmdscale$points)%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%
  mutate(magnitude=sqrt(PC1^2+PC2^2))%>%
  arrange(-magnitude)
WGS_batchBefore_plot<-
  WGS_batchBefore_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(WGS_extractSites,by=c("SampleID"="Barcode"))%>%
  left_join(wgs_full_annotations,by=c("SampleID"="Barcode"))%>%
  mutate(natgen_232 = if_else(SampleID %in% natgen_232_barcodes,"Zhang et al 2021","This study"),
         extraction_site=case_when(extraction_site2=="Yale"~"Connecticut",
                                   extraction_site2=="Harvard"~"Massachusetts",
                                   extraction_site2=="INCAN"~"Mexico City",
                                   extraction_site2=="Roswell Park"~"New York",
                                   extraction_site2=="Mayo"~"Minnesota",
                                   extraction_site2=="Moffitt"~"Florida",
                                   TRUE~extraction_site2))%>%
  ggplot(aes(x=PC1,y=PC2,group=extraction_site2,label=NA))+
  stat_ellipse(aes(color=extraction_site2),show.legend = F)+
  geom_point(aes(fill=extraction_site2,shape=natgen_232),color="gray80")+
  labs(x=paste0("PC1 (",round(WGS_batchBefore_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(WGS_batchBefore_cmdscale$eig[2],1),"%)"),
       fill="DNA extraction",
       shape="Study")+
  guides(fill = guide_legend(ncol=2,override.aes=list(shape = 21)),
         shape=guide_legend(override.aes = list(fill="black")))+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=WGS_extractColor)+
  scale_color_manual(values=WGS_extractColor)#+
#geom_point(data=loadings[1:20,],aes(x=PC1,y=PC2),inherit.aes = F)+
#geom_text_repel(data=loadings[1:20,],aes(x=PC1,y=PC2,label=name),inherit.aes = F)
WGS_batchBefore_plot

#################### BATCH CORRECT
#
natgen_232<-read_delim("extra_data/NatGen_232.txt",delim = "\t")%>%drop_na(Subject)
natgen_232_barcodes<-c(natgen_232$Tumor_Barcode,natgen_232$Normal_Barcode)

#################### COMBAT
WGS_combat_in<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% (unified_taxonomy%>%
                        filter(str_detect(taxonomy,"Bacteria"))%>%
                        filter(type=="genus")%>%
                        pull(tax_id)))%>%
  column_to_rownames("tax_id")%>%
  select(any_of(wgs_nopublic$Barcode))%>%
  select(-any_of(natgen_232_barcodes))%>%
  select_if(colSums(.)>100)%>%
  filter(rowMeans(.>0)>0.01)
WGS_combat_in_phylum<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% (unified_taxonomy%>%
                        filter(str_detect(taxonomy,"Bacteria"))%>%
                        filter(type=="phylum")%>%
                        pull(tax_id))
  )%>%
  column_to_rownames("tax_id")%>%
  select(any_of(wgs_nopublic$Barcode))%>%
  select(-any_of(natgen_232_barcodes))%>%
  select_if(colSums(.)>100)%>%
  filter(rowMeans(.>0)>0.01)

wgs.adj.var<-
  annotations_withClinical%>%
  filter(Barcode %in% colnames(WGS_combat_in))%>%
  pull(Study,Barcode)
wgs.adj.var2<-
  WGS_extractSites%>%
  filter(Barcode %in% colnames(WGS_combat_in))%>%
  filter(!extraction_site2 %in% c("Yale","Valencia","Nice"))%>%
  pull(extraction_site2,Barcode)
wgs.bio.var<-
  annotations_withClinical%>%
  filter(Barcode %in% colnames(WGS_combat_in))%>%
  filter(!Study %in% c("Yale","Valencia","Nice"))%>%
  pull(Sample_Source,Barcode)

wgs_decontamd<-
  ComBat_seq(as.matrix(WGS_combat_in)[,names(wgs.adj.var2)],
             batch = wgs.adj.var2
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")%>%
  filter(!tax_id %in% salter_list_nonHuman$tax_id)

wgs_decontamd_phylum<-
  ComBat_seq(as.matrix(WGS_combat_in_phylum)[,names(wgs.adj.var2)],
             batch = wgs.adj.var2
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")%>%
  filter(!tax_id %in% salter_list_nonHuman$tax_id)

##################### BATCH CORRECTION PLOTS - AFTER
WGS_batchAfter<-
  wgs_decontamd%>%
  select(-tax_id)%>%
  as.data.frame()%>%
  replace(is.na(.),0)%>%
  t()%>%floor()%>%
  avgdist(sample = 100,iterations=50)
WGS_batchAfter_cmdscale<-
  cmdscale(WGS_batchAfter,list. = T,eig=T,k=10)
WGS_batchAfter_cmdscale$eig<-
  (WGS_batchAfter_cmdscale$eig/sum(WGS_batchAfter_cmdscale$eig))*100
loadings2<-
  wgs_decontamd%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>0)>0.01)%>%
  replace(is.na(.),0)%>%
  filter(rownames(.) %in% all_bact_genera)%>%
  .[,labels(WGS_batchAfter)]%>%
  apply(MARGIN = 2,FUN = function(x){x/sum(x)})%>%
  apply(MARGIN = 1,FUN = function(x){x/sum(x)})%>%
  t()%>%
  `%*%`(WGS_batchAfter_cmdscale$points)%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%
  mutate(magnitude=sqrt(PC1^2+PC2^2))%>%
  arrange(-magnitude)
WGS_batchAfter_plot<-
  WGS_batchAfter_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(WGS_extractSites,by=c("SampleID"="Barcode"))%>%
  left_join(wgs_full_annotations,by=c("SampleID"="Barcode"))%>%
  mutate(natgen_232 = SampleID %in% natgen_232_barcodes)%>%
  mutate(extraction_site=case_when(extraction_site2=="Yale"~"Connecticut",
                   extraction_site2=="Harvard"~"Massachusetts",
                   extraction_site2=="INCAN"~"Mexico City",
                   extraction_site2=="Roswell Park"~"New York",
                   extraction_site2=="Mayo"~"Minnesota",
                   extraction_site2=="Moffitt"~"Florida",
                   TRUE~extraction_site2))%>%
  ggplot(aes(x=PC1,y=PC2,group=extraction_site2,label=NA))+
  stat_ellipse(aes(color=extraction_site2),show.legend = F)+
  geom_point(aes(fill=extraction_site2),shape=21,color="gray80")+
  labs(x=paste0("PC1 (",round(WGS_batchAfter_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(WGS_batchAfter_cmdscale$eig[2],1),"%)"),
       fill="DNA Extraction")+
  scale_fill_manual(values=WGS_extractColor)+
  scale_color_manual(values=WGS_extractColor)+
  theme(legend.text = element_text(size=8))+
  guides(fill=guide_legend(ncol=2))#+
  #geom_point(data=loadings2[1:20,],aes(x=PC1,y=PC2),inherit.aes = F)+
  #geom_text_repel(data=loadings2[1:20,],aes(x=PC1,y=PC2,label=name),inherit.aes = F)
WGS_batchAfter_plot

adonis2(WGS_batchAfter~extraction_site2,
       data = WGS_extractSites%>%
         filter(Barcode %in% labels(WGS_batchAfter)))

##################### DECONTAMINATE RAW READS
wgs_rawCounts_decontamd<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  select(any_of(wgs_nopublic$Barcode))%>%
  select(-any_of(natgen_232_barcodes))%>%
  filter(rowMeans(.>0)>0.01)%>%
  select_if(colSums(.)>100)%>%
  decontaminate_up(c(salter_list_nonHuman$tax_id),
                   "genus",
                   tax_table)%>%
  .[[1]]%>%
  rownames_to_column("tax_id")%>%
  replace(is.na(.),0)

#################### NATGEN 232 COMBAT
NG232_combat_in<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% (unified_taxonomy%>%
                        filter(str_detect(taxonomy,"Bacteria"))%>%
                        filter(type=="genus")%>%
                        pull(tax_id)))%>%
  column_to_rownames("tax_id")%>%
  select(any_of(natgen_232_barcodes))%>%
  select_if(colSums(.)>100)%>%
  filter(rowMeans(.>0)>0.01)
NG232_combat_in_phylum<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% (unified_taxonomy%>%
                        filter(str_detect(taxonomy,"Bacteria"))%>%
                        filter(type=="phylum")%>%
                        pull(tax_id))
  )%>%
  column_to_rownames("tax_id")%>%
  select(any_of(natgen_232_barcodes))%>%
  select_if(colSums(.)>100)%>%
  filter(rowMeans(.>0)>0.01)
NG232.adj.var2<-
  WGS_extractSites%>%
  filter(Barcode %in% colnames(NG232_combat_in))%>%
  pull(extraction_site2,Barcode)
NG232.bio.var<-
  annotations_withClinical%>%
  filter(Barcode %in% colnames(NG232_combat_in))%>%
  pull(Sample_Source,Barcode)

NG232_decontamd<-
  ComBat_seq(as.matrix(NG232_combat_in)[,names(NG232.adj.var2)],
             batch = NG232.adj.var2
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")%>%
  filter(!tax_id %in% salter_list_nonHuman$tax_id)
NG232_decontamd_phylum<-
  ComBat_seq(as.matrix(NG232_combat_in_phylum)[,names(NG232.adj.var2)],
             batch = NG232.adj.var2
  )%>%
  as.data.frame()%>%
  rownames_to_column("tax_id")%>%
  filter(!tax_id %in% salter_list_nonHuman$tax_id)

NG232_batchAfter<-
  NG232_decontamd%>%
  select(-tax_id)%>%
  as.data.frame()%>%
  replace(is.na(.),0)%>%
  t()%>%floor()%>%
  avgdist(sample = 100,iterations=50)
NG232_batchAfter_cmdscale<-
  cmdscale(NG232_batchAfter,list. = T,eig=T,k=10)
NG232_batchAfter_cmdscale$eig<-
  (NG232_batchAfter_cmdscale$eig/sum(NG232_batchAfter_cmdscale$eig))*100
loadings2<-
  NG232_decontamd%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>0)>0.01)%>%
  replace(is.na(.),0)%>%
  filter(rownames(.) %in% all_bact_genera)%>%
  .[,labels(NG232_batchAfter)]%>%
  apply(MARGIN = 2,FUN = function(x){x/sum(x)})%>%
  apply(MARGIN = 1,FUN = function(x){x/sum(x)})%>%
  t()%>%
  `%*%`(NG232_batchAfter_cmdscale$points)%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%
  mutate(magnitude=sqrt(PC1^2+PC2^2))%>%
  arrange(-magnitude)
NG232_batchAfter_plot<-
  NG232_batchAfter_cmdscale$points%>%
  as.data.frame()%>%
  `colnames<-`(paste0("PC",seq(1:10)))%>%
  rownames_to_column("SampleID")%>%
  left_join(WGS_extractSites,by=c("SampleID"="Barcode"))%>%
  left_join(wgs_full_annotations,by=c("SampleID"="Barcode"))%>%
  mutate(natgen_232 = SampleID %in% natgen_232_barcodes)%>%
  mutate(extraction_site=case_when(extraction_site2=="Yale"~"Connecticut",
                                   extraction_site2=="Harvard"~"Massachusetts",
                                   extraction_site2=="INCAN"~"Mexico City",
                                   extraction_site2=="Roswell Park"~"New York",
                                   extraction_site2=="Mayo"~"Minnesota",
                                   extraction_site2=="Moffitt"~"Florida",
                                   TRUE~extraction_site2))%>%
  ggplot(aes(x=PC1,y=PC2,group=extraction_site2,label=NA))+
  stat_ellipse(aes(color=extraction_site2),show.legend = F)+
  geom_point(aes(fill=extraction_site2),shape=24,color="gray80")+
  labs(x=paste0("PC1 (",round(NG232_batchAfter_cmdscale$eig[1],1),"%)"),
       y=paste0("PC2 (",round(NG232_batchAfter_cmdscale$eig[2],1),"%)"),
       fill="DNA Extraction")+
  scale_fill_manual(values=WGS_extractColor)+
  scale_color_manual(values=WGS_extractColor)+
  theme(legend.text = element_text(size=8))
  #guides(fill=guide_legend(ncol=2))#+
  #geom_point(data=loadings2[1:20,],aes(x=PC1,y=PC2),inherit.aes = F)+
  #geom_text_repel(data=loadings2[1:20,],aes(x=PC1,y=PC2,label=name),inherit.aes = F)
NG232_batchAfter_plot


ng232_rawCounts_decontamd<-
  wgs_sherlockOnly%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  select(any_of(wgs_nopublic$Barcode))%>%
  select(any_of(natgen_232_barcodes))%>%
  filter(rowMeans(.>0)>0.01)%>%
  select_if(colSums(.)>100)%>%
  decontaminate_up(c(salter_list_nonHuman$tax_id),
                   "genus",
                   tax_table)%>%
  .[[1]]%>%
  rownames_to_column("tax_id")%>%
  replace(is.na(.),0)



