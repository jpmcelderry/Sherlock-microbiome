###### RNA diversity calculations
# alpha-diversity
rna_diver_map<-
  map_dfr(c(500,750,1000,1500),
          function(x) map_dfr(1:100,~(rna_combatd_decontamd%>%
                                        filter(tax_id %in% all_bact_genera)%>%
                                        column_to_rownames("tax_id")%>%
                                        t()%>%as.data.frame()%>%
                                        filter(rowSums(.)>x)%>%
                                        rrarefy(sample = x)%>%
                                        diversity(index="shannon")))%>%
            apply(2,median,na.rm=TRUE),.id = "sample_depth")
RNAseq_diversity<-
  rna_diver_map%>%
  select(-1)%>%
  t()%>%as.data.frame()%>%
  `colnames<-`(c("500","750","1000","1500"))%>%
  rownames_to_column("RNAseq_SampleID")
# richness
rna_richness<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  t()%>%as.data.frame()%>%
  rarefy(sample = c(100,200,500,750,1000,1500))%>%
  as.data.frame()%>%
  mutate(.before = N100,RNAseq_SampleID=colnames(rna_combatd_decontamd)[-1])
# beta-diversity
rna_combatd_avgdist<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  filter(rowMeans(.>10)>0.01)%>%
  t()%>%
  avgdist(500,iterations=50)

###### 16S diversity calculations
# alpha-diversity
s16_diver_map<-
  map_dfr(c(100,200,250,500,750),
          function(x) map_dfr(1:100,~(s16_scrubbed%>%
                                        filter(tax_id %in% all_bact_genera)%>%
                                        column_to_rownames("tax_id")%>%
                                        t()%>%as.data.frame()%>%
                                        filter(rowSums(.)>=x)%>%
                                        rrarefy(sample = x)%>%
                                        diversity(index="shannon")))%>%
            apply(2,median,na.rm=TRUE),.id = "sample_depth")
s16_diversity<-
  s16_diver_map%>%
  select(-1)%>%
  t()%>%as.data.frame()%>%
  `colnames<-`(c("100","200","250","500","750"))%>%
  rownames_to_column("SampleID")
# richness
s16_richness<-
  s16_scrubbed%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  t()%>%as.data.frame()%>%
  rarefy(sample = c(100,200,250,500,750))%>%
  as.data.frame()%>%
  mutate(.before = N100,SampleID=colnames(s16_scrubbed)[-1])
# beta-diversity
s16_avgdist<-
  s16_scrubbed%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  filter(rowMeans(.>0)>0.05)%>%
  t()%>%
  avgdist(250,iterations=50)

###### WGS diversity calculations
# alpha-diversity
wgs_diver_map<-
  map_dfr(c(100,250,500,750),
          function(x) map_dfr(1:100,~(wgs_decontamd%>%
                                        column_to_rownames("tax_id")%>%
                                        t()%>%as.data.frame()%>%
                                        filter(rowSums(.)>x)%>%
                                        rrarefy(sample = x)%>%
                                        diversity(index="shannon")))%>%
            apply(2,median,na.rm=TRUE),.id = "sample_depth")
WGS_diversity<-
  wgs_diver_map%>%
  select(-1)%>%
  t()%>%as.data.frame()%>%
  `colnames<-`(c("100","250","500","750"))%>%
  rownames_to_column("Barcode")
# richness
WGS_richness<-
  wgs_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  t()%>%as.data.frame()%>%
  rarefy(sample = c(100,200,250,500,750))%>%
  as.data.frame()%>%
  mutate(.before = N100,Barcode=colnames(wgs_decontamd)[-1])
# beta-diversity
wgs_avgdist<-
  wgs_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  filter(rowMeans(.>0)>0.01)%>%
  t()%>%
  avgdist(100,iterations = 50)

###### NG232 diversity calculations
# alpha-diversity
NG232_diver_map<-
  map_dfr(c(100,250,500,750),
          function(x) map_dfr(1:100,~(NG232_decontamd%>%
                                        column_to_rownames("tax_id")%>%
                                        t()%>%as.data.frame()%>%
                                        filter(rowSums(.)>x)%>%
                                        rrarefy(sample = x)%>%
                                        diversity(index="shannon")))%>%
            apply(2,median,na.rm=TRUE),.id = "sample_depth")
NG232_diversity<-
  NG232_diver_map%>%
  select(-1)%>%
  t()%>%as.data.frame()%>%
  `colnames<-`(c("100","250","500","750"))%>%
  rownames_to_column("Barcode")
# richness
NG232_richness<-
  NG232_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  t()%>%as.data.frame()%>%
  rarefy(sample = c(100,200,250,500,750))%>%
  as.data.frame()%>%
  mutate(.before = N100,Barcode=colnames(NG232_decontamd)[-1])
# beta-diversity
NG232_avgdist<-
  NG232_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  select(-tax_id)%>%
  filter(rowMeans(.>0)>0.01)%>%
  t()%>%
  avgdist(100,iterations = 50)