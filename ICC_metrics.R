# 
# WGS_tumorOnly<-
#   wgs_decontamd%>%
#   select(any_of(c("tax_id",(annotations_withClinical%>%
#                               filter(Sample_Source=="Tumor")%>%
#                               pull(Barcode)))))
# 

######### FIG 2
#### PHYLUM RELATIVE ABUNDANCE CORRELATION
rna_wgs_key<-
  paired_rna_wgs%>%
  pivot_wider(names_from = "experiment",values_from = "SampleID")%>%
  drop_na()%>%
  unnest_longer("RNA")%>%
  unnest_longer("WGS")
wgs_16s_key<-
  paired_wgs_16s%>%
  pivot_wider(names_from = "experiment",values_from = "SampleID")%>%
  drop_na()%>%
  mutate_at(c("WGS","16s"),as.character)

rna_phylum_CLR<-
  rna_combatd_decontamd_phylum%>%
  filter(tax_id %in% all_bact_phyla)%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  select_if(colSums(.)>250)%>%
  t()%>%
  clr(.+.1)%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  pivot_longer(-1,names_to = "tax_id",values_to = "RNA_Reads")
s16_phylum_CLR<-
  s16_scrubbed%>%
  filter(tax_id %in% all_bact_phyla)%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  select_if(colSums(.)>250)%>%
  t()%>%
  clr(.+.1)%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  pivot_longer(-1,names_to = "tax_id",values_to = "s16_Reads")
WGS_phylum_CLR<-
  wgs_decontamd_phylum%>%
  filter(tax_id %in% all_bact_phyla)%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  select_if(colSums(.)>250)%>%
  t()%>%
  clr(.+.1)%>%
  as.data.frame()%>%
  rownames_to_column("Barcode")%>%
  pivot_longer(-1,names_to = "tax_id",values_to = "WGS_Reads")
RNA_16s_phylum_correlation<-
  rna_phylum_CLR%>%
  left_join(rna_annots)%>%
  left_join(s16_metadata,by=c("Sherlock_EAGLE_PID"="AdditionalAttributes","Type"="Tumor-NormalStatus"))%>%
  left_join(s16_phylum_CLR)%>%
  select(RNAseq_SampleID,SampleID,RNA_Reads,s16_Reads,tax_id)%>%
  filter(!s16_Reads==0,!RNA_Reads==0)%>%
  drop_na(s16_Reads)%>%
  left_join(kraken_taxonomy)%>%
  filter(name %in% c("Pseudomonadota","Actinomycetota","Bacillota","Bacteroidota"))%>%
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinoccocus-Thermus",
                        TRUE~name
  ))%>%
  ggplot(aes(x=RNA_Reads,y=s16_Reads))+
  geom_point(alpha=0.2)+
  facet_wrap(~name,scales = "free")+
  geom_smooth(method="lm")+
  labs(x="RNAseq CLR",y="16S CLR")+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black",size=3)
RNA_16s_phylum_correlation

RNA_WGS_phylum_correlation<-
  rna_phylum_CLR%>%
  left_join(rna_wgs_key,by=c("RNAseq_SampleID"="RNA"))%>%
  left_join(WGS_phylum_CLR,by=c("WGS"="Barcode","tax_id"))%>%
  select(RNAseq_SampleID,WGS,RNA_Reads,WGS_Reads,tax_id)%>%
  filter(!RNA_Reads==0,!WGS_Reads==0)%>%
  drop_na(WGS_Reads)%>%
  left_join(kraken_taxonomy)%>%
  filter(name %in% c("Pseudomonadota","Actinomycetota","Bacillota"))%>%
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinoccocus-Thermus",
                        TRUE~name
  ))%>%
  ggplot(aes(x=RNA_Reads,y=WGS_Reads))+
  geom_point(alpha=0.2)+
  facet_wrap(~name,scales = "free")+
  geom_smooth(method="lm")+
  labs(x="RNAseq CLR",y="WGS CLR")+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black",size=3)
RNA_WGS_phylum_correlation

WGS_16s_phylum_correlation<-
  WGS_phylum_CLR%>%
  filter(tax_id %in% all_bact_phyla)%>%
  left_join(wgs_16s_key,by=c("Barcode"="WGS"))%>%
  drop_na(`16s`)%>%
  left_join(s16_phylum_CLR,by=c("16s"="SampleID","tax_id"))%>%
  rename(SampleID="16s")%>%
  filter(!s16_Reads==0,!WGS_Reads==0)%>%
  drop_na(s16_Reads)%>%
  select(Barcode,SampleID,WGS_Reads,s16_Reads,tax_id)%>%
  left_join(kraken_taxonomy)%>%
  filter(name %in% c("Pseudomonadota","Actinomycetota","Bacillota"))%>%
  mutate(name=case_when(name=="Pseudomonadota"~"Proteobacteria",
                        name=="Actinomycetota"~"Actinobacteria",
                        name=="Bacillota"~"Firmicutes",
                        name=="Bacteroidota"~"Bacteroidetes",
                        name=="Deinococcota"~"Deinoccocus-Thermus",
                        TRUE~name
  ))%>%
  ggplot(aes(x=WGS_Reads,y=s16_Reads))+
  geom_point(alpha=0.2)+
  facet_wrap(~name,scales = "free")+
  geom_smooth(method="lm")+
  labs(x="WGS CLR",y="16S CLR")+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black",size=3)
WGS_16s_phylum_correlation

### ALPHA DIVERSITY
set.seed(123)
rna_overlap_diversity<-
  rna_combatd_decontamd%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  t()%>%
  as.data.frame()%>%
  filter(rowSums(.)>250)%>%
  calculate_avg_aDiversity(depth=250,index = "shannon")%>%
  as.data.frame()%>%
  rownames_to_column("RNAseq_SampleID")%>%
  rename(diversity_RNA=".")
s16_overlap_diversity<-
  s16_scrubbed%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  t()%>%
  as.data.frame()%>%
  filter(rowSums(.)>250)%>%
  calculate_avg_aDiversity(depth=250,index = "shannon")%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  rename(diversity_16s=".")
wgs_overlap_diversity<-
  wgs_decontamd%>%
  column_to_rownames("tax_id")%>%
  replace(is.na(.),0)%>%
  t()%>%
  as.data.frame()%>%
  filter(rowSums(.)>250)%>%
  calculate_avg_aDiversity(depth=250,index = "shannon")%>%
  as.data.frame()%>%
  rownames_to_column("SampleID")%>%
  rename(diversity_WGS=".")

rna_s16_diversity<-
  rna_overlap_diversity%>%
  left_join(rna_annots)%>%
  left_join(s16_metadata,by=c("Sherlock_EAGLE_PID"="AdditionalAttributes","Type"="Tumor-NormalStatus"))%>%
  left_join(s16_overlap_diversity)%>%
  drop_na(diversity_16s,diversity_RNA)%>%
  ggplot(aes(x=diversity_RNA,y=diversity_16s))+
  geom_point(alpha=0.2)+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black")+
  geom_smooth(method="lm")+
  labs(x="RNA diversity",y="16S diversity")
rna_s16_diversity

rna_wgs_diversity<-
  rna_overlap_diversity%>%
  left_join(rna_wgs_key,by=c("RNAseq_SampleID"="RNA"))%>%
  left_join(wgs_overlap_diversity,by=c("WGS"="SampleID"))%>%
  drop_na(diversity_RNA,diversity_WGS)%>%
  ggplot(aes(x=diversity_RNA,y=diversity_WGS))+
  geom_point(alpha=0.2)+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black")+
  geom_smooth(method="lm")+
  labs(x="RNA diversity",y="WGS diversity")
rna_wgs_diversity

wgs_s16_diversity<-
  wgs_overlap_diversity%>%
  left_join(wgs_16s_key,by=c("SampleID"="WGS"))%>%
  left_join(s16_overlap_diversity,by=c("16s"="SampleID"))%>%
  drop_na(diversity_WGS,diversity_16s)%>%
  ggplot(aes(x=diversity_WGS,y=diversity_16s))+
  geom_point(alpha=0.2)+
  stat_cor(method = "pearson",geom="shadowtext",bg.colour="white",color="black")+
  geom_smooth(method="lm")+
  labs(x="WGS diversity",y="16S diversity")
wgs_s16_diversity

### BETA DIVERSITY
samples_crossPlatform2<-
  samples_crossPlatform%>%
  separate(Sherlock_PID,into = c("Sherlock_PID","TNStatus"),sep="_")
all_samples_avgdist<-
  s16_scrubbed%>%
  left_join(rna_combatd_decontamd,"tax_id")%>%
  left_join(wgs_decontamd,"tax_id")%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>2)>0.05)%>%
  t()%>%
  avgdist(sample = 250,iterations=25)
still_paired_wgs16s<-
  samples_crossPlatform2%>%
  filter(platform %in% c("WGS","16s"),Barcode %in% labels(all_samples_avgdist))%>%
  filter(n()>1,
         .by=c(Sherlock_PID,TNStatus))%>%
  arrange(Sherlock_PID)
still_paired_rna16s<-
  samples_crossPlatform2%>%
  filter(platform %in% c("RNA","16s"),Barcode %in% labels(all_samples_avgdist))%>%
  filter(n()>1,
         .by=c(Sherlock_PID,TNStatus))%>%
  arrange(Sherlock_PID)
still_paired_rnaWGS<-
  samples_crossPlatform2%>%
  filter(platform %in% c("WGS","RNA"),Barcode %in% labels(all_samples_avgdist))%>%
  filter(n()>1,
         .by=c(Sherlock_PID,TNStatus))%>%
  arrange(Sherlock_PID)
all_three<-
  samples_crossPlatform2%>%
  filter(Barcode %in% labels(all_samples_avgdist))%>%
  filter(n()>2,
         .by=c(Sherlock_PID,TNStatus))%>%
  arrange(Sherlock_PID)%>%
  mutate(Patient_tissue=paste0(Sherlock_PID,"_",TNStatus))
ICC_beta_WGS_16s<-adonis2(formula=as.matrix(all_samples_avgdist)[still_paired_wgs16s$Barcode,
                                                                 still_paired_wgs16s$Barcode]~platform+Sherlock_PID,
                          data = still_paired_wgs16s,parallel = 4,by = "margin")
ICC_beta_RNA_WGS<-adonis2(formula=as.matrix(all_samples_avgdist)[still_paired_rnaWGS$Barcode,
                                                                 still_paired_rnaWGS$Barcode]~platform+Sherlock_PID,
                          data = still_paired_rnaWGS,parallel = 4,by = "margin")
ICC_beta_RNA_16s<-adonis2(formula=as.matrix(all_samples_avgdist)[still_paired_rna16s$Barcode,
                                                                 still_paired_rna16s$Barcode]~platform+TNStatus+Sherlock_PID,
                          data = still_paired_rna16s,parallel = 4,by = "margin")
ICC_beta_AllThree<-adonis2(formula=as.matrix(all_samples_avgdist)[all_three$Barcode,
                                                                  all_three$Barcode]~platform+Patient_tissue,
                           data = all_three,parallel = 4,by = "margin")

ICC_beta_results<-
  ICC_beta_WGS_16s%>%tidy()%>%
  filter(term=="Sherlock_PID")%>%
  mutate()%>%
  mutate(Experiment="WGS v 16S")%>%
  rbind(ICC_beta_RNA_16s%>%
          tidy()%>%
          filter(term=="Sherlock_PID")%>%
          mutate(Experiment="RNA v 16S"))%>%
  rbind(ICC_beta_RNA_WGS%>%
          tidy()%>%
          filter(term=="Sherlock_PID")%>%
          mutate(Experiment="RNA v WGS"))%>%
  rbind(ICC_beta_AllThree%>%
          tidy()%>%
          filter(term=="Patient_tissue")%>%
          mutate(Experiment="RNA v WGS v 16S"))%>%
  relocate(Experiment,.before = term)%>%
  select(Experiment,df,SumOfSqs,R2,p.value)%>%
  nice_table()%>%
  flextable::gen_grob()%>%
  as_ggplot()+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Mantel 16S-WGS correlation
mantel_16S1<-
  still_paired_wgs16s%>%
  filter(platform=="16s")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel_WGS1<-
  still_paired_wgs16s%>%
  filter(platform=="WGS")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel(xdis = as.matrix(all_samples_avgdist)[rownames(mantel_WGS1),
                                             rownames(mantel_WGS1)],
       ydis = as.matrix(all_samples_avgdist)[rownames(mantel_16S1),
                                             rownames(mantel_16S1)])
as.matrix(all_samples_avgdist)[rownames(mantel_WGS1),
                               rownames(mantel_WGS1)]%>%
  as.data.frame()%>%
  rownames_to_column("sample1")%>%
  pivot_longer(-1,names_to = "sample2",values_to = "WGS_dist")%>%
  cbind(as.matrix(all_samples_avgdist)[rownames(mantel_16S1),
                                       rownames(mantel_16S1)]%>%
          as.data.frame()%>%
          rownames_to_column("sample1")%>%
          pivot_longer(-1,names_to = "sample2",values_to = "s16_dist")%>%
          select(s16_dist)
        )%>%
  filter(WGS_dist!=0,s16_dist!=0)%>%
  ggplot(aes(x=-log(WGS_dist),y=-log(s16_dist)))+  
  geom_point(alpha=0.1)+
  geom_smooth(method="lm")+
  stat_cor()

# Mantel 16S-RNA correlation
mantel_16S2<-
  still_paired_rna16s%>%
  filter(platform=="16s")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel_RNA1<-
  still_paired_rna16s%>%
  filter(platform=="RNA")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel(xdis = as.matrix(all_samples_avgdist)[rownames(mantel_RNA1),
                                             rownames(mantel_RNA1)],
       ydis = as.matrix(all_samples_avgdist)[rownames(mantel_16S2),
                                             rownames(mantel_16S2)],method="spearman")

as.matrix(all_samples_avgdist)[rownames(mantel_RNA1),
                               rownames(mantel_RNA1)]%>%
  as.data.frame()%>%
  rownames_to_column("sample1")%>%
  pivot_longer(-1,names_to = "sample2",values_to = "RNA_dist")%>%
  cbind(as.matrix(all_samples_avgdist)[rownames(mantel_16S2),
                                       rownames(mantel_16S2)]%>%
          as.data.frame()%>%
          rownames_to_column("sample1")%>%
          pivot_longer(-1,names_to = "sample2",values_to = "s16_dist")%>%
          select(s16_dist)
  )%>%
  filter(RNA_dist!=0,s16_dist!=0)%>%
  ggplot(aes(x=-log(RNA_dist),y=-log(s16_dist)))+
  geom_point(alpha=0.1)+
  geom_smooth(method="lm")+
  stat_cor()

# Mantel WGS-RNA correlation
mantel_WGS2<-
  still_paired_rnaWGS%>%
  filter(platform=="WGS")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel_RNA2<-
  still_paired_rnaWGS%>%
  filter(platform=="RNA")%>%
  mutate(PID_TN=paste0(Sherlock_PID,"_",TNStatus))%>%
  column_to_rownames("Barcode")%>%
  arrange(PID_TN)
mantel(xdis = as.matrix(all_samples_avgdist)[rownames(mantel_RNA2),
                                             rownames(mantel_RNA2)],
       ydis = as.matrix(all_samples_avgdist)[rownames(mantel_WGS2),
                                             rownames(mantel_WGS2)],method="spearman")

as.matrix(all_samples_avgdist)[rownames(mantel_RNA2),
                               rownames(mantel_RNA2)]%>%
  as.data.frame()%>%
  rownames_to_column("sample1")%>%
  pivot_longer(-1,names_to = "sample2",values_to = "RNA_dist")%>%
  cbind(as.matrix(all_samples_avgdist)[rownames(mantel_WGS2),
                                       rownames(mantel_WGS2)]%>%
          as.data.frame()%>%
          rownames_to_column("sample1")%>%
          pivot_longer(-1,names_to = "sample2",values_to = "WGS_dist")%>%
          select(WGS_dist)
  )%>%
  filter(RNA_dist!=0,WGS_dist!=0)%>%
  ggplot(aes(x=-log(RNA_dist),y=-log(WGS_dist)))+
  geom_point(alpha=0.1)+
  geom_smooth(method="lm")+
  stat_cor()

