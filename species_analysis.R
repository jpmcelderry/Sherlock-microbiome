relabund_byExperiment3<-
  rna%>% # Join data
  mutate(experiment="RNA-seq",.after="tax_id")%>%
  full_join(read_delim("/Volumes/Sherlock_Lung/JohnMce/Kraken_WGS3/kraken0.1-collated.txt")%>%
              rename(tax_id=`4`)%>%
              mutate(tax_id=as.character(tax_id))%>%
              mutate(experiment="WGS",.after="tax_id"))%>%
  right_join(unified_taxonomy,.)%>%
  filter(type=="species",str_detect(taxonomy,"Bacteria"))%>%
  mutate(genus=str_split(name," ")%>%sapply(`[`,1),.after = name)%>%
  filter(!genus %in% c(salter_list_nonHuman$Genus,"Cutibacterium"))%>%
  pivot_longer(-c(1:6),names_to = "Sample",values_to = "Reads")%>% # Pivot data, minus taxonomy info
  replace_na(list(Reads=0))%>%
  mutate(Relabund=Reads/sum(Reads),.by=Sample)%>% # Transform to relative abundance
  rename(Barcode="Sample")%>%filter(Barcode %in% samples_crossPlatform$Barcode)%>% # Filter to only cross-platform samples
  left_join(all_sample_TNstatus[,c("Barcode","Source")])%>% # Join Tumor-Normal info
  drop_na()%>%
  filter(!(Source=="Blood"))%>%
  summarise(Relabund=sum(Relabund)/n(),.by=c(experiment,Source,name,tax_id))%>% # Determine Platform/Tissue Type mean relAbund
  mutate(rank=order(order(Relabund, decreasing=TRUE)),.by=c(experiment,Source))%>% # Determine relAbund ranks per sample group
  mutate(min_rank=min(rank),.by=name) # Calculate the lowest rank per bacteria among all platforms/tissues 
species_lvl<-
  relabund_byExperiment3%>%
  mutate(name=if_else(min_rank>25,"Other",name))%>% # Lump any bacteria that's never top 15 abundant
  summarise(Relabund=sum(Relabund),
            .by=c(name,Source,experiment))%>% # Sum the "Other"s into 1 datapoint 
  mutate(name=fct_relevel(name,"Other",after=Inf))%>%# Place Other as last category
  ggplot(aes(x=Source,y=Relabund,fill=name))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual("Genus",values=unname(SBScolor))+
  labs(y="Mean Relative Abundance",x="n=798")+
  guides(fill=guide_legend(ncol=2))+
  facet_grid2(~experiment, space="free_x", scales="free_x",
              strip = strip_themed(background_x = list(element_rect(fill=s16_color),
                                                       element_rect(fill=RNA_color),
                                                       element_rect(fill=WGS_color))))+
  theme(strip.placement = "outside",
        panel.spacing=unit(0,"cm"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=45))
species_lvl

aldex_species_list<-
  unified_taxonomy%>%
  filter(type=="species")%>%
  filter(str_detect(name,"Pseudomonas")|str_detect(name,"Corynebacterium")|str_detect(name,"Staphylococcus")|str_detect(name,"Acinetobacter")|str_detect(name,"Streptococcus"))%>%
  pull(tax_id)
aldex_in<-
  rna%>%
  select(any_of(colnames(rna_combatd_decontamd)))%>%
  filter(tax_id %in% aldex_species_list)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  select_if(colSums(.)>=50)%>%
  filter(rowMeans(.>0)>0.01)
aldex_pairedSample_filter<-
  rna_annots%>%
  filter(RNAseq_SampleID %in% colnames(aldex_in))%>%
  filter(nlevels(as.factor(Type))>1,
         .by=`Subject ID`)%>%
  pull(RNAseq_SampleID)
aldex_N_count<-
  rna_annots%>%
  filter(RNAseq_SampleID %in% colnames(aldex_in))%>%
  filter(nlevels(as.factor(Type))>1,
         .by=`Subject ID`)%>%
  count(`Subject ID`)%>%
  nrow()
aldex_in<-
  aldex_in%>%
  select(any_of(aldex_pairedSample_filter))
aldex_out_rna<-
  ALDEx2::aldex(reads=aldex_in,
                conditions = (rna_annots%>%
                                filter(RNAseq_SampleID %in% aldex_pairedSample_filter)%>%
                                pull(Type,RNAseq_SampleID)))
aldex_species<-
  aldex_out_rna%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%
  mutate(color=case_when(wi.eBH>0.05~"FDR>0.5",
                         diff.btw>0~"tumor\nenriched",
                         diff.btw<0~"normal\nenriched"),
         experiment="RNA-seq")%>%
  ggplot(aes(x=diff.btw,y=-log10(wi.ep),
             label=if_else(wi.ep<=0.05,name,NA)))+
  geom_point(alpha=0.6)+
  geom_text_repel(show.legend = FALSE,size=2,min.segment.length = 0,force = 3)+
  labs(x=paste0("Difference between means\nn=",aldex_N_count ," pairs"),y="-log10(p.value)",color="")+
  geom_vline(xintercept = 0,lty=3)+
  geom_texthline(yintercept=-log10(0.00315),color="blue",
                 label="FDR=0.05",lty=2,family="Roboto Condensed",fontface="italic",
                 hjust=.99,size=3)+
  geom_texthline(yintercept=-log10(0.05),color="red",
                 label="p=0.05",lty=2,family="Roboto Condensed",fontface="italic",
                 hjust=.99,size=3)+
  annotate("text",x = -Inf, y = -Inf, label="Normal enriched",hjust=-0.1,vjust=-1,
           family="Roboto Condensed",fontface="italic",size=3.25)+
  annotate("text",x = Inf, y = -Inf, label="Tumor enriched",hjust=1.1,vjust=-1,
           family="Roboto Condensed",fontface="italic",size=3.25)+
  xlim(-3,3)+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)
  #ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))

