################################################
## ALDEX
aldex_in<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  select_if(colSums(.)>=500)%>%
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
  count(`Subject ID`)%>%nrow()
aldex_in<-
  aldex_in%>%
  select(any_of(aldex_pairedSample_filter))
aldex_out_rna<-
  ALDEx2::aldex(reads=aldex_in,
                conditions = (rna_annots%>%
                                filter(RNAseq_SampleID %in% aldex_pairedSample_filter)%>%
                                pull(Type,RNAseq_SampleID)))

fastancom_out<-
  fastANCOM::fastANCOM(t(aldex_in),
                       x=(rna_annots%>%
                            filter(RNAseq_SampleID %in% colnames(aldex_in))%>%
                            pull(Type,RNAseq_SampleID))
  )

differential_abundance_RNA<-
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
  labs(x=paste0("Difference between means\n(500 read cutoff, n=",aldex_N_count ," pairs)"),y="-log10(p.value)",color="")+
  geom_vline(xintercept = 0,lty=3)+
  # geom_texthline(yintercept=-log10(0.00315),color="blue",
  #                label="FDR=0.05",lty=2,family="Roboto Condensed",fontface="italic",
  #                hjust=.99,size=3)+
  geom_texthline(yintercept=-log10(0.05),color="red",
                 label="p=0.05",lty=2,family="Roboto Condensed",fontface="italic",
                 hjust=.99,size=3)+
  annotate("text",x = -Inf, y = -Inf, label="Normal enriched",hjust=-0.1,vjust=-1,
           family="Roboto Condensed",fontface="italic",size=3.25)+
  annotate("text",x = Inf, y = -Inf, label="Tumor enriched",hjust=1.1,vjust=-1,
           family="Roboto Condensed",fontface="italic",size=3.25)+
  xlim(-3,3)+
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
differential_abundance_RNA