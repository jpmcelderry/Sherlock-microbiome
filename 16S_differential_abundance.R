################################################
## ALDEX
aldex_in16s<-
  s16_scrubbed%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  select_if(colSums(.)>250)%>%
  filter(rowMeans(.>0)>0.01)
aldex_pairedSample_filter16s<-
  s16_metadata%>%
  filter(SampleID %in% colnames(aldex_in16s),`Tumor-NormalStatus` %in% c("Tumor","Normal"))%>%
  filter(nlevels(as.factor(`Tumor-NormalStatus`))>1,
         .by=AdditionalAttributes)%>%
  pull(SampleID)
aldex_in16s<-
  aldex_in16s%>%
  select(any_of(aldex_pairedSample_filter16s))
aldex_N_count16S<-
  s16_metadata%>%
  filter(SampleID %in% colnames(aldex_in16s))%>%
  filter(`Tumor-NormalStatus` %in% c("Tumor","Normal"),
         nlevels(as.factor(`Tumor-NormalStatus`))>1,
         .by=AdditionalAttributes)%>%
  count(AdditionalAttributes)%>%
  nrow()
aldex_out_16s<-
  ALDEx2::aldex(reads=aldex_in16s,
                conditions = (s16_metadata%>%
                                filter(SampleID %in% aldex_pairedSample_filter16s)%>%
                                pull(`Tumor-NormalStatus`,SampleID))
                )
aldex_out_16s%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%View()

differential_abundance_16S<-
  aldex_out_16s%>%
  rownames_to_column("tax_id")%>%
  left_join(unified_taxonomy)%>%
  mutate(color=case_when(wi.eBH>0.05~"FDR>0.5",
                         diff.btw>0~"tumor\nenriched",
                         diff.btw<0~"normal\nenriched"),
         experiment="16S")%>%
  ggplot(aes(x=diff.btw,
             y=-log10(wi.ep),
             label=if_else(wi.ep<=0.05,name,NA)))+
  geom_point(alpha=0.6)+
  geom_text_repel(show.legend = FALSE,size=2,min.segment.length = 0,force = 3)+
  labs(x=paste0("Difference between means\n(250 read cutoff, n=",aldex_N_count16S ," pairs)"),
       y="-log10(p.value)",color="")+
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
  ggh4x::facet_wrap2(~experiment,scales = "free_x",
                     strip = strip_themed(background_x = element_rect(fill=s16_color)))
differential_abundance_16S