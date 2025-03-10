aldex_in16s<-
  s16_scrubbed%>%
  replace(is.na(.),0)%>%
  filter(tax_id %in% all_bact_genera)%>%
  column_to_rownames("tax_id")%>%
  filter(rowMeans(.>0)>0.01)%>%
  select_if(colSums(.)>300)
aldex_pairedSample_filter16s<-
  s16_metadata%>%
  filter(SampleID %in% colnames(aldex_in16s),`Tumor-NormalStatus` %in% c("Tumor","Normal"))%>%
  filter(nlevels(as.factor(`Tumor-NormalStatus`))>1,
         .by=AdditionalAttributes)%>%
  pull(SampleID)
aldex_in16s<-
  aldex_in16s%>%
  select(any_of(aldex_pairedSample_filter16s))
aldex_out_16s<-ALDEx2::aldex(reads=aldex_in16s,
                             conditions = (s16_metadata%>%
                                             filter(SampleID %in% aldex_pairedSample_filter16s)%>%
                                             pull(`Tumor-NormalStatus`,SampleID)))
aldex_out_16s%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%View()

differential_abundance16s<-aldex_out_16s%>%
  rownames_to_column("tax_id")%>%
  left_join(kraken_taxonomy)%>%
  ggplot(aes(x=rab.win.Tumor-rab.win.Normal,y=-log10(wi.eBH),
             color=wi.eBH<0.05,label=if_else(wi.eBH<0.05,name,"")))+
  geom_point()+
  geom_text_repel(show.legend = FALSE,size=3)+
  labs(x="mean CLR difference (Tumor - Normal)",y="-log10(FDR)",color="FDR<0.05")+
  xlim(-3,3)+
  geom_vline(xintercept = 0,lty=3)