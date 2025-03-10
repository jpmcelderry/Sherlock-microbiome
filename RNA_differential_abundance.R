################################################
## ALDEX
################################################
aldex_in<-
  rna_combatd_decontamd%>%
  filter(tax_id %in% all_bact_genera)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("tax_id")%>%
  select_if(colSums(.)>=500)
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

fastancom_out<-fastANCOM::fastANCOM(t(aldex_in),
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
  geom_vline(xintercept = 0,lty=3,alpha=0.5)+
  ggh4x::facet_wrap2(~experiment,scales = "free_x",strip = strip_themed(background_x = element_rect(fill=RNA_color)))
# co-associated bacteria
# sparcc_rna_combatd<-
#             SparCC.count(t(rna_combatd_decontamd%>%
#                    column_to_rownames("tax_id")%>%
#                    filter(rowMeans(.>5)>0.1)%>%
#                    filter(rownames(.) %in% all_bact_genera)))
# rownames(sparcc_rna_combatd$cor.w)<-
#             rownames(rna_combatd_decontamd%>%
#              column_to_rownames("tax_id")%>%
#              filter(rowMeans(.>5)>0.1)%>%
#              filter(rownames(.) %in% all_bact_genera))
# colnames(sparcc_rna_combatd$cor.w)<-
#             rownames(rna_combatd_decontamd%>%
#              column_to_rownames("tax_id")%>%
#              filter(rowMeans(.>5)>0.1)%>%
#              filter(rownames(.) %in% all_bact_genera))
# diag(sparcc_rna_combatd$cor.w)<-NA
# paired_rnas<-
#   rna_annots%>%
#   filter(RNAseq_SampleID %in% names(rna.adj.vars))%>%
#   select(`Subject ID`,Type,RNAseq_SampleID)%>%
#   pivot_wider(names_from = "Type",values_from = "RNAseq_SampleID")%>%
#   drop_na()%>%
#   pivot_longer(c("Tumor","Normal"),names_to = "Type",values_to = "ID")%>%
#   unlist()
# Tumor_Normal_tile<-
#   rna_combatd_decontamd%>%
#   column_to_rownames("tax_id")%>%
#   filter(rownames(.) %in% all_bact_genera)%>%
#   filter(rowMeans(.>0)>0.05)%>%
#   t()%>%
#   `+`(0.5)%>%
#   clr()%>%
#   as.data.frame()%>%
#   rownames_to_column("RNAseq_SampleID")%>%
#   pivot_longer(-1,names_to="tax_id",values_to = "counts")%>%
#   left_join(kraken_microbeTaxonomy)%>%
#   left_join(rna_annots)%>%
#   group_by(Sherlock_PID)%>%
#   filter(nlevels(as.factor(Type))>1)%>%
#   ungroup()%>%
#   group_by(name,Type)%>%
#   summarise(mean=mean(counts))%>%
#   pivot_wider(names_from = "Type",values_from = "mean")%>%
#   summarise(Difference=Tumor-Normal)
# 
# genera_correlation_matrix<-
#   as.data.frame(sparcc_rna_combatd$cor.w)%>%
#   .[apply(.,MARGIN = 1,function(x){max(abs(x),na.rm = T)})>0.2,
#     apply(.,MARGIN = 2,function(x){max(abs(x),na.rm = T)})>0.2]%>%
#   replace(is.na(.),1)%>%
#   rownames_to_column("tax_id1")%>%
#   pivot_longer(-1,names_to = "tax_id2",values_to="corr")%>%
#   drop_na()%>%
#   left_join(kraken_microbeTaxonomy[,c("tax_id","name")]%>%rename(name1="name"),by=c("tax_id1"="tax_id"))%>%
#   left_join(kraken_microbeTaxonomy[,c("tax_id","name")]%>%rename(name2="name"),by=c("tax_id2"="tax_id"))%>%
#   left_join(Tumor_Normal_tile,by=c("name2"="name"))%>%
#   tidyHeatmap::heatmap(.row = name1,.column = name2,.value = corr,
#                        palette_value = circlize::colorRamp2(
#                          seq(-1, 1, length.out = 11),
#                          RColorBrewer::brewer.pal(11, "RdBu")),
#                        row_names_side=c("left"),
#                        row_names_gp = grid::gpar(fontsize = 6.5),
#                        column_names_gp = grid::gpar(fontsize = 7),
#                        column_names_rot=45,
#                        column_title=NULL,row_title=NULL,
#                        clustering_method_columns="ward.D2",
#                        clustering_method_rows="ward.D2",column_km=3,show_row_dend=F,show_column_dend=F,
#                        column_km_repeats=2000,km=3,row_km_repeats=2000,
#                        heatmap_legend_param=list(title="Correlation\nCoefficient"))%>%
#   tidyHeatmap::add_tile(.column=Difference,
#                         annotation_label="Abundance Diff.",
#                         annotation_name_side="right",
#                         annotation_legend_param=list(title="Abundance Difference\n(Tumor minus Normal)"),
#                         palette = circlize::colorRamp2(
#                           seq(1, -1, length.out = 11),
#                           RColorBrewer::brewer.pal(11, "BrBG")))
