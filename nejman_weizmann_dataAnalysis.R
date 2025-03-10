nejman_lung<-read_delim("~/nejman_lung_data.txt")
nejman_lung_annos<-read_delim("~/nejman_lung_metadata.txt")

nejman_lung<-nejman_lung%>%
  mutate(`Prevalence in DNA extraction/NTC negative controls`=str_replace(`Prevalence in DNA extraction/NTC negative controls`,pattern = "%",replacement = ""),
         `Prevalence in Paraf. Controls`=str_replace(`Prevalence in Paraf. Controls`,pattern = "%",replacement = ""))%>%
  mutate_at(c("Prevalence in DNA extraction/NTC negative controls","Prevalence in Paraf. Controls"),as.numeric)%>%
  rowwise()%>%
  mutate(taxonomy=paste(domain,phylum,class,order,family,genus,sep = "_"))%>%
  ungroup()%>%
  filter(`Prevalence in DNA extraction/NTC negative controls`<7.5,`Prevalence in Paraf. Controls`<7.5)%>%
  select(taxonomy,starts_with("1"))%>%
  # filter(!grepl("Unknown",genus))%>%
  #do(aggregate(.~genus,data=.,FUN=sum))
  do(aggregate(.~taxonomy,data=.,FUN=sum))
nejman_lung%>%
  .[,-1]%>%
  apply(MARGIN = 2,function(x){x/sum(x)})%>%
  rowMeans()%>%
  bind_cols(nejman_lung$taxonomy,.)%>%
  arrange(-`...2`)
nejman_avgdist<-nejman_lung%>%column_to_rownames("genus")%>%
  t()%>%
  round()%>%
  avgdist(1000)



fast.adonis(formula=as.matrix(nejman_avgdist)~Material+Center+`DNA Extraction batch`+`PCR batch`+`Tissue type`,
            data=(nejman_lung_annos%>%column_to_rownames("Sample_ID (WIS)")%>%.[labels(nejman_avgdist),]),
            by="terms",boot.times = 100,parallel = 4)$aov.tab
