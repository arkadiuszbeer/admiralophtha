oe <- rice::rice_read("root/clinical_studies/RO7200220/CDT70193/BP40899/data_analysis/BASE/data/sdtmv/oe.sas7bdat") %>%
  filter(OETESTCD %in% c("CPTTHM", "CSUBTH", "ASCAN", "BSCAN")) %>%
  filter(!is.na(OELAT))

adsl_sas <- rice::rice_read("root/clinical_studies/RO7200220/CDT70193/BP40899/data_analysis/BASE/qa/outdata_vad/adsl.sas7bdat")
adoe_sas <- rice::rice_read("root/clinical_studies/RO7200220/CDT70193/BP40899/data_analysis/BASE/qa/outdata_vad/adoe.sas7bdat")
