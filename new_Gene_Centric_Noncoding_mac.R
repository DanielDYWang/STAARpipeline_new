new_Gene_Centric_Noncoding_mac<- function (chr, gene_name, category = c("all_categories", "downstream", 
                                       "upstream", "UTR", "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", 
                                       "enhancer_DHS"), genofile, obj_nullmodel, rare_maf_cutoff = 0.01, 
          rv_num_cutoff = 2, QC_label = "annotation/filter", variant_type = c("SNV", 
                                                                              "Indel", "variant"), geno_missing_imputation = c("mean", 
                                                                                                                               "minor"), Annotation_dir = "annotation/info/FunctionalAnnotation", 
          Annotation_name_catalog, Use_annotation_weights = c(TRUE, 
                                                              FALSE), Annotation_name = NULL, silent = FALSE, MAC_super_rare=MAC_super_rare, cMAC_super_rare = cMAC_super_rare) 
{
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  if (category == "all_categories") {
    results <- new_noncoding_mac(chr, gene_name, genofile, obj_nullmodel, 
                         rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                         QC_label = QC_label, variant_type = variant_type, 
                         geno_missing_imputation = geno_missing_imputation, 
                         Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                         Use_annotation_weights = Use_annotation_weights, 
                         Annotation_name = Annotation_name, silent = silent, MAC_super_rare=MAC_super_rare, cMAC_super_rare = cMAC_super_rare)
  }
  if (category == "downstream") {
    results <- downstream(chr, gene_name, genofile, obj_nullmodel, 
                          rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                          QC_label = QC_label, variant_type = variant_type, 
                          geno_missing_imputation = geno_missing_imputation, 
                          Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                          Use_annotation_weights = Use_annotation_weights, 
                          Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "upstream") {
    results <- upstream(chr, gene_name, genofile, obj_nullmodel, 
                        rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                        QC_label = QC_label, variant_type = variant_type, 
                        geno_missing_imputation = geno_missing_imputation, 
                        Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                        Use_annotation_weights = Use_annotation_weights, 
                        Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "UTR") {
    results <- UTR(chr, gene_name, genofile, obj_nullmodel, 
                   rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                   QC_label = QC_label, variant_type = variant_type, 
                   geno_missing_imputation = geno_missing_imputation, 
                   Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                   Use_annotation_weights = Use_annotation_weights, 
                   Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "promoter_CAGE") {
    results <- promoter_CAGE(chr, gene_name, genofile, obj_nullmodel, 
                             rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                             QC_label = QC_label, variant_type = variant_type, 
                             geno_missing_imputation = geno_missing_imputation, 
                             Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                             Use_annotation_weights = Use_annotation_weights, 
                             Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "promoter_DHS") {
    results <- promoter_DHS(chr, gene_name, genofile, obj_nullmodel, 
                            rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                            QC_label = QC_label, variant_type = variant_type, 
                            geno_missing_imputation = geno_missing_imputation, 
                            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                            Use_annotation_weights = Use_annotation_weights, 
                            Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "enhancer_CAGE") {
    results <- enhancer_CAGE(chr, gene_name, genofile, obj_nullmodel, 
                             rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                             QC_label = QC_label, variant_type = variant_type, 
                             geno_missing_imputation = geno_missing_imputation, 
                             Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                             Use_annotation_weights = Use_annotation_weights, 
                             Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "enhancer_DHS") {
    results <- enhancer_DHS(chr, gene_name, genofile, obj_nullmodel, 
                            rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                            QC_label = QC_label, variant_type = variant_type, 
                            geno_missing_imputation = geno_missing_imputation, 
                            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                            Use_annotation_weights = Use_annotation_weights, 
                            Annotation_name = Annotation_name, silent = silent)
  }
  return(results)
}
