new_Gene_Centric_Coding_mac <- function (chr, gene_name, category = c("all_categories", "plof", 
                                       "plof_ds", "missense", "disruptive_missense", "synonymous"), 
          genofile, obj_nullmodel, rare_maf_cutoff = 0.01, rv_num_cutoff = 2, 
          QC_label = "annotation/filter", variant_type = c("SNV", "Indel", 
                                                           "variant"), geno_missing_imputation = c("mean", "minor"), 
          Annotation_dir = "annotation/info/FunctionalAnnotation", 
          Annotation_name_catalog, Use_annotation_weights = c(TRUE, 
                                                              FALSE), Annotation_name = NULL, silent = FALSE, MAC_super_rare=MAC_super_rare, cMAC_super_rare = cMAC_super_rare) 
{
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  genes <- genes_info[genes_info[, 2] == chr, ]
  if (category == "all_categories") {
    results <- new_coding_mac(chr, gene_name, genofile, obj_nullmodel, 
                      genes, rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                      QC_label = QC_label, variant_type = variant_type, 
                      geno_missing_imputation = geno_missing_imputation, 
                      Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                      Use_annotation_weights = Use_annotation_weights, 
                      Annotation_name = Annotation_name, silent = silent, MAC_super_rare=MAC_super_rare, cMAC_super_rare = cMAC_super_rare)
  }
  if (category == "plof") {
    results <- plof(chr, gene_name, genofile, obj_nullmodel, 
                    genes, rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                    QC_label = QC_label, variant_type = variant_type, 
                    geno_missing_imputation = geno_missing_imputation, 
                    Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                    Use_annotation_weights = Use_annotation_weights, 
                    Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "plof_ds") {
    results <- plof_ds(chr, gene_name, genofile, obj_nullmodel, 
                       genes, rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                       QC_label = QC_label, variant_type = variant_type, 
                       geno_missing_imputation = geno_missing_imputation, 
                       Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                       Use_annotation_weights = Use_annotation_weights, 
                       Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "missense") {
    results <- missense(chr, gene_name, genofile, obj_nullmodel, 
                        genes, rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                        QC_label = QC_label, variant_type = variant_type, 
                        geno_missing_imputation = geno_missing_imputation, 
                        Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                        Use_annotation_weights = Use_annotation_weights, 
                        Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "disruptive_missense") {
    results <- disruptive_missense(chr, gene_name, genofile, 
                                   obj_nullmodel, genes, rare_maf_cutoff = rare_maf_cutoff, 
                                   rv_num_cutoff = rv_num_cutoff, QC_label = QC_label, 
                                   variant_type = variant_type, geno_missing_imputation = geno_missing_imputation, 
                                   Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                                   Use_annotation_weights = Use_annotation_weights, 
                                   Annotation_name = Annotation_name, silent = silent)
  }
  if (category == "synonymous") {
    results <- synonymous(chr, gene_name, genofile, obj_nullmodel, 
                          genes, rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, 
                          QC_label = QC_label, variant_type = variant_type, 
                          geno_missing_imputation = geno_missing_imputation, 
                          Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog, 
                          Use_annotation_weights = Use_annotation_weights, 
                          Annotation_name = Annotation_name, silent = silent)
  }
  return(results)
}