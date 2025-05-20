new_Gene_Centric_Noncoding_Results_Summary <- function (agds_dir, gene_centric_noncoding_jobs_num, input_path, 
          output_path, gene_centric_results_name, ncRNA_jobs_num, ncRNA_input_path, 
          ncRNA_output_path, ncRNA_results_name, obj_nullmodel, known_loci = NULL, 
          method_cond = c("optimal", "naive"), QC_label = "annotation/filter", 
          geno_missing_imputation = c("mean", "minor"), variant_type = c("SNV", 
                                                                         "Indel", "variant"), Annotation_dir = "annotation/info/FunctionalAnnotation", 
          Annotation_name_catalog, Use_annotation_weights = FALSE, 
          Annotation_name = NULL, alpha = 2.5e-06, suggest.alpha = NA, manhattan_plot = FALSE, 
          QQ_plot = FALSE, cMAC = 10, gene_name = TRUE, min_y = NA) 
{
  method_cond <- match.arg(method_cond)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  results_noncoding_genome <- c()
  for (kk in 1:gene_centric_noncoding_jobs_num) {
    print(kk)
    results_noncoding <- get(load(paste0(input_path, gene_centric_results_name, 
                                         "_", kk, ".Rdata")))
    results_noncoding_genome <- c(results_noncoding_genome, 
                                  results_noncoding)
  }
  results_UTR_genome <- c()
  results_upstream_genome <- c()
  results_downstream_genome <- c()
  results_promoter_CAGE_genome <- c()
  results_promoter_DHS_genome <- c()
  results_enhancer_CAGE_genome <- c()
  results_enhancer_DHS_genome <- c()
  for (kk in 1:length(results_noncoding_genome)) {
    results <- results_noncoding_genome[[kk]]
    if (is.null(results) == FALSE) {
      if (results[3] == "UTR" & results[5] >= cMAC) {
        results_UTR_genome <- rbind(results_UTR_genome, 
                                    results)
      }
      if (results[3] == "upstream" & results[5] >= cMAC) {
        results_upstream_genome <- rbind(results_upstream_genome, 
                                         results)
      }
      if (results[3] == "downstream" & results[5] >= cMAC) {
        results_downstream_genome <- rbind(results_downstream_genome, 
                                           results)
      }
      if (results[3] == "promoter_CAGE" & results[5] >= cMAC) {
        results_promoter_CAGE_genome <- rbind(results_promoter_CAGE_genome, 
                                              results)
      }
      if (results[3] == "promoter_DHS" & results[5] >= cMAC) {
        results_promoter_DHS_genome <- rbind(results_promoter_DHS_genome, 
                                             results)
      }
      if (results[3] == "enhancer_CAGE" & results[5] >= cMAC) {
        results_enhancer_CAGE_genome <- rbind(results_enhancer_CAGE_genome, 
                                              results)
      }
      if (results[3] == "enhancer_DHS" & results[5] >= cMAC) {
        results_enhancer_DHS_genome <- rbind(results_enhancer_DHS_genome, 
                                             results)
      }
    }
    if (kk%%1000 == 0) {
      print(kk)
    }
  }
  save(results_UTR_genome, file = paste0(output_path, "UTR.Rdata"))
  save(results_upstream_genome, file = paste0(output_path, 
                                              "upstream.Rdata"))
  save(results_downstream_genome, file = paste0(output_path, 
                                                "downstream.Rdata"))
  save(results_promoter_CAGE_genome, file = paste0(output_path, 
                                                   "promoter_CAGE.Rdata"))
  save(results_promoter_DHS_genome, file = paste0(output_path, 
                                                  "promoter_DHS.Rdata"))
  save(results_enhancer_CAGE_genome, file = paste0(output_path, 
                                                   "enhancer_CAGE.Rdata"))
  save(results_enhancer_DHS_genome, file = paste0(output_path, 
                                                  "enhancer_DHS.Rdata"))
  results_ncRNA_genome <- c()
  for (kk in 1:ncRNA_jobs_num) {
    print(kk)
    results_ncRNA <- get(load(paste0(ncRNA_input_path, ncRNA_results_name, 
                                     "_", kk, ".Rdata")))
    results_ncRNA_genome <- rbind(results_ncRNA_genome, results_ncRNA) 
  }
  
  results_ncRNA_genome <- data.frame(results_ncRNA_genome)
  results_ncRNA_genome$cMAC <- as.numeric( results_ncRNA_genome$cMAC)
  results_ncRNA_genome_2 <-c()
  for (kk in 1:dim(results_ncRNA_genome)[1]) {
    if (results_ncRNA_genome$cMAC[[kk]] > cMAC){
    results_ncRNA_genome_2 <- rbind(results_ncRNA_genome_2, results_ncRNA_genome[kk,])
    }
  }
  
  save(results_ncRNA_genome_2, file = paste0(ncRNA_output_path, 
                                           "results_ncRNA_genome.Rdata"))
  UTR_sig <- results_UTR_genome[results_UTR_genome[, "STAAR-O"] < 
                                  alpha, , drop = FALSE]
  write.csv(UTR_sig, file = paste0(output_path, "UTR_sig.csv"))
  noncoding_sig <- c()
  noncoding_sig <- rbind(noncoding_sig, UTR_sig)
  upstream_sig <- results_upstream_genome[results_upstream_genome[, 
                                                                  "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(upstream_sig, file = paste0(output_path, "upstream_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, upstream_sig)
  downstream_sig <- results_downstream_genome[results_downstream_genome[, 
                                                                        "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(downstream_sig, file = paste0(output_path, "downstream_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, downstream_sig)
  promoter_CAGE_sig <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[, 
                                                                                 "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(promoter_CAGE_sig, file = paste0(output_path, "promoter_CAGE_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, promoter_CAGE_sig)
  promoter_DHS_sig <- results_promoter_DHS_genome[results_promoter_DHS_genome[, 
                                                                              "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(promoter_DHS_sig, file = paste0(output_path, "promoter_DHS_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, promoter_DHS_sig)
  enhancer_CAGE_sig <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[, 
                                                                                 "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(enhancer_CAGE_sig, file = paste0(output_path, "enhancer_CAGE_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, enhancer_CAGE_sig)
  enhancer_DHS_sig <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[, 
                                                                              "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(enhancer_DHS_sig, file = paste0(output_path, "enhancer_DHS_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, enhancer_DHS_sig)
  ncRNA_sig <- results_ncRNA_genome_2[as.numeric(results_ncRNA_genome_2[, 
                                                         "STAAR.O"]) < alpha, , drop = FALSE]
  ncRNA_sig <- lapply(ncRNA_sig,as.character)
  write.csv(ncRNA_sig, file = paste0(ncRNA_output_path, "ncRNA_sig.csv"))
  noncoding_sig <- rbind(noncoding_sig, ncRNA_sig)
  write.csv(noncoding_sig, file = paste0(output_path, "noncoding_sig.csv"))
  if (length(known_loci) != 0) {
    UTR_sig_cond <- c()
    if (length(UTR_sig) != 0) {
      if ((class(UTR_sig)[1] != "matrix") & (class(UTR_sig)[1] != 
                                             "data.frame")) {
        UTR_sig <- matrix(UTR_sig, nrow = 1)
      }
      for (k in 1:dim(UTR_sig)[1]) {
        chr <- as.numeric(UTR_sig[k, 2])
        gene_name <- as.character(UTR_sig[k, 1])
        category <- as.character(UTR_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        UTR_sig_cond <- rbind(UTR_sig_cond, res_cond)
        seqClose(genofile)
      }
    }
    write.csv(UTR_sig_cond, file = paste0(output_path, "UTR_sig_cond.csv"))
    noncoding_sig_cond <- c()
    noncoding_sig_cond <- rbind(noncoding_sig_cond, UTR_sig_cond)
    upstream_sig_cond <- c()
    if (length(upstream_sig) != 0) {
      if ((class(upstream_sig)[1] != "matrix") & (class(upstream_sig)[1] != 
                                                  "data.frame")) {
        upstream_sig <- matrix(upstream_sig, nrow = 1)
      }
      for (k in 1:dim(upstream_sig)[1]) {
        chr <- as.numeric(upstream_sig[k, 2])
        gene_name <- as.character(upstream_sig[k, 1])
        category <- as.character(upstream_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        upstream_sig_cond <- rbind(upstream_sig_cond, 
                                   res_cond)
        seqClose(genofile)
      }
    }
    write.csv(upstream_sig_cond, file = paste0(output_path, 
                                               "upstream_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, upstream_sig_cond)
    downstream_sig_cond <- c()
    if (length(downstream_sig) != 0) {
      if ((class(downstream_sig)[1] != "matrix") & (class(downstream_sig)[1] != 
                                                    "data.frame")) {
        downstream_sig <- matrix(downstream_sig, nrow = 1)
      }
      for (k in 1:dim(downstream_sig)[1]) {
        chr <- as.numeric(downstream_sig[k, 2])
        gene_name <- as.character(downstream_sig[k, 1])
        category <- as.character(downstream_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        downstream_sig_cond <- rbind(downstream_sig_cond, 
                                     res_cond)
        seqClose(genofile)
      }
    }
    write.csv(downstream_sig_cond, file = paste0(output_path, 
                                                 "downstream_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, downstream_sig_cond)
    promoter_CAGE_sig_cond <- c()
    if (length(promoter_CAGE_sig) != 0) {
      if ((class(promoter_CAGE_sig)[1] != "matrix") & (class(promoter_CAGE_sig)[1] != 
                                                       "data.frame")) {
        promoter_CAGE_sig <- matrix(promoter_CAGE_sig, 
                                    nrow = 1)
      }
      for (k in 1:dim(promoter_CAGE_sig)[1]) {
        chr <- as.numeric(promoter_CAGE_sig[k, 2])
        gene_name <- as.character(promoter_CAGE_sig[k, 
                                                    1])
        category <- as.character(promoter_CAGE_sig[k, 
                                                   3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        promoter_CAGE_sig_cond <- rbind(promoter_CAGE_sig_cond, 
                                        res_cond)
        seqClose(genofile)
      }
    }
    write.csv(promoter_CAGE_sig_cond, file = paste0(output_path, 
                                                    "promoter_CAGE_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, promoter_CAGE_sig_cond)
    promoter_DHS_sig_cond <- c()
    if (length(promoter_DHS_sig) != 0) {
      if ((class(promoter_DHS_sig)[1] != "matrix") & (class(promoter_DHS_sig)[1] != 
                                                      "data.frame")) {
        promoter_DHS_sig <- matrix(promoter_DHS_sig, 
                                   nrow = 1)
      }
      for (k in 1:dim(promoter_DHS_sig)[1]) {
        chr <- as.numeric(promoter_DHS_sig[k, 2])
        gene_name <- as.character(promoter_DHS_sig[k, 
                                                   1])
        category <- as.character(promoter_DHS_sig[k, 
                                                  3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        promoter_DHS_sig_cond <- rbind(promoter_DHS_sig_cond, 
                                       res_cond)
        seqClose(genofile)
      }
    }
    write.csv(promoter_DHS_sig_cond, file = paste0(output_path, 
                                                   "promoter_DHS_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, promoter_DHS_sig_cond)
    enhancer_CAGE_sig_cond <- c()
    if (length(enhancer_CAGE_sig) != 0) {
      if ((class(enhancer_CAGE_sig)[1] != "matrix") & (class(enhancer_CAGE_sig)[1] != 
                                                       "data.frame")) {
        enhancer_CAGE_sig <- matrix(enhancer_CAGE_sig, 
                                    nrow = 1)
      }
      for (k in 1:dim(enhancer_CAGE_sig)[1]) {
        chr <- as.numeric(enhancer_CAGE_sig[k, 2])
        gene_name <- as.character(enhancer_CAGE_sig[k, 
                                                    1])
        category <- as.character(enhancer_CAGE_sig[k, 
                                                   3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        enhancer_CAGE_sig_cond <- rbind(enhancer_CAGE_sig_cond, 
                                        res_cond)
        seqClose(genofile)
      }
    }
    write.csv(enhancer_CAGE_sig_cond, file = paste0(output_path, 
                                                    "enhancer_CAGE_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, enhancer_CAGE_sig_cond)
    enhancer_DHS_sig_cond <- c()
    if (length(enhancer_DHS_sig) != 0) {
      if ((class(enhancer_DHS_sig)[1] != "matrix") & (class(enhancer_DHS_sig)[1] != 
                                                      "data.frame")) {
        enhancer_DHS_sig <- matrix(enhancer_DHS_sig, 
                                   nrow = 1)
      }
      for (k in 1:dim(enhancer_DHS_sig)[1]) {
        chr <- as.numeric(enhancer_DHS_sig[k, 2])
        gene_name <- as.character(enhancer_DHS_sig[k, 
                                                   1])
        category <- as.character(enhancer_DHS_sig[k, 
                                                  3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Noncoding_cond(chr = chr, 
                                                gene_name = gene_name, category = category, 
                                                genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                                known_loci = known_loci, variant_type = variant_type, 
                                                geno_missing_imputation = geno_missing_imputation, 
                                                QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                                Use_annotation_weights = Use_annotation_weights, 
                                                Annotation_name_catalog = Annotation_name_catalog, 
                                                Annotation_name = Annotation_name)
        enhancer_DHS_sig_cond <- rbind(enhancer_DHS_sig_cond, 
                                       res_cond)
        seqClose(genofile)
      }
    }
    write.csv(enhancer_DHS_sig_cond, file = paste0(output_path, 
                                                   "enhancer_DHS_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, enhancer_DHS_sig_cond)
    ncRNA_sig_cond <- c()
    if (length(ncRNA_sig) != 0) {
      if ((class(ncRNA_sig)[1] != "matrix") & (class(ncRNA_sig)[1] != 
                                               "data.frame")) {
        ncRNA_sig <- matrix(ncRNA_sig, nrow = 1)
      }
      for (k in 1:dim(ncRNA_sig)[1]) {
        chr <- as.numeric(ncRNA_sig[k, 2])
        gene_name <- as.character(ncRNA_sig[k, 1])
        category <- as.character(ncRNA_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- ncRNA_cond(chr = chr, gene_name = gene_name, 
                               genofile = genofile, obj_nullmodel = obj_nullmodel, 
                               known_loci = known_loci, variant_type = variant_type, 
                               geno_missing_imputation = geno_missing_imputation, 
                               QC_label = QC_label, Annotation_dir = Annotation_dir, 
                               Use_annotation_weights = Use_annotation_weights, 
                               Annotation_name_catalog = Annotation_name_catalog, 
                               Annotation_name = Annotation_name)
        ncRNA_sig_cond <- rbind(ncRNA_sig_cond, res_cond)
        seqClose(genofile)
      }
    }
    write.csv(ncRNA_sig_cond, file = paste0(ncRNA_output_path, 
                                            "ncRNA_sig_cond.csv"))
    noncoding_sig_cond <- rbind(noncoding_sig_cond, ncRNA_sig_cond)
    write.csv(noncoding_sig_cond, file = paste0(output_path, 
                                                "noncoding_sig_cond.csv"))
  }
  if (manhattan_plot) {
    results_STAAR <- results_UTR_genome[, c(1, 2, dim(results_UTR_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "UTR"
    results_STAAR <- results_upstream_genome[, c(1, 2, dim(results_upstream_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "upstream"
    results_STAAR <- results_downstream_genome[, c(1, 2, 
                                                   dim(results_downstream_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "downstream"
    results_STAAR <- results_promoter_CAGE_genome[, c(1, 
                                                      2, dim(results_promoter_CAGE_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_CAGE"
    results_STAAR <- results_promoter_DHS_genome[, c(1, 2, 
                                                     dim(results_promoter_DHS_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_DHS"
    results_STAAR <- results_enhancer_CAGE_genome[, c(1, 
                                                      2, dim(results_enhancer_CAGE_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_CAGE"
    results_STAAR <- results_enhancer_DHS_genome[, c(1, 2, 
                                                     dim(results_enhancer_DHS_genome)[2])]
    results_m <- c()
    for (i in 1:dim(results_STAAR)[2]) {
      results_m <- cbind(results_m, unlist(results_STAAR[, 
                                                         i]))
    }
    colnames(results_m) <- colnames(results_STAAR)
    results_m <- data.frame(results_m, stringsAsFactors = FALSE)
    results_m[, 2] <- as.numeric(results_m[, 2])
    results_m[, 3] <- as.numeric(results_m[, 3])
    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan, 
                                             results_m, by = c(chromosome_name = "Chr", hgnc_symbol = "Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_DHS"
    noncoding_minp <- min(genes_info_manhattan[, (dim(genes_info_manhattan)[2] - 
                                                    6):dim(genes_info_manhattan)[2]])
    if(is.na(min_y)){
      min_y <- ceiling(-log10(noncoding_minp)) + 1
    }else{
      min_y <- min_y
    }
    pch <- c(0, 1, 2, 3, 4, 5, 6)
    if (!gene_name){
    figure1 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$UTR, sig.level = alpha, 
                              pch = pch[1], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure2 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$upstream, sig.level = alpha, 
                              pch = pch[2], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure3 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$downstream, sig.level = alpha, 
                              pch = pch[3], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure4 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$promoter_CAGE, sig.level = alpha, 
                              pch = pch[4], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure5 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$promoter_DHS, sig.level = alpha, 
                              pch = pch[5], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure6 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$enhancer_CAGE, sig.level = alpha, 
                              pch = pch[6], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
    figure7 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$enhancer_DHS, sig.level = alpha, 
                              pch = pch[7], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                           "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                           "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
                                                                                                    
     print("Manhattan plot")
    png(paste0(output_path, "gene_centric_noncoding_manhattan.png"), 
        width = 9, height = 6, units = "in", res = 600)
    print(figure1)
    print(figure2, newpage = FALSE)
    print(figure3, newpage = FALSE)
    print(figure4, newpage = FALSE)
    print(figure5, newpage = FALSE)
    print(figure6, newpage = FALSE)
    print(figure7, newpage = FALSE)
    dev.off()}
    else {
      genes_info_manhattan$label1 <- ifelse(genes_info_manhattan$UTR < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure1 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$UTR, genes_info_manhattan$label1,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[1], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label2 <- ifelse(genes_info_manhattan$upstream < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure2 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$upstream, genes_info_manhattan$label2,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[2], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label3 <- ifelse(genes_info_manhattan$downstream < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure3 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$downstream, genes_info_manhattan$label3,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[3], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label4 <- ifelse(genes_info_manhattan$promoter_CAGE < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure4 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$promoter_CAGE, genes_info_manhattan$label4,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[4], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label5 <- ifelse(genes_info_manhattan$promoter_DHS < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure5 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$promoter_DHS, genes_info_manhattan$label5,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[5], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label6 <- ifelse(genes_info_manhattan$enhancer_CAGE < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure6 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$enhancer_CAGE, genes_info_manhattan$label6,sig.level = alpha, suggest.level = suggest.alpha, 
                                pch = pch[6], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      genes_info_manhattan$label7 <- ifelse(genes_info_manhattan$enhancer_DHS < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure7 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$enhancer_DHS, genes_info_manhattan$label7,sig.level = alpha, suggest.level = suggest.alpha,
                                pch = pch[7], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("UTR", 
                                                                                                                                                                      "upstream", "downstream", "promoter_CAGE", 
                                                                                                                                                                      "promoter_DHS", "enhancer_CAGE", "enhancer_DHS"))))
      
      print("Manhattan plot")
      png(paste0(output_path, "gene_centric_noncoding_manhattan.png"), 
          width = 9, height = 6, units = "in", res = 600)
      print(figure1)
      print(figure2, newpage = FALSE)
      print(figure3, newpage = FALSE)
      print(figure4, newpage = FALSE)
      print(figure5, newpage = FALSE)
      print(figure6, newpage = FALSE)
      print(figure7, newpage = FALSE)
      dev.off()}
    }
  if (QQ_plot) {
    print("Q-Q plot")
    cex_point <- 1
    png(paste0(output_path, "gene_centric_noncoding_qqplot.png"), 
        width = 9, height = 9, units = "in", res = 600)
    observed <- sort(genes_info_manhattan$UTR[genes_info_manhattan$UTR < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    plot(lexp, lobs, pch = 0, col='green',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$upstream[genes_info_manhattan$upstream < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 1, col='blue',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$downstream[genes_info_manhattan$downstream < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 2, col='brown',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$promoter_CAGE[genes_info_manhattan$promoter_CAGE < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 3, col='orange',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$promoter_DHS[genes_info_manhattan$promoter_DHS < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 4, col='black',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$enhancer_CAGE[genes_info_manhattan$enhancer_CAGE < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 5, col='slateblue',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$enhancer_DHS[genes_info_manhattan$enhancer_DHS < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    par(new = T)
    plot(lexp, lobs, pch = 6, col='cyan',cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    legend("topleft", legend = c("UTR", "upstream", "downstream", 
                                 "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", 
                                 "enhancer_DHS"), ncol = 1, bty = "o", box.lwd = 1, 
           pch = 0:6,col = c('green','blue','brown','orange','black','slateblue','cyan'), cex = 1, text.font = 2)
    dev.off()
  }
}
