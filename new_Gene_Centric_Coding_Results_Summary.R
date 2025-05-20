new_Gene_Centric_Coding_Results_Summary <- function (agds_dir, gene_centric_coding_jobs_num, input_path, 
                                                     output_path, gene_centric_results_name, obj_nullmodel, known_loci = NULL, 
                                                     method_cond = c("optimal", "naive"), QC_label = "annotation/filter", 
                                                     geno_missing_imputation = c("mean", "minor"), variant_type = c("SNV", 
                                                                                                                    "Indel", "variant"), Annotation_dir = "annotation/info/FunctionalAnnotation", 
                                                     Annotation_name_catalog, Use_annotation_weights = FALSE, 
                                                     Annotation_name = NULL, alpha = 2.5e-06, manhattan_plot = FALSE, 
                                                     QQ_plot = FALSE, cMAC = 10,gene_name=FALSE,min_y=NA) 
{
  method_cond <- match.arg(method_cond)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  results_coding_genome <- c()
  for (kk in 1:gene_centric_coding_jobs_num) {
    print(kk)
    results_coding <- get(load(paste0(input_path, gene_centric_results_name, 
                                      "_", kk, ".Rdata")))
    results_coding_genome <- c(results_coding_genome, results_coding)
  }
  results_plof_genome <- c()
  results_plof_ds_genome <- c()
  results_missense_genome <- c()
  results_disruptive_missense_genome <- c()
  results_synonymous_genome <- c()
  for (kk in 1:length(results_coding_genome)) {
    results <- results_coding_genome[[kk]]
    if (is.null(results) == FALSE) {
      if (results[3] == "plof" & results[5] >= cMAC) {
        results_plof_genome <- rbind(results_plof_genome, 
                                     results)
      }
      if (results[3] == "plof_ds" & results[5] >= cMAC) {
        results_plof_ds_genome <- rbind(results_plof_ds_genome, 
                                        results)
      }
      if (results[3] == "missense" & results[5] >= cMAC) {
        results_missense_genome <- rbind(results_missense_genome, 
                                         results)
      }
      if (results[3] == "disruptive_missense" & results[5] >= cMAC) {
        results_disruptive_missense_genome <- rbind(results_disruptive_missense_genome, 
                                                    results)
      }
      if (results[3] == "synonymous" & results[5] >= cMAC) {
        results_synonymous_genome <- rbind(results_synonymous_genome, 
                                           results)
      }
    }
    if (kk%%1000 == 0) {
      print(kk)
    }
  }
  save(results_plof_genome, file = paste0(output_path, "plof.Rdata"))
  save(results_plof_ds_genome, file = paste0(output_path, "plof_ds.Rdata"))
  save(results_missense_genome, file = paste0(output_path, 
                                              "missense.Rdata"))
  save(results_disruptive_missense_genome, file = paste0(output_path, 
                                                         "disruptive_missense.Rdata"))
  save(results_synonymous_genome, file = paste0(output_path, 
                                                "synonymous.Rdata"))
  plof_sig <- results_plof_genome[results_plof_genome[, "STAAR-O"] < 
                                    alpha, , drop = FALSE]
  write.csv(plof_sig, file = paste0(output_path, "plof_sig.csv"))
  missense_sig <- results_missense_genome[results_missense_genome[, 
                                                                  "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(missense_sig, file = paste0(output_path, "missense_sig.csv"))
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[, 
                                                                        "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(synonymous_sig, file = paste0(output_path, "synonymous_sig.csv"))
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[, 
                                                               "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(plof_ds_sig, file = paste0(output_path, "plof_ds_sig.csv"))
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[, 
                                                                                                   "STAAR-O"] < alpha, , drop = FALSE]
  write.csv(disruptive_missense_sig, file = paste0(output_path, 
                                                   "disruptive_missense_sig.csv"))
  library(readr)
  library(plyr)
  plof_sig <- read_csv(file = paste0(output_path, "plof_sig.csv"))
  missense_sig <- read_csv(file = paste0(output_path, "missense_sig.csv"))
  plof_ds_sig <- read_csv(file = paste0(output_path, "plof_ds_sig.csv"))
  disruptive_missense_sig <- read_csv(file = paste0(output_path, "disruptive_missense_sig.csv"))
  synonymous_sig <- read_csv(file = paste0(output_path, "synonymous_sig.csv"))
  
  coding_sig <- rbind.fill(missense_sig,plof_sig)
  coding_sig <- rbind.fill(coding_sig, synonymous_sig)
  coding_sig <- rbind.fill(coding_sig, plof_ds_sig)
  coding_sig <- rbind.fill(coding_sig, disruptive_missense_sig)
  
  write.csv(coding_sig, file = paste0(output_path, "coding_sig.csv"))
  
  if (length(known_loci) != 0) {
    plof_sig_cond <- c()
    if (length(plof_sig) != 0) {
      if ((class(plof_sig)[1] != "matrix") & (class(plof_sig)[1] != 
                                              "data.frame")) {
        plof_sig <- matrix(plof_sig, nrow = 1)
      }
      for (k in 1:dim(plof_sig)[1]) {
        chr <- as.numeric(plof_sig[k, 2])
        gene_name <- as.character(plof_sig[k, 1])
        category <- as.character(plof_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Coding_cond(chr = chr, 
                                             gene_name = gene_name, category = category, 
                                             genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                             known_loci = known_loci, variant_type = variant_type, 
                                             geno_missing_imputation = geno_missing_imputation, 
                                             QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                             Annotation_name_catalog = Annotation_name_catalog, 
                                             Use_annotation_weights = Use_annotation_weights, 
                                             Annotation_name = Annotation_name)
        plof_sig_cond <- rbind(plof_sig_cond, res_cond)
        seqClose(genofile)
      }
    }
    write.csv(plof_sig_cond, file = paste0(output_path, "plof_sig_cond.csv"))
    missense_sig_cond <- c()
    if (length(missense_sig) != 0) {
      if ((class(missense_sig)[1] != "matrix") & (class(missense_sig)[1] != 
                                                  "data.frame")) {
        missense_sig <- matrix(missense_sig, nrow = 1)
      }
      for (k in 1:dim(missense_sig)[1]) {
        chr <- as.numeric(missense_sig[k, 2])
        gene_name <- as.character(missense_sig[k, 1])
        category <- as.character(missense_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Coding_cond(chr = chr, 
                                             gene_name = gene_name, category = category, 
                                             genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                             known_loci = known_loci, variant_type = variant_type, 
                                             geno_missing_imputation = geno_missing_imputation, 
                                             QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                             Annotation_name_catalog = Annotation_name_catalog, 
                                             Use_annotation_weights = Use_annotation_weights, 
                                             Annotation_name = Annotation_name)
        missense_sig_cond <- rbind(missense_sig_cond, 
                                   res_cond)
        seqClose(genofile)
      }
    }
    write.csv(missense_sig_cond, file = paste0(output_path, 
                                               "missense_sig_cond.csv"))
    synonymous_sig_cond <- c()
    if (length(synonymous_sig) != 0) {
      if ((class(synonymous_sig)[1] != "matrix") & (class(synonymous_sig)[1] != 
                                                    "data.frame")) {
        synonymous_sig <- matrix(synonymous_sig, nrow = 1)
      }
      for (k in 1:dim(synonymous_sig)[1]) {
        chr <- as.numeric(synonymous_sig[k, 2])
        gene_name <- as.character(synonymous_sig[k, 1])
        category <- as.character(synonymous_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Coding_cond(chr = chr, 
                                             gene_name = gene_name, category = category, 
                                             genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                             known_loci = known_loci, variant_type = variant_type, 
                                             geno_missing_imputation = geno_missing_imputation, 
                                             QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                             Annotation_name_catalog = Annotation_name_catalog, 
                                             Use_annotation_weights = Use_annotation_weights, 
                                             Annotation_name = Annotation_name)
        synonymous_sig_cond <- rbind(synonymous_sig_cond, 
                                     res_cond)
        seqClose(genofile)
      }
    }
    write.csv(synonymous_sig_cond, file = paste0(output_path, 
                                                 "synonymous_sig_cond.csv"))
    plof_ds_sig_cond <- c()
    if (length(plof_ds_sig) != 0) {
      if ((class(plof_ds_sig)[1] != "matrix") & (class(plof_ds_sig)[1] != 
                                                 "data.frame")) {
        plof_ds_sig <- matrix(plof_ds_sig, nrow = 1)
      }
      for (k in 1:dim(plof_ds_sig)[1]) {
        chr <- as.numeric(plof_ds_sig[k, 2])
        gene_name <- as.character(plof_ds_sig[k, 1])
        category <- as.character(plof_ds_sig[k, 3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Coding_cond(chr = chr, 
                                             gene_name = gene_name, category = category, 
                                             genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                             known_loci = known_loci, variant_type = variant_type, 
                                             geno_missing_imputation = geno_missing_imputation, 
                                             QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                             Annotation_name_catalog = Annotation_name_catalog, 
                                             Use_annotation_weights = Use_annotation_weights, 
                                             Annotation_name = Annotation_name)
        plof_ds_sig_cond <- rbind(plof_ds_sig_cond, res_cond)
        seqClose(genofile)
      }
    }
    write.csv(plof_ds_sig_cond, file = paste0(output_path, 
                                              "plof_ds_sig_cond.csv"))
    disruptive_missense_sig_cond <- c()
    if (length(disruptive_missense_sig) != 0) {
      if ((class(disruptive_missense_sig)[1] != "matrix") & 
          (class(disruptive_missense_sig)[1] != "data.frame")) {
        disruptive_missense_sig <- matrix(disruptive_missense_sig, 
                                          nrow = 1)
      }
      for (k in 1:dim(disruptive_missense_sig)[1]) {
        chr <- as.numeric(disruptive_missense_sig[k, 
                                                  2])
        gene_name <- as.character(disruptive_missense_sig[k, 
                                                          1])
        category <- as.character(disruptive_missense_sig[k, 
                                                         3])
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        res_cond <- Gene_Centric_Coding_cond(chr = chr, 
                                             gene_name = gene_name, category = category, 
                                             genofile = genofile, obj_nullmodel = obj_nullmodel, 
                                             known_loci = known_loci, variant_type = variant_type, 
                                             geno_missing_imputation = geno_missing_imputation, 
                                             QC_label = QC_label, Annotation_dir = Annotation_dir, 
                                             Annotation_name_catalog = Annotation_name_catalog, 
                                             Use_annotation_weights = Use_annotation_weights, 
                                             Annotation_name = Annotation_name)
        disruptive_missense_sig_cond <- rbind(disruptive_missense_sig_cond, 
                                              res_cond)
        seqClose(genofile)
      }
    }
    write.csv(disruptive_missense_sig_cond, file = paste0(output_path, 
                                                          "disruptive_missense_sig_cond.csv"))
    coding_sig_cond <- rbind(plof_sig_cond, missense_sig_cond)
    coding_sig_cond <- rbind(coding_sig_cond, synonymous_sig_cond)
    coding_sig_cond <- rbind(coding_sig_cond, plof_ds_sig_cond)
    coding_sig_cond <- rbind(coding_sig_cond, disruptive_missense_sig_cond)
    write.csv(coding_sig_cond, file = paste0(output_path, 
                                             "coding_sig_cond.csv"))
  }
  if (manhattan_plot) {
    results_STAAR <- results_plof_genome[, c(1, 2, dim(results_plof_genome)[2])]
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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof"
    results_STAAR <- results_plof_ds_genome[, c(1, 2, dim(results_plof_ds_genome)[2])]
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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"
    results_STAAR <- results_missense_genome[, c(1, 2, dim(results_missense_genome)[2])]
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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"
    results_STAAR <- results_disruptive_missense_genome[, 
                                                        c(1, 2, dim(results_missense_genome)[2])]
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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"
    results_STAAR <- results_synonymous_genome[, c(1, 2, 
                                                   dim(results_synonymous_genome)[2])]
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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"
    coding_minp <- min(genes_info_manhattan[, (dim(genes_info_manhattan)[2] - 
                                                 4):dim(genes_info_manhattan)[2]])
    if(is.na(min_y)){
    min_y <- ceiling(-log10(coding_minp)) + 1
    }else{
      min_y <- min_y
    }
    pch <- c(0, 1, 2, 3, 4)
    if (!gene_name){
    figure1 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$plof, sig.level = alpha, 
                              pch = pch[1], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                           "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                           "Synonymous"))))
    figure2 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$plof_ds, sig.level = alpha, 
                              pch = pch[2], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                           "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                           "Synonymous"))))
    figure3 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$missense, sig.level = alpha, 
                              pch = pch[3], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                           "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                           "Synonymous"))))
    figure4 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$disruptive_missense, 
                              sig.level = alpha, pch = pch[4], col = c("lightseagreen","chocolate"), ylim = c(0, min_y), auto.key = T, 
                              key = list(space = "top", columns = 5, title = "Functional Category", 
                                         cex.title = 1, points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                "Synonymous"))))
    figure5 <- manhattan_plot(genes_info_manhattan[, 2], 
                              (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                4])/2, genes_info_manhattan$synonymous, sig.level = alpha, 
                              pch = pch[5], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                  min_y), auto.key = T, key = list(space = "top", 
                                                                                                                   columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                   points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                           "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                           "Synonymous"))))
                                                                                                                                  
    
    print("Manhattan plot")
    png(paste0(output_path, "gene_centric_coding_manhattan.png"), 
        width = 9, height = 6, units = "in", res = 600)
    print(figure1)
    print(figure2, newpage = FALSE)
    print(figure3, newpage = FALSE)
    print(figure4, newpage = FALSE)
    print(figure5, newpage = FALSE)
    dev.off()}
    else{
      genes_info_manhattan$label1 <- ifelse(genes_info_manhattan$plof < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure1 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$plof, genes_info_manhattan$label1,sig.level = alpha, 
                                pch = pch[1], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                                      "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                                      "Synonymous"))))
      genes_info_manhattan$label2 <- ifelse(genes_info_manhattan$plof_ds < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure2 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$plof_ds, genes_info_manhattan$label2,sig.level = alpha, 
                                pch = pch[2], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                                      "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                                      "Synonymous"))))
      genes_info_manhattan$label3 <- ifelse(genes_info_manhattan$missense < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure3 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$missense, genes_info_manhattan$label3,sig.level = alpha, 
                                pch = pch[3], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                                      "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                                      "Synonymous"))))
      genes_info_manhattan$label4 <- ifelse(genes_info_manhattan$disruptive_missense < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure4 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$disruptive_missense,genes_info_manhattan$label4, 
                                sig.level = alpha, pch = pch[4], col = c("lightseagreen","chocolate"), ylim = c(0, min_y), auto.key = T, 
                                key = list(space = "top", columns = 5, title = "Functional Category", 
                                           cex.title = 1, points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                  "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                  "Synonymous"))))
      genes_info_manhattan$label5 <- ifelse(genes_info_manhattan$synonymous < alpha, genes_info_manhattan$hgnc_symbol, '')
      figure5 <- manhattan_plot_new(genes_info_manhattan[, 2], 
                                (genes_info_manhattan[, 3] + genes_info_manhattan[, 
                                                                                  4])/2, genes_info_manhattan$synonymous, genes_info_manhattan$label5,sig.level = alpha, 
                                pch = pch[5], col = c("lightseagreen","chocolate"), ylim = c(0, 
                                                                                             min_y), auto.key = T, key = list(space = "top", 
                                                                                                                              columns = 5, title = "Functional Category", cex.title = 1, 
                                                                                                                              points = TRUE, pch = pch, text = list(c("pLoF", 
                                                                                                                                                                      "pLoF+D", "Missense", "Disruptive Missense", 
                                                                                                                                                                      "Synonymous"))))
      print("Manhattan plot")
      png(paste0(output_path, "gene_centric_coding_manhattan.png"), 
          width = 9, height = 6, units = "in", res = 600)
      print(figure1)
      print(figure2, newpage = FALSE)
      print(figure3, newpage = FALSE)
      print(figure4, newpage = FALSE)
      print(figure5, newpage = FALSE)
      dev.off()
    }
  }
  if (QQ_plot) {
    print("Q-Q plot")
    cex_point <- 1
    png(paste0(output_path, "gene_centric_coding_qqplot.png"), 
        width = 9, height = 9, units = "in", res = 600)
    observed <- sort(genes_info_manhattan$plof[genes_info_manhattan$plof < 1])
    lobs <- -(log10(observed))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected/(length(expected) + 1)))
    plot(lexp, lobs, pch = 0, col='green', cex = cex_point, xlim = c(0, 
                                                        5), ylim = c(0, min_y), xlab = expression(Expected ~ 
                                                                                                    ~-log[10](italic(p))), ylab = expression(Observed ~ 
                                                                                                                                               ~-log[10](italic(p))), font.lab = 2, cex.lab = 1, 
         cex.axis = 1, font.axis = 2)
    abline(0, 1, col = "red", lwd = 1)
    observed <- sort(genes_info_manhattan$plof_ds[genes_info_manhattan$plof_ds < 1])
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
    observed <- sort(genes_info_manhattan$missense[genes_info_manhattan$missense < 1])
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
    observed <- sort(genes_info_manhattan$disruptive_missense[genes_info_manhattan$disruptive_missense < 1])
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
    observed <- sort(genes_info_manhattan$synonymous[genes_info_manhattan$synonymous < 1])
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
    legend("topleft", legend = c("pLoF", "pLoF+D", "Missense", 
                                 "Disruptive Missense", "Synonymous"), ncol = 1, bty = "o", 
           box.lwd = 1, pch = 0:4, col=c('green','blue','brown','orange','black'),cex = 1, text.font = 2)
    dev.off()
  }
}