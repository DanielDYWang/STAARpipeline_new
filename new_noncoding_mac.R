new_noncoding_mac <- function (chr, gene_name, genofile, obj_nullmodel, rare_maf_cutoff = 0.01, 
          rv_num_cutoff = 2, QC_label = "annotation/filter", variant_type = c("SNV", 
                                                                              "Indel", "variant"), geno_missing_imputation = c("mean", 
                                                                                                                               "minor"), Annotation_dir = "annotation/info/FunctionalAnnotation", 
          Annotation_name_catalog, Use_annotation_weights = c(TRUE, 
                                                              FALSE), Annotation_name = NULL, silent = FALSE, MAC_super_rare=MAC_super_rare, cMAC_super_rare = cMAC_super_rare) 
{
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  phenotype.id <- as.character(obj_nullmodel$id_include)
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter)
  gc()
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir, 
                                                  Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                      "GENCODE.Category")]))
  is.in <- (GENCODE.Category == "downstream") & (SNVlist)
  variant.id.downstream <- variant.id[is.in]
  seqSetFilter(genofile, variant.id = variant.id.downstream, 
               sample.id = phenotype.id)
  rm(variant.id.downstream)
  gc()
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir, 
                                              Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                  "GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split, function(z) length(z))
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV, variant_gene_num)
  rm(GENCODE.Info)
  gc()
  rm(variant_gene_num)
  gc()
  Gene <- as.character(unlist(GENCODE.Info.split))
  rm(GENCODE.Info.split)
  gc()
  seqResetFilter(genofile)
  is.in <- which(Gene == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_downstream <- c()
  if (class(pvalues) == "list") {
    results_temp <- rep(NA, 5)
    results_temp[3] <- "downstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_downstream <- rbind(results_downstream, results_temp)
  }
  if (!is.null(results_downstream)) {
    colnames(results_downstream) <- colnames(results_downstream, 
                                             do.NULL = FALSE, prefix = "col")
    colnames(results_downstream)[1:5] <- c("Gene name", "Chr", 
                                           "Category", "#SNV","cMAC")
    colnames(results_downstream)[(dim(results_downstream)[2] - 
                                    1):dim(results_downstream)[2]] <- c("ACAT-O", "STAAR-O")
  }
  seqResetFilter(genofile)
  is.in <- (GENCODE.Category == "upstream") & (SNVlist)
  variant.id.upstream <- variant.id[is.in]
  seqSetFilter(genofile, variant.id = variant.id.upstream, 
               sample.id = phenotype.id)
  rm(variant.id.upstream)
  gc()
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir, 
                                              Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                  "GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split, function(z) length(z))
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV, variant_gene_num)
  rm(GENCODE.Info)
  gc()
  rm(variant_gene_num)
  gc()
  Gene <- as.character(unlist(GENCODE.Info.split))
  rm(GENCODE.Info.split)
  gc()
  seqResetFilter(genofile)
  is.in <- which(Gene == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_upstream <- c()
  if (class(pvalues) == "list") {
    results_temp <- rep(NA, 5)
    results_temp[3] <- "upstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_upstream <- rbind(results_upstream, results_temp)
  }
  if (!is.null(results_upstream)) {
    colnames(results_upstream) <- colnames(results_upstream, 
                                           do.NULL = FALSE, prefix = "col")
    colnames(results_upstream)[1:5] <- c("Gene name", "Chr", 
                                         "Category", "#SNV","cMAC")
    colnames(results_upstream)[(dim(results_upstream)[2] - 
                                  1):dim(results_upstream)[2]] <- c("ACAT-O", "STAAR-O")
  }
  seqResetFilter(genofile)
  is.in <- ((GENCODE.Category == "UTR3") | (GENCODE.Category == 
                                              "UTR5") | (GENCODE.Category == "UTR5;UTR3")) & (SNVlist)
  variant.id.UTR <- variant.id[is.in]
  rm(GENCODE.Category)
  gc()
  seqSetFilter(genofile, variant.id = variant.id.UTR, sample.id = phenotype.id)
  rm(variant.id.UTR)
  gc()
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir, 
                                              Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                  "GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
  rm(GENCODE.Info)
  gc()
  Gene <- as.character(sapply(GENCODE.Info.split, function(z) z[1]))
  rm(GENCODE.Info.split)
  gc()
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  seqResetFilter(genofile)
  is.in <- which(Gene == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_UTR <- c()
  if (class(pvalues) == "list") {
    results_temp <- rep(NA, 5)
    results_temp[3] <- "UTR"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_UTR <- rbind(results_UTR, results_temp)
  }
  if (!is.null(results_UTR)) {
    colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, 
                                      prefix = "col")
    colnames(results_UTR)[1:5] <- c("Gene name", "Chr", "Category", 
                                    "#SNV","cMAC")
    colnames(results_UTR)[(dim(results_UTR)[2] - 1):dim(results_UTR)[2]] <- c("ACAT-O", 
                                                                              "STAAR-O")
  }
  seqResetFilter(genofile)
  varid <- seqGetData(genofile, "variant.id")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              "CAGE")]))
  CAGEBvt <- CAGEAnno != ""
  CAGEidx <- which(CAGEBvt, useNames = TRUE)
  seqSetFilter(genofile, variant.id = varid[CAGEidx])
  seqSetFilter(genofile, promGobj, intersect = TRUE)
  CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir, 
                                              Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                  "GENCODE.Info")]))
  CAGEGene <- unlist(lapply(strsplit(CAGEpromgene, "\\(|\\,|;|-"), 
                            `[[`, 1))
  CAGEvchr <- as.numeric(seqGetData(genofile, "chromosome"))
  CAGEvpos <- as.numeric(seqGetData(genofile, "position"))
  CAGEvref <- as.character(seqGetData(genofile, "$ref"))
  CAGEvalt <- as.character(seqGetData(genofile, "$alt"))
  dfPromCAGEVarGene <- data.frame(CAGEvchr, CAGEvpos, CAGEvref, 
                                  CAGEvalt, CAGEGene)
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]
  dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist, ]
  dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
  dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
  dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)
  seqResetFilter(genofile)
  rm(dfPromCAGEVarGene)
  gc()
  is.in <- which(dfPromCAGEVarGene.SNV[, 5] == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_promoter_CAGE <- c()
  if (class(pvalues) == "list") {
    results_temp <- dfPromCAGEVarGene.SNV[1, 1:5]
    results_temp[3] <- "promoter_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_promoter_CAGE <- rbind(results_promoter_CAGE, 
                                   results_temp)
  }
  if (!is.null(results_promoter_CAGE)) {
    colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, 
                                                do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_CAGE)[1:5] <- c("Gene name", 
                                              "Chr", "Category", "#SNV","cMAC")
    colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2] - 
                                       1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O", 
                                                                              "STAAR-O")
  }
  seqResetFilter(genofile)
  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir, 
                                           Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                               "DHS")]))
  rOCRsBvt <- rOCRsAnno != ""
  rOCRsidx <- which(rOCRsBvt, useNames = TRUE)
  seqSetFilter(genofile, variant.id = varid[rOCRsidx])
  seqSetFilter(genofile, promGobj, intersect = TRUE)
  rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir, 
                                               Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                   "GENCODE.Info")]))
  rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene, "\\(|\\,|;|-"), 
                             `[[`, 1))
  rOCRsvchr <- as.numeric(seqGetData(genofile, "chromosome"))
  rOCRsvpos <- as.numeric(seqGetData(genofile, "position"))
  rOCRsvref <- as.character(seqGetData(genofile, "$ref"))
  rOCRsvalt <- as.character(seqGetData(genofile, "$alt"))
  dfPromrOCRsVarGene <- data.frame(rOCRsvchr, rOCRsvpos, rOCRsvref, 
                                   rOCRsvalt, rOCRsGene)
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]
  dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist, ]
  dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
  dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
  dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)
  seqResetFilter(genofile)
  rm(dfPromrOCRsVarGene)
  gc()
  is.in <- which(dfPromrOCRsVarGene.SNV[, 5] == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_promoter_DHS <- c()
  if (class(pvalues) == "list") {
    results_temp <- dfPromrOCRsVarGene.SNV[1, 1:5]
    results_temp[3] <- "promoter_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_promoter_DHS <- rbind(results_promoter_DHS, results_temp)
  }
  if (!is.null(results_promoter_DHS)) {
    colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, 
                                               do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_DHS)[1:5] <- c("Gene name", 
                                             "Chr", "Category", "#SNV","cMAC")
    colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2] - 
                                      1):dim(results_promoter_DHS)[2]] <- c("ACAT-O", "STAAR-O")
  }
  seqResetFilter(genofile)
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir, 
                                                Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                    "GeneHancer")]))
  genehancer <- genehancerAnno != ""
  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              "CAGE")]))
  CAGE <- CAGEAnno != ""
  CAGEGeneHancervt <- CAGEAnno != "" & genehancerAnno != ""
  CAGEGeneHanceridx <- which(CAGEGeneHancervt, useNames = TRUE)
  seqSetFilter(genofile, variant.id = varid[CAGEGeneHanceridx])
  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir, 
                                               Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                   "GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet, "="), 
                                `[[`, 4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene, ";"), 
                                 `[[`, 1))
  enhancervchr <- as.numeric(seqGetData(genofile, "chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile, "position"))
  enhancervref <- as.character(seqGetData(genofile, "$ref"))
  enhancervalt <- as.character(seqGetData(genofile, "$alt"))
  dfHancerVarGene <- data.frame(enhancervchr, enhancervpos, 
                                enhancervref, enhancervalt, enhancer2GENE)
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]
  dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist, ]
  dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
  dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
  dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)
  seqResetFilter(genofile)
  rm(dfHancerVarGene)
  gc()
  is.in <- which(dfHancerVarGene.SNV[, 5] == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_enhancer_CAGE <- c()
  if (class(pvalues) == "list") {
    results_temp <- dfHancerVarGene.SNV[1, 1:5]
    results_temp[3] <- "enhancer_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_enhancer_CAGE <- rbind(results_enhancer_CAGE, 
                                   results_temp)
  }
  if (!is.null(results_enhancer_CAGE)) {
    colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, 
                                                do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_CAGE)[1:5] <- c("Gene name", 
                                              "Chr", "Category", "#SNV","cMAC")
    colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2] - 
                                       1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O", 
                                                                              "STAAR-O")
  }
  seqResetFilter(genofile)
  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir, 
                                           Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                               "DHS")]))
  rOCRs <- rOCRsAnno != ""
  rOCRsGeneHancervt <- rOCRsAnno != "" & genehancerAnno != 
    ""
  rOCRsGeneHanceridx <- which(rOCRsGeneHancervt, useNames = TRUE)
  seqSetFilter(genofile, variant.id = varid[rOCRsGeneHanceridx])
  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir, 
                                               Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                   "GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet, "="), 
                                `[[`, 4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene, ";"), 
                                 `[[`, 1))
  enhancervchr <- as.numeric(seqGetData(genofile, "chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile, "position"))
  enhancervref <- as.character(seqGetData(genofile, "$ref"))
  enhancervalt <- as.character(seqGetData(genofile, "$alt"))
  dfHancerVarGene <- data.frame(enhancervchr, enhancervpos, 
                                enhancervref, enhancervalt, enhancer2GENE)
  rm(varid)
  gc()
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]
  dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist, ]
  dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
  dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
  dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)
  seqResetFilter(genofile)
  rm(dfHancerVarGene)
  gc()
  is.in <- which(dfHancerVarGene.SNV[, 5] == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  try(pvalues <- STAAR_new(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff, MAC_super_rare = MAC_super_rare, cMAC_super_rare = cMAC_super_rare), 
      silent = silent)
  results_enhancer_DHS <- c()
  if (class(pvalues) == "list") {
    results_temp <- dfHancerVarGene.SNV[1, 1:5]
    results_temp[3] <- "enhancer_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    results_temp[5] <- pvalues$cMAC
    results_temp <- c(results_temp, pvalues$results_STAAR_S_1_25, 
                      pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                      pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                      pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                      pvalues$results_STAAR_O)
    results_enhancer_DHS <- rbind(results_enhancer_DHS, results_temp)
  }
  if (!is.null(results_enhancer_DHS)) {
    colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, 
                                               do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_DHS)[1:5] <- c("Gene name", 
                                             "Chr", "Category", "#SNV","cMAC")
    colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2] - 
                                      1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O", "STAAR-O")
  }
  seqResetFilter(genofile)
  results_noncoding <- list(upstream = results_upstream, downstream = results_downstream, 
                            UTR = results_UTR, promoter_CAGE = results_promoter_CAGE, 
                            promoter_DHS = results_promoter_DHS, enhancer_CAGE = results_enhancer_CAGE, 
                            enhancer_DHS = results_enhancer_DHS)
  seqResetFilter(genofile)
  return(results_noncoding)
}
