# Functions used in the analysis part of the project 

# accumulation_curve

# This functions calculates the accumulation curve of the eDNA dataset
# You input a dataframe (x) and need to have at least three columns: 
# a station column which you can choose (column_station: default is sample_name_all_pcr)
# a species_unit column which you can choose (species_unit: defaut is sequence)
# a count_reads columns (abundance of reads at the finest level, i.e. PCR level)
# you can also specify the method used for the accumulation curve, from the specaccum function in {vegan}
# 
# the output is a list with each element of the list an output of parts of the function 
# matrix is the species/station matrix
# accumulation is the output of specaccum
# accumulation_plot is a dataframe from accumulation for easy plotting 
# asymptote is the asymptote calculated from the multimodel inference 
# complete_samples: list of all the samples necessary for the accumulation curve. If some samples are discarded, use it to provide the entire list of samples to have a proper accumulation. 

accumulation_curve_all_samples <- function(x, column_station = "sample_name_all_pcr", complete_samples,
                               species_unit = "sequence", method_accumulation = "exact"){
  
  # lib 
  library(vegan)
  
  # Matrix species / station
  df_acc <- x %>%
    group_by(.dots = c(column_station, species_unit)) %>% 
    summarise(n_reads = sum(count_reads)) %>% 
    ungroup() %>% 
    tidyr::spread(species_unit, n_reads) %>% 
    as.data.frame() %>% 
    replace(is.na(.), 0)
  
  # Add the missing samples 
  add_samples <- complete_samples[complete_samples %ni% df_acc$sample_name_all_pcr]
  
  # Create the add matrix & convert it to numeric 
  add_matrix <- matrix(ncol= ncol(df_acc), nrow=length(add_samples))
  colnames(add_matrix) <- colnames(df_acc)
  add_matrix[,1] <- add_samples
  add_matrix[is.na(add_matrix)] <- 0
  add_matrix <- data.frame(add_matrix, stringsAsFactors = F)
  add_matrix <- add_matrix %>% mutate_at(vars(colnames(add_matrix)[-1]), ~as.numeric(as.character(.)))  
  
  # bind both 
  df_acc <- rbind(df_acc, add_matrix)
  
  # accumulation curve result
  acc_edna <- specaccum(df_acc[,2:ncol(df_acc)], 
                        method = method_accumulation, permutations = 100,
                        conditioned =TRUE, gamma = "jack1")
  
  # Accumulation for plot 
  acc_edna_df <- data.frame(richness = acc_edna$richness, sd = acc_edna$sd, sites = acc_edna$sites)
  
  # Asymptote 
  akaike.weights <- function(x)
  {
    x <- x[!is.na(x)]
    delta.aic <- x - min(x, na.rm = TRUE)
    rel.LL <- exp(-0.5 * delta.aic)
    sum.LL <- sum(rel.LL, na.rm = TRUE)
    weights.aic <- rel.LL/sum.LL
    return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
  }
  
  # models asymptotes
  lomo_edna <- fitspecaccum(acc_edna, "lomolino")
  aic_lomo_edna <-AIC(lomo_edna)
  mm_edna <- fitspecaccum(acc_edna, "michaelis-menten")
  aic_mm_edna <- AIC(mm_edna)
  gom_edna <- fitspecaccum(acc_edna, "gompertz")
  aic_gom_edna <- AIC(gom_edna)
  asy_edna <- fitspecaccum(acc_edna, "asymp")
  aic_asy_edna <- AIC(asy_edna)
  gis_edna <- fitspecaccum(acc_edna, "logis")
  aic_gis_edna <- AIC(gis_edna)
  
  # compute results
  res_edna <- matrix(NA,nrow = 5, ncol = 3)
  rownames(res_edna) <- c("lomolino", "michaelis-menten", "gompertz", "asymp", "logis")
  colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
  res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna)
  res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna))$weights
  res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(gom_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])
  
  # Calculation asymptote value
  asymp_edna <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])
  
  # Return
  return(list(matrix = df_acc, 
              accumulation = acc_edna,
              accumulation_plot = acc_edna_df, 
              asymptote = asymp_edna))
  
}

# Only the accumulation for plot 
accumulation_curve_df <- function(x, column_station = "sample_name_all_pcr", 
                                   species_unit = "sequence", method_accumulation = "exact"){
  
  # lib 
  library(vegan)
  
  # Matrix species / station
  df_acc <- x %>%
    group_by(.dots = c(column_station, species_unit)) %>% 
    summarise(n_reads = sum(count_reads)) %>% 
    ungroup() %>% 
    tidyr::spread(species_unit, n_reads) %>% 
    as.data.frame() %>% 
    replace(is.na(.), 0)
  
  # accumulation curve result
  acc_edna <- specaccum(df_acc[,2:ncol(df_acc)], 
                        method = method_accumulation, permutations = 100,
                        conditioned =TRUE, gamma = "jack1")
  
  # Accumulation for plot 
  acc_edna_df <- data.frame(richness = acc_edna$richness, sd = acc_edna$sd, sites = acc_edna$sites)
  
  # return 
  return(acc_edna_df)
  
}


# asymptote multi model
asymptote_mm <- function(x, column_station = "sample_name_all_pcr", 
                                  species_unit = "sequence", method_accumulation = "exact"){
  
  # lib
  require(vegan)
  
  # Matrix species / station
  df_acc <- x %>%
    group_by(.dots = c(column_station, species_unit)) %>% 
    summarise(n_reads = sum(count_reads)) %>% 
    ungroup() %>% 
    tidyr::spread(species_unit, n_reads) %>% 
    as.data.frame() %>% 
    replace(is.na(.), 0)
  
  # accumulation curve result
  acc_edna <- specaccum(df_acc[,2:ncol(df_acc)], 
                        method = method_accumulation, permutations = 100,
                        conditioned =TRUE, gamma = "jack1")
  
  # Accumulation for plot 
  acc_edna_df <- data.frame(richness = acc_edna$richness, sd = acc_edna$sd, sites = acc_edna$sites)
  
  # Asymptote 
  akaike.weights <- function(x)
  {
    x <- x[!is.na(x)]
    delta.aic <- x - min(x, na.rm = TRUE)
    rel.LL <- exp(-0.5 * delta.aic)
    sum.LL <- sum(rel.LL, na.rm = TRUE)
    weights.aic <- rel.LL/sum.LL
    return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
  }
  
  # models asymptotes
  lomo_edna <- fitspecaccum(acc_edna, "lomolino")
  aic_lomo_edna <-AIC(lomo_edna)
  mm_edna <- fitspecaccum(acc_edna, "michaelis-menten")
  aic_mm_edna <- AIC(mm_edna)
  gom_edna <- fitspecaccum(acc_edna, "gompertz")
  aic_gom_edna <- AIC(gom_edna)
  asy_edna <- fitspecaccum(acc_edna, "asymp")
  aic_asy_edna <- AIC(asy_edna)
  gis_edna <- fitspecaccum(acc_edna, "logis")
  aic_gis_edna <- AIC(gis_edna)
  
  # compute results
  res_edna <- matrix(NA,nrow = 5, ncol = 3)
  rownames(res_edna) <- c("lomolino", "michaelis-menten", "gompertz", "asymp", "logis")
  colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
  res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna)
  res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna))$weights
  res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(gom_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])
  
  # Calculation asymptote value
  asymp_edna <- data.frame(asymptote = weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"]), slope= coef(lomo_edna)[[3]])
  
  # Return 
  return(asymp_edna)
  
}


