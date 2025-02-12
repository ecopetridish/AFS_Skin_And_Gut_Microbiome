
##################################################################

# Functions needed for the pilot study data analysis #

##################################################################



# (1) From A. d, 01: Remove samples with low number of reads

plot_seq_depth<-function(df, plot_title){
  # number of reads per sample (sequencing depth)
  ggplot(df, aes(x=Sample_ID, y=sum)) +
    geom_bar(stat = "identity") +
    theme_classic()+
    labs(y=paste0("Number of reads (",plot_title,")")) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(labels = comma)
}

plot_reads_freq<-function(df, plot_title){
  # frequency of read counts
  ggplot(df, aes(x = sum)) + 
    geom_histogram(color = "black", fill = "firebrick3", binwidth = 150) +
    ggtitle(paste0("Distribution of sample sequencing depth (",plot_title,")")) + 
    xlab("Read counts") +
    ylab("# of Samples")
}




# (2) From B. ### 01. Create abundance tables 

make_abundance_tables<-function(level, Level) {
  
  # aggregate data frames dynamically using the specified 'level'
  TotalLevelCounts<-aggregate(TotalCount ~ ., data = asv_taxa[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_g<-aggregate(TotalCount ~ ., data = asv_taxa_g[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_s<-aggregate(TotalCount ~ ., data = asv_taxa_s[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_ec<-aggregate(TotalCount ~ ., data = asv_taxa_ec[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_gloves<-aggregate(TotalCount ~ ., data = asv_taxa_gloves[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_hands<-aggregate(TotalCount ~ ., data = asv_taxa_hands[, c(level, "TotalCount")], FUN = sum)
  TotalLevelCounts_others<-aggregate(TotalCount ~ ., data = asv_taxa_others[, c(level, "TotalCount")], FUN = sum)
  
  
  # rename columns for consistency
  colnames(TotalLevelCounts)[1]<-"Level"
  colnames(TotalLevelCounts_g)[1]<-"Level"
  colnames(TotalLevelCounts_s)[1]<-"Level"
  colnames(TotalLevelCounts_ec)[1]<-"Level"
  colnames(TotalLevelCounts_gloves)[1]<-"Level"
  colnames(TotalLevelCounts_hands)[1]<-"Level"
  colnames(TotalLevelCounts_others)[1]<-"Level"
  
  
  # calculate percentages
  TotalLevelCounts$Abundance <- round((TotalLevelCounts$TotalCount * 100) / sum(TotalLevelCounts$TotalCount), digits = 3)
  TotalLevelCounts <- dplyr::arrange(TotalLevelCounts, desc(TotalCount))
  
  TotalLevelCounts_g$Abundance <- round((TotalLevelCounts_g$TotalCount * 100) / sum(TotalLevelCounts_g$TotalCount), digits = 3)
  TotalLevelCounts_g <- dplyr::arrange(TotalLevelCounts_g, desc(TotalCount))
  
  TotalLevelCounts_s$Abundance <- round((TotalLevelCounts_s$TotalCount * 100) / sum(TotalLevelCounts_s$TotalCount), digits = 3)
  TotalLevelCounts_s <- dplyr::arrange(TotalLevelCounts_s, desc(TotalCount))
  
  TotalLevelCounts_ec$Abundance <- round((TotalLevelCounts_ec$TotalCount * 100) / sum(TotalLevelCounts_ec$TotalCount), digits = 3)
  TotalLevelCounts_ec <- dplyr::arrange(TotalLevelCounts_ec, desc(TotalCount))
  
  TotalLevelCounts_gloves$Abundance <- round((TotalLevelCounts_gloves$TotalCount * 100) / sum(TotalLevelCounts_gloves$TotalCount), digits = 3)
  TotalLevelCounts_gloves <- dplyr::arrange(TotalLevelCounts_gloves, desc(TotalCount))
  
  TotalLevelCounts_hands$Abundance <- round((TotalLevelCounts_hands$TotalCount * 100) / sum(TotalLevelCounts_hands$TotalCount), digits = 3)
  TotalLevelCounts_hands <- dplyr::arrange(TotalLevelCounts_hands, desc(TotalCount))
  
  TotalLevelCounts_others$Abundance <- round((TotalLevelCounts_others$TotalCount * 100) / sum(TotalLevelCounts_others$TotalCount), digits = 3)
  TotalLevelCounts_others <- dplyr::arrange(TotalLevelCounts_others, desc(TotalCount))
  
  # Count the number of ASVs for each level
  count_asvs <- function(data, level_col) {
    df <- data.frame(Level = character(), noASVs = integer(), stringsAsFactors = FALSE)
    unique_levels <- unique(data[[level_col]])
    
    for (i in unique_levels) {
      #len <- length(which(data[[level_col]] == i))
      len <- length(unique(data[data[[level_col]] == i, "ASV"]))
      df2 <- data.frame(Level = i, noASVs = len, stringsAsFactors = FALSE)
      df <- rbind(df, df2)
    }
    df$noASVs <- as.numeric(df$noASVs)
    return(df)
  }
  
  df_all <- count_asvs(asv_taxa, level)
  df_g <- count_asvs(asv_taxa_g, level)
  df_s <- count_asvs(asv_taxa_s, level)
  df_ec <- count_asvs(asv_taxa_ec, level)
  df_gloves <- count_asvs(asv_taxa_gloves, level)
  df_hands <- count_asvs(asv_taxa_hands, level)
  df_others <- count_asvs(asv_taxa_others, level)
  
  
  # Add this information to the abundance table
  TotalLevelCounts <- dplyr::left_join(TotalLevelCounts, df_all, by = "Level")
  TotalLevelCounts_g <- dplyr::left_join(TotalLevelCounts_g, df_g, by = "Level")
  TotalLevelCounts_s <- dplyr::left_join(TotalLevelCounts_s, df_s, by = "Level")
  TotalLevelCounts_ec <- dplyr::left_join(TotalLevelCounts_ec, df_ec, by = "Level")
  TotalLevelCounts_gloves <- dplyr::left_join(TotalLevelCounts_gloves, df_gloves, by = "Level")
  TotalLevelCounts_hands <- dplyr::left_join(TotalLevelCounts_hands, df_hands , by = "Level")
  TotalLevelCounts_others <- dplyr::left_join(TotalLevelCounts_others, df_others, by = "Level")
  
  
  # Rename columns after joining
  colnames(TotalLevelCounts) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_g) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_s) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_ec) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_gloves) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_hands) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  colnames(TotalLevelCounts_others) <- c("Level", "Read Count", "Abundance (%)", "No of ASVs")
  
  
  # Assign dynamically named variables to the global environment
  colnames(TotalLevelCounts)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts"), TotalLevelCounts, envir = .GlobalEnv)
  colnames(TotalLevelCounts_g)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_g"), TotalLevelCounts_g, envir = .GlobalEnv)
  colnames(TotalLevelCounts_s)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_s"), TotalLevelCounts_s, envir = .GlobalEnv)
  colnames(TotalLevelCounts_ec)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_ec"), TotalLevelCounts_ec, envir = .GlobalEnv)
  colnames(TotalLevelCounts_gloves)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_gloves"), TotalLevelCounts_gloves, envir = .GlobalEnv)
  colnames(TotalLevelCounts_hands)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_hands"), TotalLevelCounts_hands, envir = .GlobalEnv)
  colnames(TotalLevelCounts_others)[1]<-paste0(Level)
  assign(paste0("Total", Level, "Counts_others"), TotalLevelCounts_others, envir = .GlobalEnv)
  
  # Save as CSV files
  write.csv(TotalLevelCounts, file.path(output_dir, paste0("Total", Level, "Counts.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_g, file.path(output_dir, paste0("Total", Level, "Counts_g.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_s, file.path(output_dir, paste0("Total", Level, "Counts_s.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_ec, file.path(output_dir, paste0("Total", Level, "Counts_ec.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_gloves, file.path(output_dir, paste0("Total", Level, "Counts_gloves.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_hands, file.path(output_dir, paste0("Total", Level, "Counts_hands.csv")), row.names = FALSE)
  write.csv(TotalLevelCounts_others, file.path(output_dir, paste0("Total", Level, "Counts_others.csv")), row.names = FALSE)
  
}


# (3) From F. ### 02. Make forest plots

# function to calculate coefficients

#note: confidence intervals are calculated with beta +/- qt(0.975, df)*std.error

calculate_coefficients <- function(model_out, df, n) {
  for (i in 1:length(model_out)) {
    df[i, "beta"] <- round(model_out[[i]]$coefficients[2, 1], digits = 3)
    df[i, "std_err"] <- round(model_out[[i]]$coefficients[2, 2], digits = 3)
    df[i, "confint.low"] <- round(df[i, "beta"] - qt(0.975, n - 3) * df[i, "std_err"], digits = 3)
    df[i, "confint.up"] <- round(df[i, "beta"] + qt(0.975, n - 3) * df[i, "std_err"], digits = 3)
  }
  return(df)
}


# function to create forest plots
plot_forest <- function(df, title_left, title_right,col1,col2) {
  df$color <- ifelse(df$beta < 0, title_left, title_right)
  
  ggplot(df, aes(x = beta, y = taxa, color = color)) +
    geom_point(size = 3.5, shape=15) +
    geom_errorbarh(aes(xmin = confint.low, xmax = confint.up), linewidth = 1.2, height = 0) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#636363") +
    labs(x = "Beta coefficient (95% CI)", y = "", fill = "") +
    scale_color_manual(values = c(col1, col2)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(size = 0.6, linetype = 1, colour = "black"),
          axis.ticks = element_line(size = 0.5, linetype = 1, colour = "black"),
          axis.text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(size = 21, colour = "black", margin = margin(t=15))) +
    ggtitle("")
    #annotate("text", x = -5, y = Inf, hjust = 0, label = paste0("abundant in ",title_left), vjust = 1, size = 4) +
    #annotate("text", x = 2, y = Inf, hjust = 1, label = paste0("abundant in ",title_right), vjust = 1, size = 4, family = "Arial") +
    
    #theme(plot.title = element_text(size = 10, face = "bold", hjust = 0, lineheight = 1.5, family = "Arial"))
}
