rm(list = ls())
gc()

## real qpcr data
library(data.table)
library(dplyr)
library(ggbreak)
library(ggpubr)
library(ggprism)

## check data form of qpcr, should be a data.frame with needed columns
check_qpcr_data <- function(x) {
  # Check if the data is a data.frame
  if (!is.data.frame(x)) {
    stop("Error: The input is not a data.frame.")
  }
  
  # Check if the required columns are present
  required_columns <- c("Well", "Target", "Sample", "Cq")
  missing_columns <- setdiff(required_columns, colnames(x))
  if (length(missing_columns) > 0) {
    stop(paste("Error: The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }
  
  # Check if the Cq column is numeric
  if (!is.numeric(x$Cq)) {
    stop("Error: The 'Cq' column is not numeric.")
  }
}

tidy_qpcr <- function(x, Housekeeping_gene, Control, Exp){
  
  ## check data form of qpcr, should be a data.frame with needed columns
  check_qpcr_data(x)
  
  # Separate housekeeping gene (ACTIN) and other genes
  housekeeping <- x %>%
    filter(Target == Housekeeping_gene) %>%
    summarise(housekeeping_Ct = mean(Cq, na.rm = TRUE), .by = Sample)  
  # Calculate mean Cq for ACTIN
  
  ## △Ct =Gene△Ct- Housekeeping gene Ct）
  qpcr_with_deltaCt <- x %>%
    filter(Target != c(Housekeeping_gene)) %>%
    left_join(housekeeping, by = "Sample") %>%  # Join by Sample
    mutate(deltaCt = Cq - housekeeping_Ct)  # Calculate ΔCt
  
  ## 2^-△△Ct = Exp△Ct - Control△Ct(reference Ct)
  Con <- qpcr_with_deltaCt %>%
    filter(Sample == c(Control)) 
  Con <- Con %>% summarise(ref_deltaCt = mean(deltaCt, na.rm = TRUE), .by = Target) %>%
    right_join(Con, by = "Target") %>%
    mutate(deltadelta_Ct = deltaCt - ref_deltaCt,
           Fold_change = 2^-deltadelta_Ct)
  
  Exp <- qpcr_with_deltaCt %>%
    filter(Sample == c(Exp)) 
  
  
  if(all(Con$Target == Exp$Target)){
    Exp$deltadelta_Ct <- c(Exp$deltaCt - Con$ref_deltaCt)
    Exp$Fold_change <- 2^-(Exp$deltadelta_Ct)
  } else {
    stop("Target in Control is not equal to Target in Experiment!")
  }
  
  # Find overlapping columns
  overlap_cols <- intersect(names(Con), names(Exp))
  
  # Select only the overlapping columns and bind them
  combined_df <- rbind(
    Con[, overlap_cols, drop = FALSE], ## When you set drop = FALSE, the dimensions are preserved, and the result remains a dataframe or matrix.But not a vector
    Exp[, overlap_cols, drop = FALSE]
  )
  
  ## barplot with error bar
  combined_df$Fold_change <- round(combined_df$Fold_change, 4)
  
  return(combined_df)
}

plot_qpcr_data <- function(combined_df, facet = TRUE, color = "jco",
                           target_name = NULL, use_break = FALSE,
                           lowerbreak = NULL, upperbreak = NULL) {
  # Check if target_name is provided when facet = FALSE
  if (!facet && is.null(target_name)) {
    stop("Error: 'target_name' must be provided when facet = FALSE.")
  }
  
  # Filter the data if facet = FALSE and target_name is provided
  if (!facet) {
    combined_df <- combined_df %>%
      filter(Target == target_name)
  }
  
  # Create the basic plot
  p <- ggbarplot(combined_df,
                 x = 'Sample', 
                 y = 'Fold_change',
                 fill = "Sample",
                 palette = color,
                 add = c("mean_sd","jitter"),
                 add.params = list(color = "black", shape = 2),
                 ylab = 'Relative mRNA expression',
                 legend = 'none',
                 ggtheme = theme_prism()) + 
    stat_compare_means(aes(label = ..p.signif..),  # Add p-values
                       comparisons = list(unique(combined_df$Sample)), 
                       method = 't.test',
                       bracket.size = 0.5)
 # ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001
  
  # Apply y-axis break if use_break is TRUE and both lowerbreak and upperbreak are provided
  if (use_break) {
    # Ensure that both lowerbreak and upperbreak are provided
    if (is.null(lowerbreak) || is.null(upperbreak)) {
      stop("Error: 'lowerbreak' and 'upperbreak' must be provided when 'use_break' is TRUE.")
    }
    
    p <- p + scale_y_break(breaks = c(lowerbreak, upperbreak))
  }
  
  # Conditionally add facetting
  if (facet) {
    p <- p + facet_wrap(~ Target)
  } else {
    p <- p + ggtitle(target_name)
  }
  
  return(p)
}


## mimic qpcR data
qpcr <- data.frame(Well = c(paste0("A", seq(1:12)),
                            paste0("B", seq(1:12)),
                            paste0("C", seq(1:12))),
                   Target = c(rep("ACTIN",12),
                              rep("BCL2",12),
                              rep("TFRC",12)),
                   Sample = rep(c("TF-1", "TF-1_IFNa", "HEL", "HEL_IFNa"),times = 3, each =3),
                   Cq = c(16.09, 15.93, 15.91, 14.56, 14.67, 14.82, 15.18, 15.36, 15.28, 15.40,15.30, 15.38,
                          26.19,26.65,28.06,28.33,26.04,26.91,27.90,28.01,25.78,26.63,28.35,28.28,
                          19.06,19.05,19.15,17.17,17.15,17.31,19.67,20.35,19.15,17.31,19.76,20.26)
)
