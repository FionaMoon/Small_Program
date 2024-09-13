## load packages
library(ggalluvial)
library(ggnewscale)
library(plyr);library(dplyr)

## select color for final group
c4 <- c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d')
c_light12 <- c("#FFD966","#4666FF","#b33939","#CCE593","#FF7F00","#E06666", "#A8CAEA", "#B29F66","#FCE5CD", "#F4CCCC", "#AC93E5", "#6AA84F")
c15 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "skyblue2","gold1", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
c35 <- c(
  '#FFC312','#12CBC4','#ED4C67','#A3CB38',
  '#F79F1F','#1289A7','#D980FA','#B53471','#cd6133',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471', 
  '#EA2027','#006266','#5758BB','#d1ccc0','#C4E538',
  '#40407a','#706fd3','#34ace0','#FDA7DF',
  '#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#ffda79','#ccae62',
  '#b33939','#84817a','#cc8e35','#33d9b2')

## pt_Alluvial
theme_no_background <- function(base.theme = theme_bw(), pos = "right", legend_size = 12, axis_size = 10, title_size = 12) {
  base.theme %+replace%
    theme(legend.position = pos,
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = (axis_size + 1), face = "bold", angle = 90, hjust = 0.5),
            axis.text.x = element_text(size = (axis_size + 1), face = "bold", ),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_size),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = axis_size),
            axis.title = element_text(size = title_size, face = "bold"))
}

## x_key, the x axis of your sankey plot
## stratum_group, the stack barplot group of your each barplot
## df should be a data.frame which contains the colunm of `x_key` and `stratum_group`

sc_Alluvial_pl <- function(df , x_key, order_x, stratum_group, percent = NULL, color_use = NULL){

    if(is.null(percent) == T){
        x <- c(x_key, stratum_group)
        df_plot <- as.data.frame(table(df[,x])) ## calculate Freq
        colnames(df_plot) <- c("grp1", "grp2", "Freq")
        df_plot <- ddply(df_plot,.(grp1),transform,percent=Freq/sum(Freq)*100) 
    } else {
        x <- c(x_key, stratum_group, percent)
        df_plot <- as.data.frame(df[,x]) 
        colnames(df_plot) <- c("grp1", "grp2", "percent")
    }

    df_plot$label = paste0(sprintf("%.1f", df_plot$percent), "%")

    if(all(is.na(color_use)) == T){
            n_types <- length(unique(df_plot$grp2))
        if(n_types < 5){
                color_list = c4[1:n_types]
        } else if(n_types < 13){
            color_list = c_light12[1:n_types]
        } else if(n_types < 16){
            color_list = c15[1:n_types]
        } else if(n_types < 37){
            color_list = c36[1:n_types]    
        } else {
            stop("Out of color boundery")
        }
    col_define = color_list
    names(col_define) = unique(df_plot$grp2)
    } else {
        if(all(is.na(names(color_use))) == T){
            stop("color_use should be given as a named vector! \n named by stratum_group")
        } else {
            col_define = color_use
        }
    }


    x <- is_alluvia_form(as.data.frame(df_plot), axes = 1:(ncol(df_plot)-1), silent = TRUE)


    df_plot$grp1 <-  factor(df_plot$grp1, levels = order_x)


    if(x == T){
        # Alluvial plot without legend
    p <- ggplot(df_plot,
        aes(x = grp1, stratum = grp2, alluvium = grp2,
           y = percent, fill = grp2)) +
    geom_flow(aes(fill = grp2), alpha = .3) +
    geom_stratum(aes(color = grp2), alpha = .9) + 
    # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(values = col_define, breaks = names(col_define), aesthetics = c("color", "fill")) +
    theme_no_background()
    }

    return(p)
}



# # # example for you to learn
# data(majors)
# majors$curriculum <- as.factor(majors$curriculum)
# p <- sc_Alluvial_pl(majors, x_key = "semester", 
#     order_x = c("CURR1", "CURR3", "CURR5", "CURR7", "CURR9", "CURR11", "CURR13", "CURR15"),
# stratum_group = "curriculum")

# # if we already calculate percentage
# df <- read.table(text = "6  0.5 4   0.95    0.001   0.001   0.0012  0.01    0.0013  0.0025  91  0.75    0.95
# 8.1 2.1 2.9 0.7 5   2   4.5 1   0.5 0.7 90  0.9 0.2
# 8.1 2.1 2.9 0.7 5   2   4.5 1   0.5 0.7 90  0.9 0.2
# 5.3 2.9 2.2 0.3 3.2 3.8 2.1 25  16  7   56  0.03    0.001
# ")
# rownames(df) <- c("7w", "8w", "11w", "17w")
# colnames(df) <- c("CD33", "CD14", "Granulocytes", "Neutrophils", "DCs", "M1", "M2", "CD3", "CD4", "CD8", "B cell", "NK cell", "Plasma cell")
# library(tidyr)
# library(dplyr)
# df <- df %>% mutate(time = rownames(df))
# huNSGS <- df %>%
#   pivot_longer(!time, names_to = "Celltypes", values_to = "percentage")

# p <- sc_Alluvial_pl(huNSGS,x_key = "time", order_x = c("7w", "8w", "11w", "17w"),
# stratum_group = "Celltypes", percent = "percentage")

## if we define color to stratum
# data(majors)
# majors$curriculum <- factor(majors$curriculum)
# colors <- c15[1:length(unique(majors$curriculum))]
# names(colors) <- levels(majors$curriculum)
# color
#
# p <- sc_Alluvial_pl(majors, x_key = "semester", 
#     order_x = c("CURR1", "CURR3", "CURR5", "CURR7", "CURR9", "CURR11", "CURR13", "CURR15"),
# stratum_group = "curriculum",
# color_use = colors)

## This function can be used with ggplot2 to enhance the appearance of plots.
