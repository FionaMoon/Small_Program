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

sc_Alluvial_pl <- function(df , x_key, order_x, stratum_group, percent = NULL){

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
