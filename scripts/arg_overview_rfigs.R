###############################################
###############################################
### CREATING PANELS FOR ARG OVERVIEW FIGURE ###
###############################################
###############################################

####################
### SCRIPT SETUP ###
####################
library(here)
library(dplyr)
library(ggplot2)



############################################
### PANEL C:  LINEAGE COUNT THROUGH TIME ###
############################################
arg_attribute_df <- data.frame(time = 1:8,
                               lineage_count = c(4, 5, 4, 3, 2, 3, 2, 1),
                               point_color = NA,
                               local_tree_node_count = c(12, 12, 11, 9, 6, 6, 5, 3),
                               shared_node_count = as.integer(c(4, 5, 4, 3, 2, 3, 2, 1)))

for (i in seq_len(nrow(arg_attribute_df))) {
  if (i == 1) {
    arg_attribute_df[i, 'point_color'] <- '#e3e3e3'
  } else if ((arg_attribute_df[i, 'lineage_count'] < arg_attribute_df[i - 1, 'lineage_count'])) {
    arg_attribute_df[i, 'point_color'] <- '#5c5c5c'
  } else {
    arg_attribute_df[i, 'point_color'] <- '#fe0000'
  }
}

lineage_count_plot <- arg_attribute_df %>% 
  ggplot() +
  geom_line(aes(y = lineage_count,
                x = time),
            linewidth = 2) +
  geom_point(aes(y = lineage_count,
                 x = time),
             color = arg_attribute_df$point_color,
             size = 10) +
  ylab('Number of lineages') +
  xlab('Time') +
  theme_classic() +
  theme(panel.background = element_rect(fill ='transparent'),
        plot.background = element_rect(fill ='transparent', color = NA),
        axis.title.x = element_text(size = 29),
        axis.text.x = element_text(size = 22),
        #panel.grid.major.x = element_line(linetype = 'solid', size = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 1.5),
        axis.ticks.x = element_line(linewidth = 1.5)
        #axis.line.y = element_blank(),
        #axis.ticks.y = element_blank()
  ) +
  coord_flip(clip = "off")

ggsave(filename = here('figures', 'pdf', 'lineage_count_fig1.png'),
       plot = lineage_count_plot,
       bg = 'transparent',
       width = 3,
       height = 10,
       units = 'in',
       device = "png")



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

### PANEL C: VISUALIZATION OF EDGES  ###

# edge_df <- data.frame(edge = as.character(1:11),
#                       genome_start = c(1, 1, 1, 3, 1, 1, 1, 1, 1, 2, 1),
#                       genome_end =   c(4, 3, 4, 4, 4, 4, 4, 2, 4, 4, 4),
#                       time_start =   c(0, 4, 6, 4, 4, 2, 0, 0, 2, 0, 0),
#                       time_end =     c(4, 7, 7, 6, 6, 4, 2, 2, 4, 3, 3)) %>% 
#   arrange(time_start, time_end, genome_start)
# 
# coord_offset <- 0
# edge_vis_list <- list()
# for (i in 1:nrow(edge_df)) {
#   edge_vis_list[[i]] <- data.frame(xstart = seq(edge_df$genome_start[i], 
#                                                 edge_df$genome_end[i] - 0.01, 
#                                                 by = 0.01),
#                                    xstartcoord = seq(edge_df$genome_start[i] - edge_df$genome_start[i], 
#                                                      edge_df$genome_end[i] - edge_df$genome_start[i] - 0.01, 
#                                                      by = 0.01) + coord_offset,
#                                    xend = seq(edge_df$genome_start[i] - edge_df$genome_start[i], 
#                                               edge_df$genome_end[i] - edge_df$genome_start[i] - 0.01, 
#                                               by = 0.01) + 0.01 + coord_offset,
#                                    ystart = edge_df$time_start[i],
#                                    yend = edge_df$time_end[i])
#   
#   #coord_offset <- coord_offset + edge_df$genome_end[i] - edge_df$genome_start[i] 
#   coord_offset <- sum(edge_df$genome_end[1:i] - edge_df$genome_start[1:i])
# }
# 
# 
# edge_viz_plot <- edge_vis_list %>% 
#   bind_rows() %>% 
#   ggplot() +
#   geom_rect(aes(xmin = xstartcoord, xmax = xend,
#                 ymin = ystart, ymax = yend, fill = xstart),
#             alpha = 1) +
#   scale_fill_gradientn(name = 'Genome\nposition',
#                        colours = c("#efedfb", "#dfdcf8", "#cfcaf4", "#afa7ed", "#8f84e6", "#6f61df", "#5648c6", "#43389a", "#262058"),
#                        limits = c(0, 4),
#                        breaks = c(0, 4),
#                        labels = c('Start', 'End')) +
#   #scale_fill_gradient(low = "#efedfb", high = "#393084", 
#   #                    limits = c(01, 4)) +
#   theme_void() +
#   #theme(legend.position = c(0.82, 0.25))
#   theme(legend.position = c(0.25, 0.77),
#         legend.background = element_rect(fill="white",
#                                          linetype = "solid", 
#                                          colour ="white"),
#         legend.direction = "horizontal",
#         legend.title = element_text(size = 16, hjust = 0.5),
#         legend.text = element_text(size = 12)
#   ) +
#   guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
# 
# 
# ggsave(filename = here('figures', 'pdf', 'edge_viz_plot_fig1.png'),
#        plot = edge_viz_plot,
#        bg = 'transparent',
#        width = 4,
#        height = 8,
#        units = 'in',
#        device = "png")


# fragment_size <- data.frame(
#   time = c(1, 1, 1, 1,
#            2, 2, 2, 2, 2,
#            3, 3, 3, 3,
#            4, 4, 4,
#            5, 5, 
#            6, 6, 6,
#            7, 7, 
#            8),
#   fragment_size = c(3, 3, 3, 3,
#                     3, 3, 1, 2, 3,
#                     3, 3, 2, 3,
#                     3, 3, 3,
#                     3, 3, 
#                     2, 1, 3,
#                     2, 3, 
#                     3)
# )

#edge_df <- data.frame(edge = as.character(1:11),
#                      genome_start = c(1, 1, 1, 3, 1, 1, 1, 1, 1, 2, 1),
#                      genome_end =   c(4, 3, 4, 4, 4, 4, 4, 2, 4, 4, 4),
#                      time_start =   c(0, 4, 6, 4, 4, 2, 0, 0, 2, 0, 0),
#                      time_end =     c(4, 7, 7, 6, 6, 4, 2, 2, 4, 3, 3))


# node_count_plot <- arg_attribute_df %>% 
#   ggplot() +
#   geom_segment(aes(x = time, xend = time, 
#                    y = shared_node_count, yend = local_tree_node_count),
#                size = 2) +
#   geom_point(aes(y = local_tree_node_count,
#                  x = time),
#              size = 9, color = '#CA4B9B') +
#   geom_point(aes(y = shared_node_count,
#                  x = time),
#              size = 9,
#              color = '#2ACAEA') +
#   ylab('Number of nodes') +
#   xlab('Time') +
#   coord_flip() +
#   theme_classic() +
#   theme(panel.background = element_rect(fill ='transparent'),
#         plot.background = element_rect(fill ='transparent', color = NA),
#         axis.title.x = element_text(size = 17),
#         axis.text.x = element_text(size = 14),
#         #panel.grid.major.x = element_line(linetype = 'solid'),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank()#,
#         #axis.line.y = element_blank(),
#         #axis.ticks.y = element_blank()
#   ) +
#   scale_y_continuous(breaks = seq(1, 12, 1))
# 
# ggsave(filename = '/Users/alexlewanski/Documents/michigan_state/research/arg_review/figures/pdf/node_count_fig1.png',
#        plot = node_count_plot,
#        bg = 'transparent',
#        width = 3.5,
#        height = 8,
#        units = 'in',
#        device = "png")
# 

# edge_df %>% 
#   ggplot() +
#   geom_rect(aes(xmin = genome_start, 
#                xmax = genome_end, 
#                ymin = time_start, 
#                ymax = time_end),
#             fill = 'gray',
#             alpha = 0.2) +
#   theme_classic()

# fragment_size %>% 
#   ggplot() +
#   #geom_line(aes(y = fragment_size,
#   #              x = time),
#   #          size = 2) +
#   geom_point(aes(y = fragment_size,
#                  x = time),
#              color = 'gray',
#              size = 8) +
#   coord_flip() +
#   ylab('Number of lineages') +
#   xlab('Time') +
#   theme_classic() +
#   theme(panel.background = element_rect(fill ='transparent'),
#         plot.background = element_rect(fill ='transparent', color = NA),
#         axis.title.x = element_text(size = 17),
#         axis.text.x = element_text(size = 14),
#         #panel.grid.major.x = element_line(linetype = 'solid', size = 1),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank()#,
#         #axis.line.y = element_blank(),
#         #axis.ticks.y = element_blank()
#   )
