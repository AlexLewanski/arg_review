######################################################
######################################################
### PROCESSING AND VISUALIZING MSPRIME SIMULATIONS ###
######################################################
######################################################

####################
### SCRIPT SETUP ###
####################

### packages ###
library(here)
library(dplyr)
library(ggplot2)
library(ape)
library(TreeDist)
library(ggforce)
library(egg)
library(ggtree)
library(ggpubr)
library(ggridges)
library(dichromat)
library(cowplot)


### custom functions ###

interval_midpt <- function(interval_string) {
  number_list <- lapply(strsplit(gsub('\\[([0-9e+\\.]*),([0-9e+\\.]*).*$', '\\1 \\2', interval_string), split = " "), function(x) as.numeric(x))
  return(sapply(number_list, mean))
}

tree_dist_wrapper <- function(tree_list, include_self_comparison = FALSE) {
  output_df <- as.data.frame(t(combn(seq_along(tree_list), 2))) %>% 
    rename(tree1 = V1,
           tree2 = V2)
  
  if (isTRUE(include_self_comparison)) {
    output_df <- rbind(output_df, 
                       data.frame(tree1 = seq_along(tree_list),
                                  tree2 = seq_along(tree_list)))
  }
  
  output_df$distance <- NA
  
  prog_bar <- txtProgressBar(min = 1, max = nrow(output_df), style = 3, char = "*")
  
  for (COMPARISON_INDEX in seq_len(nrow(output_df))) {
    #output_df[COMPARISON_INDEX,'distance'] <- TreeDistance(tree_list[[output_df[COMPARISON_INDEX,'tree1']]], 
    #                                                       tree_list[[output_df[COMPARISON_INDEX,'tree2']]])
    output_df[COMPARISON_INDEX,'distance'] <- InfoRobinsonFoulds(
      tree1 = tree_list[[output_df[COMPARISON_INDEX,'tree1']]],
      tree2 = tree_list[[output_df[COMPARISON_INDEX,'tree2']]],
      similarity = TRUE,
      normalize = FALSE,
      reportMatching = FALSE
    )
    
    setTxtProgressBar(prog_bar, COMPARISON_INDEX)
  }
  
  return(output_df)
}

spr_wrapper <- function(tree_list, include_self_comparison = FALSE) {
  output_df <- as.data.frame(t(combn(seq_along(tree_list), 2))) %>% 
    rename(tree1 = V1,
           tree2 = V2)
  
  if (isTRUE(include_self_comparison)) {
    output_df <- rbind(output_df, 
                       data.frame(tree1 = seq_along(tree_list),
                                  tree2 = seq_along(tree_list)))
  }
  
  output_df$distance <- NA
  
  prog_bar <- txtProgressBar(min = 1, max = nrow(output_df), style = 3, char = "*")
  
  for (COMPARISON_INDEX in seq_len(nrow(output_df))) {
    #output_df[COMPARISON_INDEX,'distance'] <- TreeDistance(tree_list[[output_df[COMPARISON_INDEX,'tree1']]], 
    #                                                       tree_list[[output_df[COMPARISON_INDEX,'tree2']]])
    output_df[COMPARISON_INDEX,'distance'] <- SPRDist(
      tree1 = tree_list[[output_df[COMPARISON_INDEX,'tree1']]],
      tree2 = tree_list[[output_df[COMPARISON_INDEX,'tree2']]],
      symmetric = TRUE
    )
    
    setTxtProgressBar(prog_bar, COMPARISON_INDEX)
  }
  
  return(output_df)
}


id_index_breakpoints <- function(x, output_type = c('index', 'val')) {
  
  store_list <- list()
  init <- 1
  for (i in seq_along(x)) {
    if (i == 1) next
    
    if (x[i] != (x[i - 1] + 1)) {
      if (output_type == 'val') {
        store_list <- append(store_list, list(c(x[init], x[i - 1])))
      } else {
        store_list <- append(store_list, list(c(init, i - 1)))
      }
      
      init <- i
    }
  }
  
  if (output_type == 'val') {
    return(append(store_list, list(c(x[init], x[length(x)]))))
  } else {
    return(append(store_list, list(c(init, length(x)))))
  }
}


### load data ###
tree_span_df <- read.csv(here('figures', 'sim_material', 'output', 'tree_span_df.csv.gz'))
tree_height_df <- read.csv(here('figures', 'sim_material', 'output', 'tree_height_df.csv.gz'))

node_composition_df <- read.csv(here('figures', 'sim_material', 'output', 'node_composition_df.csv.gz'))
node_membership_df <- read.csv(here('figures', 'sim_material', 'output', 'node_membership_df.csv.gz'))
node_info_full_arg_df <- read.csv(here('figures', 'sim_material', 'output', 'node_info_full_arg_df.csv.gz'))

main_treeseq_list <- read.nexus(file = here('figures', 'sim_material', 'output', 'sim_trees_review.nexus'))

### initial processing of some of the input ###
tree_info_merged <- tree_span_df %>% 
  left_join(., tree_height_df, by = 'tree_index') %>% 
  mutate(tree_span = right - left)

node_membership_df_wtreeinfo <- node_membership_df %>% 
  left_join(., tree_info_merged, 
            by = 'tree_index')


### VIZ MISC. ###
purple_vec <- c("#dfdcf8", "#cfcaf4", "#8f84e6", "#6f61df", "#5648c6", "#43389a")



###########################
### TREE HEIGHT AND AGE ###
###########################

#tree_indices <- c(46, 438, 439, 576)
tree_indices <- c(45, 437, 438, 575)

tree_height_plot <- tree_info_merged %>% 
  ggplot() +
  geom_step(aes(x = left, y = tree_height),
            size = 1.5,
            color = '#6C757D'
  ) +
  geom_point(data = . %>% 
               slice(tree_indices) %>% 
               #filter(tree_index %in% tree_indices) %>% 
               #group_by(tree_index) %>% 
               rowwise() %>% 
               mutate(mdpt = mean( c(left, right) )) %>% 
               ungroup(),
             aes(x = mdpt, 
                 y = tree_height, 
                 color = factor(tree_index, levels = tree_indices - 1)),
             size = 6) +
  scale_color_manual(values = c("#cfcaf4", "#8f84e6", "#5648c6", "#43389a")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(linewidth = 0.35, linetype = 'dashed', color = "#cfcaf4"),
        legend.position = 'none') +
  xlab("Genome position") +
  ylab('Tree height (TMRCA)')

ggsave(filename = here('figures', 'pdf', 'tree_height_plot.png'),
       plot = tree_height_plot,
       bg = 'transparent',
       width = 10,
       height = 2,
       units = 'in',
       device = "png")

### tree height summaries reported in main text ###
round(range(tree_info_merged$tree_height), 2)
round(mean(tree_info_merged$tree_height), 2)
round(sd(tree_info_merged$tree_height), 2)



############################
### VIZ OF EXAMPLE TREES ###
############################
tree_indices <- c(45, 437, 438, 575)

color_df <- data.frame(id = main_treeseq_list[[tree_indices[1]]]$tip.label,
                       index = as.character(1:length(main_treeseq_list[[tree_indices[1]]]$tip.label)) )


colfunc <- colorRampPalette(c("black", "white"))
gray_scale <- colfunc(24)
# color_vec <- c(gray_scale[16], gray_scale[3], gray_scale[4], gray_scale[1], gray_scale[21],
#                gray_scale[14], gray_scale[8], gray_scale[5], gray_scale[6],
#                gray_scale[23], gray_scale[19], gray_scale[17], gray_scale[20], gray_scale[11],
#                gray_scale[22],gray_scale[2], gray_scale[10], gray_scale[15], gray_scale[18],
#                gray_scale[7], gray_scale[9], gray_scale[12],gray_scale[24], gray_scale[13])
# 

# color_vec <- c(gray_scale[24], gray_scale[20], gray_scale[19], gray_scale[9], gray_scale[8],
#                gray_scale[7], gray_scale[6], gray_scale[5], gray_scale[4],
#                gray_scale[3], gray_scale[2], gray_scale[23], gray_scale[1], gray_scale[13],
#                gray_scale[12], gray_scale[11], gray_scale[10], gray_scale[22], gray_scale[21],
#                gray_scale[18], gray_scale[17], gray_scale[16], gray_scale[15], gray_scale[14])

color_vec <- c(gray_scale[24], gray_scale[20], gray_scale[19], gray_scale[13], gray_scale[12],
               gray_scale[11], gray_scale[10], gray_scale[9], gray_scale[8],
               gray_scale[7], gray_scale[6], gray_scale[23], gray_scale[5], gray_scale[4],
               gray_scale[3], gray_scale[2], gray_scale[1], gray_scale[22], gray_scale[21],
               gray_scale[18], gray_scale[17], gray_scale[16], gray_scale[15],gray_scale[14])


tree_scale_x <- 500
for (tree_ind in tree_indices) {
  
  #tree_plot <- ggtree(main_treeseq_list[[tree_ind - 1]]) %<+%
  tree_plot <- ggtree(main_treeseq_list[[tree_ind]]) %<+%
    color_df + 
    geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
    scale_fill_manual(values = color_vec) +
    layout_dendrogram() +
    theme(legend.position = 'none') +
    geom_treescale(x = tree_scale_x)
  
  ggsave(filename = here('figures', 'pdf', paste0('treeplot_', tree_ind,'.png')),
         plot = tree_plot,
         bg = 'transparent',
         width = 4,
         height = 4,
         units = 'in',
         device = "png")
  
}



#################################
### COMPARING TREE TOPOLOGIES ###
#################################

### CALCUlATE RF SIMILARITY AND SPR DISTANCE ###
rf_tree_sim_df <- tree_dist_wrapper(tree_list = main_treeseq_list, 
                                    include_self_comparison = TRUE)

spr_tree_sim_df <- spr_wrapper(tree_list = main_treeseq_list, 
                               include_self_comparison = TRUE)

### PLOTS FOR ROBINSON-FOULDS SIMILARITY ###
rf_tree_matrix <- rf_tree_sim_df %>% 
  rbind(.,
        rf_tree_sim_df %>% 
          select(tree2, tree1, distance) %>% 
          rename(tree1 = tree2, tree2 = tree1)) %>% 
  rename(`RF similarity` = distance) %>% 
  ggplot() +
  geom_tile(aes(x = tree1, y = tree2, fill = `RF similarity`)) +
  #theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = -3, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = -3, r = 0, b = 0, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white')) +
  xlab('Tree index') +
  ylab('Tree index') +
  scale_fill_gradientn(colors = purple_vec[-1])

index_vs_rf_plot <- rf_tree_sim_df %>% 
  mutate(tree_separation = abs(tree1 - tree2)) %>% 
  group_by(tree_separation) %>% 
  summarize(mean_size = mean(distance),
            lower0 = quantile(distance, prob = 0),
            lower10 = quantile(distance, prob = 0.1),
            lower20 = quantile(distance, prob = 0.20),
            lower30 = quantile(distance, prob = 0.20),
            lower40 = quantile(distance, prob = 0.20),
            mdpt = quantile(distance, prob = 0.50),
            upper60 = quantile(distance, prob = 0.60),
            upper70 = quantile(distance, prob = 0.70),
            upper80 = quantile(distance, prob = 0.80),
            upper90 = quantile(distance, prob = 0.90),
            upper100 = quantile(distance, prob = 1)) %>% 
  ggplot() +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower0, ymax = upper100),
              color = '#e1e3e5', fill = '#e1e3e5') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower10, ymax = upper90),
              color = '#caccce', fill = '#caccce') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower20, ymax = upper80),
              color = '#b4b5b7', fill = '#b4b5b7') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower30, ymax = upper70),
              color = '#989ea4', fill = '#989ea4') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower40, ymax = upper60),
              color = '#616970', fill = '#616970') +
  geom_line(aes(x = tree_separation, 
                y = mdpt), color = 'black', linewidth = 1.5, linetype = 'solid') +
  theme_classic() +
  theme(panel.grid.major.y = element_line(linewidth = 0.35, linetype = 'dashed', color = "#cfcaf4")) +
  xlab("Tree separation (count of intervening trees)") +
  ylab('Robinson–Foulds (RF) similarity')


### PLOTS FOR SPR DISTANCE ###
spr_tree_matrix <- spr_tree_sim_df %>% 
  rbind(.,
        spr_tree_sim_df %>% 
          select(tree2, tree1, distance) %>% 
          rename(tree1 = tree2, tree2 = tree1)) %>% 
  rename(`SPR distance` = distance) %>% 
  ggplot() +
  geom_tile(aes(x = tree1, y = tree2, fill = `SPR distance`)) +
  #theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = -3, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = -3, r = 0, b = 0, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white')) +
  xlab('Tree index') +
  ylab('Tree index') +
  scale_fill_gradientn(colors = rev(purple_vec),
                       breaks = c(0, 2, 4, 6, 8, 10), 
                       labels = as.character(c(0, 2, 4, 6, 8, 10)))

scale_color_continuous(breaks = c(100, 200, 300), labels = c("low", "med", "high"))
index_vs_spr_plot <- spr_tree_sim_df %>% 
  mutate(tree_separation = abs(tree1 - tree2)) %>% 
  group_by(tree_separation) %>% 
  summarize(mean_size = mean(distance),
            lower0 = quantile(distance, prob = 0),
            lower10 = quantile(distance, prob = 0.1),
            lower20 = quantile(distance, prob = 0.20),
            lower30 = quantile(distance, prob = 0.20),
            lower40 = quantile(distance, prob = 0.20),
            mdpt = quantile(distance, prob = 0.50),
            upper60 = quantile(distance, prob = 0.60),
            upper70 = quantile(distance, prob = 0.70),
            upper80 = quantile(distance, prob = 0.80),
            upper90 = quantile(distance, prob = 0.90),
            upper100 = quantile(distance, prob = 1)) %>% 
  ggplot() +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower0, ymax = upper100),
              color = '#e1e3e5', fill = '#e1e3e5') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower10, ymax = upper90),
              color = '#caccce', fill = '#caccce') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower20, ymax = upper80),
              color = '#b4b5b7', fill = '#b4b5b7') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower30, ymax = upper70),
              color = '#989ea4', fill = '#989ea4') +
  geom_ribbon(aes(x = tree_separation,
                  ymin = lower40, ymax = upper60),
              color = '#616970', fill = '#616970') +
  geom_line(aes(x = tree_separation, 
                y = mdpt), color = 'black', linewidth = 1.5, linetype = 'solid') +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), 
                     labels = as.character(c(0, 2, 4, 6, 8, 10))) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(linewidth = 0.35, linetype = 'dashed', color = "#cfcaf4")) +
  xlab("Tree separation (count of intervening trees)") +
  ylab('Subtree-prune-and-regraft (SPR) distance')


tree_rf_dist_multipanel <- cowplot::plot_grid(index_vs_rf_plot, rf_tree_matrix, nrow = 1)
tree_spr_dist_multipanel <- cowplot::plot_grid(index_vs_spr_plot, spr_tree_matrix, nrow = 1)

ggsave(filename = paste0(here('figures', 'pdf', 'tree_rf_dist_multipanel.png')),
       plot = tree_rf_dist_multipanel,
       bg = 'transparent',
       width = 10.5,
       height = 4,
       units = 'in',
       device = "png")

ggsave(filename = paste0(here('figures', 'pdf', 'tree_spr_dist_multipanel.png')),
       plot = tree_spr_dist_multipanel,
       bg = 'transparent',
       width = 10.5,
       height = 4,
       units = 'in',
       device = "png")



########################################
### ANCESTRY TRACTS AND NODE SHARING ###
########################################
merged_tract_list <- list()

for (SAMPLE in unique(node_membership_df_wtreeinfo$sample)) {
  
  focal_sample_df <- node_membership_df_wtreeinfo %>% 
    filter(sample == SAMPLE)
  
  for (NODE_INDEX in unique(focal_sample_df$node_index)) {
    focal_sample_df_nodeindex <- focal_sample_df %>% 
      filter(node_index == NODE_INDEX) %>% 
      arrange(tree_index)
    
    breakpoint_indices <- id_index_breakpoints(x = focal_sample_df_nodeindex$tree_index, 
                                               output_type = c('index', 'val')[1])
    
    for (i in breakpoint_indices) {
      focal_processing_df <- focal_sample_df_nodeindex[c(i[1]:i[2]),]
      
      subsetdf <- focal_processing_df[1,]
      subsetdf[,'right'] <- focal_processing_df[nrow(focal_processing_df),'right']
      
      merged_tract_list <- append(merged_tract_list, list(subsetdf))
    }
  }
}

merged_tract_df <- do.call(rbind, merged_tract_list)


node_composition_df$total_trees <- nrow(tree_info_merged)
node_composition_df_decompose <- split(node_composition_df, node_composition_df$node_index)

node_composition_df_decompose_update <- lapply(node_composition_df_decompose, function(x, tree_info) {
  x$available_trees <- sum(tree_info_merged$tree_height >= x$node_time[1])
  return(x)
}, tree_info = tree_info_merged) %>% 
  bind_rows()

node_composition_df_decompose_update1 <- node_composition_df_decompose_update %>% 
  mutate(time_interval = cut(node_time,
                             seq(0, 2000, by = 10), 
                             include.lowest = TRUE, 
                             right = FALSE)
  ) %>% 
  group_by(node_index) %>% 
  summarize(time_interval = first(time_interval),
            prop_total = n()/total_trees[1],
            prop_available = n()/available_trees[1],
            node_time = first(node_time)) %>% 
  ungroup() %>% 
  mutate(total_prop_quantile = cut(prop_total,
                                   seq(0, 1, by = 0.1), 
                                   include.lowest = TRUE, 
                                   right = FALSE),
         available_prop_quantile = cut(prop_available,
                                       seq(0, 1, by = 0.1), 
                                       include.lowest = TRUE, 
                                       right = FALSE)
  )


interval_summary_size <- 20

segment_summary <- merged_tract_df %>% 
  mutate(size = right - left) %>% 
  mutate(time_interval = cut(node_time,
                             seq(min(merged_tract_df$node_time), 
                                 ceiling(max(merged_tract_df$node_time)) + interval_summary_size, 
                                 by = interval_summary_size), 
                             include.lowest = TRUE, 
                             right = FALSE)
  ) %>% 
  group_by(time_interval) %>% 
  summarize(mean_size = mean(size),
            lower0 = quantile(size, prob = 0),
            lower5 = quantile(size, prob = 0.05),
            lower10 = quantile(size, prob = 0.1),
            lower25 = quantile(size, prob = 0.25),
            midpoint = quantile(size, prob = 0.50),
            upper90 = quantile(size, prob = 0.90),
            upper95 = quantile(size, prob = 0.95),
            upper75 = quantile(size, prob = 0.75),
            upper100 = quantile(size, prob = 1)) %>% 
  mutate(time_midpoint = interval_midpt(time_interval))

node_info_full_arg_df_sumstep <- node_info_full_arg_df %>%
  filter(flags == 131072) %>%
  group_by(time) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(time_interval_one = cut(time,
                                 seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
                                 include.lowest = TRUE,
                                 right = FALSE)
  ) %>%
  group_by(time_interval_one) %>%
  summarize(time_interval_sum = n(),
            .groups = 'drop')

recombo_df_cumsum <- data.frame(time_interval_one = cut(seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
                                   seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
                                   include.lowest = TRUE,
                                   right = FALSE),
           time_interval_sum = 0) %>%
  rbind(., node_info_full_arg_df_sumstep) %>%
  mutate(time_interval_one = factor(time_interval_one, levels = levels(cut(seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
                                                                           seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
                                                                           include.lowest = TRUE,
                                                                           right = FALSE))
  )) %>%
  group_by(time_interval_one) %>%
  summarize(time_interval_sum = sum(time_interval_sum), .groups = 'drop') %>%
  mutate(recombo_cumsum = cumsum(time_interval_sum),
         lower_bound = as.numeric(gsub('\\[([0-9e+\\.]*),([0-9e+\\.]*).*$', '\\1', time_interval_one)),
         upper_bound = as.numeric(gsub('\\[([0-9e+\\.]*),([0-9e+\\.]*).*$', '\\2', time_interval_one)))




tract_vis_plot <- merged_tract_df %>% 
  ggplot() +
  geom_segment(aes(x = left, xend = right, 
                   y = node_time, yend = node_time,
                   color = as.character(sample)),
               size = 1.5, alpha = 1) +
  theme_classic() +
  theme(plot.margin = margin(5.5, 0, 5.5, 5.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.grid.major.y = element_line(linewidth = 0.45, linetype = 'dashed', color = "#cfcaf4")) +
  scale_color_manual(values = c('#DEE2E6', '#6C757D', '#343A40')) +
  xlab('Genome position') +
  ylab('Time since the present') +
  ylim(0, 1225)

segment_size_summary_plot <- segment_summary %>%
  ggplot() +
  geom_segment(aes(x = lower25, xend = upper75, y = time_midpoint, yend = time_midpoint), 
               color = 'gray', size = 1.5) +
  geom_point(aes(x = midpoint, y = time_midpoint), 
             shape = 21, color = '#616970', fill = 'white', stroke = 1, size = 2) +
  theme_classic() +
  theme(plot.margin = margin(5.5, 0, 5.5, 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.45, linetype = 'dashed', color = "#cfcaf4")) +
  xlab('Tract size') +
  ylim(0, 1225)

recombo_cumsum_plot <- recombo_df_cumsum %>% 
  ggplot() +
  geom_rect(aes(xmin = 0, xmax = recombo_cumsum,
                ymin = lower_bound, ymax = upper_bound, fill = recombo_cumsum),
            alpha = 1) +
  scale_fill_gradientn(name = 'Genome\nposition',
                       colors = c('#e7e8ea', '#b4b5b7', '#5a5a5b'),
                       limits = c(0, max(recombo_df_cumsum$recombo_cumsum)),
                       breaks = c(0, max(recombo_df_cumsum$recombo_cumsum)),
                       labels = c('Start', 'End')) +
  theme_classic() +
  theme(plot.margin = margin(5.5, 0, 5.5, 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.45, linetype = 'dashed', color = "#cfcaf4"),
        legend.position = 'none') +
  xlab('Cumsum recombination') +
  ylim(0, 1225)


node_sharing_plot <- node_composition_df_decompose_update1 %>% 
  ggplot() +
  geom_point(aes(x = prop_total, y = node_time),
             color= "#8f84e6", alpha = 0.5, size = 2.5) +
  theme_classic() +
  theme(plot.margin = margin(5.5, -1, 5.5, 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.line.y = element_blank(),
        #axis.text.x = element_text(color = 'white'),
        panel.grid.major.y = element_line(linewidth = 0.45, linetype = 'dashed', color = "#cfcaf4"),
        legend.key.height= unit(0.35, 'cm'),
        legend.key.width= unit(0.35, 'cm'),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9)
  ) +
  ylim(0, 1225) +
  xlab('Proportion of trees with node') +
  scale_x_continuous(labels = c('0', '0.25', '0.50', '0.75', '1'))



density_plot_node_time <- node_composition_df_decompose_update1 %>% 
  ggplot() +
  geom_density(aes(y = node_time), 
               fill = "#8f84e6",
               color = "#8f84e6") +
  #theme_void() +
  theme(
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = 'transparent'),
    axis.title.x = element_text(color = 'transparent'),
    axis.ticks = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, -2)
    ) +
  ylim(0, 1225)
  
node_sharing_plot_with_marg <- ggarrange(node_sharing_plot,
                                         density_plot_node_time, 
                                         widths = c(0.8, 0.2))

segment_multipanel <- ggarrange(tract_vis_plot, 
                                segment_size_summary_plot, 
                                recombo_cumsum_plot, 
                                node_sharing_plot_with_marg,
                                ncol = 4, 
                                widths = c(0.5, 0.2, 0.15, 0.25))


ggsave(filename = here('figures', 'pdf', 'segment_multipanel.png'),
       plot = segment_multipanel,
       bg = 'white',
       width = 12,
       height = 4,
       units = 'in',
       device = "png")



###############################
### MULTI-SIMS: GENE FLOW #####
###############################

mig_tree_height_df <- read.csv(here('figures', 'sim_material', 'output', 'mig_tree_height_df_combined.csv.gz'))
mig_tree_span_df <- read.csv(here('figures', 'sim_material', 'output', 'mig_tree_span_combined.csv.gz'))

mig_tree_count_df <- mig_tree_height_df %>% 
  mutate(#mig_character = as.character(mig_rate),
         #sim_index_character = as.character(sim_index),
         mig_ind = paste0(mig_rate, '_', sim_index)) %>% 
  group_by(mig_ind) %>% 
  summarize(mig = first(mig_rate),
            tree_count = n(), .groups = 'drop')

mig_tree_span_count_combined <- mig_tree_span_df %>% 
  mutate(span = right - left,
         mig_factor = factor(mig_rate),
         mig_ind = paste0(mig_rate, '_', sim_index)) %>% 
  group_by(mig_ind) %>% 
  summarize(mean_span = mean(span)) %>% 
  left_join(., mig_tree_count_df, 
            by = 'mig_ind')

mig_tree_height_plot <- mig_tree_height_df %>% 
  mutate(mig_character = as.character(mig_rate),
         sim_index_character = as.character(sim_index)) %>%
  #group_by(mig_character, sim_index_character) %>% 
  filter(tree_height < quantile(tree_height, prob = 0.95)) %>% 
  mutate(mig = factor(mig_rate, levels = sort(unique(mig_tree_height_df$mig_rate)))) %>% 
  ggplot(aes(x = tree_height, y = mig)) +
  geom_density_ridges(aes(fill = sim_index_character), 
                      scale = 0.6, 
                      alpha = 0.4, 
                      size = 0.2, 
                      rel_min_height = 0.01,
                      color = "#8f84e6",
                      #fill = '#989ea4'
  ) + 
  scale_fill_manual(values = rep('#989ea4', 30)) +
  geom_point(data = . %>% 
               group_by(mig_character, sim_index_character) %>% 
               summarize(mig = first(mig),
                         tree_height = mean(tree_height)),
             aes(x = tree_height, y = mig),
             position=position_jitter(width = 0, height = 0.05),
             size = 1.2, alpha = 0.6, shape = 21,
             color = '#363a3e', stroke = 0.1,
             fill = "#8f84e6") +
  theme_ridges() +
  #theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 0, 5, 5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4")) +
  xlab('Tree height (TMRCA of local trees)') +
  ylab('Migration rate')

mig_tree_count_plot <- mig_tree_span_count_combined %>% 
  ggplot() +
  geom_sina(aes(y = factor(mig), x = tree_count, fill = mean_span),
            shape = 21, color = '#363a3e', size = 2, stroke = 0.1) +
  scale_fill_gradientn(name = 'Mean non-recomb.\nregion size',
                       colors = c('#e1e3e5','#989ea4', '#616970', '#363a3e', '#2d2d2d', '#000000')[-5]) +
  theme_ridges() +
  theme(axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        plot.margin = margin(5, 5, 5, 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4"),
        legend.key.height= unit(0.35, 'cm'),
        legend.key.width= unit(0.35, 'cm'),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9)) +
  xlab('Local tree count')

mig_sim_multipanel <- egg::ggarrange(mig_tree_height_plot, mig_tree_count_plot, ncol = 2,
                                             widths = c(0.7, 0.3))

mig_multipanel_fig <- annotate_figure(mig_sim_multipanel, 
                                              top = text_grob("Gene flow",
                                                              color = "Black",
                                                              face = "plain",
                                                              size = 18))

ggsave(filename = here('figures', 'pdf', 'mig_multipanel_fig.png'),
       plot = mig_multipanel_fig,
       bg = 'transparent',
       width = 8.25,
       height = 5,
       units = 'in',
       device = "png")



###############################
### MULTI-SIMS: SAMPLE SIZE ###
###############################

sample_size_tree_height_df <- read.csv(here('figures', 'sim_material', 'output', 'samp_size_tree_height_combined.csv.gz'))
sample_size_tree_span_df <- read.csv(here('figures', 'sim_material', 'output', 'samp_size_tree_span_combined.csv.gz'))

sample_size_tree_count_df <- sample_size_tree_height_df %>% 
  mutate(sample_size_character = as.character(samp_size),
         sim_index_character = as.character(sim_index)) %>% 
  group_by(sample_size_character, sim_index_character) %>% 
  summarize(sample_size = first(samp_size),
            tree_count = n(), .groups = 'drop')

sample_size_tree_span_count_combined <- sample_size_tree_span_df %>% 
  mutate(span = right - left,
         sample_size_factor = factor(samp_size)) %>% 
  group_by(sample_size_factor) %>% 
  summarize(mean_span = mean(span)) %>% 
  left_join(., sample_size_tree_count_df %>% mutate(sample_size_factor = factor(sample_size)), 
            by = 'sample_size_factor')

sample_size_tree_height_plot <- sample_size_tree_height_df %>% 
  mutate(sample_size_character = as.character(samp_size),
         sim_index_character = as.character(sim_index)) %>%
  #group_by(sample_size_character, sim_index_character) %>% 
  filter(tree_height < quantile(tree_height, prob = 0.95)) %>% 
  mutate(sample_size = factor(samp_size, levels = sort(unique(sample_size_tree_height_df$samp_size)))) %>% 
  ggplot(aes(x = tree_height, y = sample_size)) +
  geom_density_ridges(aes(fill = sim_index_character), 
                      scale = 0.6, 
                      alpha = 0.4, 
                      size = 0.2, 
                      rel_min_height = 0.01,
                      color = "#8f84e6",
                      #fill = '#989ea4'
                      ) + 
  scale_fill_manual(values = rep('#989ea4', 30)) +
  geom_point(data = . %>% 
               group_by(sample_size_character, sim_index_character) %>% 
               summarize(sample_size = first(sample_size),
                         tree_height = mean(tree_height)),
             aes(x = tree_height, y = sample_size),
             position=position_jitter(width = 0, height = 0.05),
             size = 1.2, alpha = 0.6, shape = 21,
             color = '#363a3e', stroke = 0.1,
             fill = "#8f84e6") +
  theme_ridges() +
  #theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 0, 5, 5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4")) +
  xlab('Tree height (TMRCA of local trees)') +
  ylab('Sample size')

sample_size_tree_count_plot <- sample_size_tree_span_count_combined %>% 
  ggplot() +
  geom_sina(aes(y = factor(sample_size), x = tree_count, fill = mean_span),
            shape = 21, color = '#363a3e', size = 2, stroke = 0.1) +
  scale_fill_gradientn(name = 'Mean non-recomb.\nregion size',
                       colors = c('#e1e3e5','#989ea4', '#616970', '#363a3e', '#2d2d2d', '#000000')) +
  theme_ridges() +
  theme(axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        plot.margin = margin(5, 5, 5, 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4"),
        legend.key.height= unit(0.35, 'cm'),
        legend.key.width= unit(0.35, 'cm'),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9)) +
  xlab('Local tree count')

sample_size_sim_multipanel <- egg::ggarrange(sample_size_tree_height_plot, sample_size_tree_count_plot, ncol = 2,
                                         widths = c(0.7, 0.3))

sample_size_multipanel_fig <- annotate_figure(sample_size_sim_multipanel, 
                                              top = text_grob("Sample size",
                                                              color = "Black",
                                                              face = "plain",
                                                              size = 18))

ggsave(filename = here('figures', 'pdf', 'sample_size_multipanel_fig.png'),
       plot = sample_size_multipanel_fig,
       bg = 'transparent',
       width = 8.25,
       height = 5,
       units = 'in',
       device = "png")



############################
### MULTI-SIMS: POP SIZE ###
############################

pop_size_tree_height_df <- read.csv(here('figures', 'sim_material', 'output', 'pop_size_tree_height_df_combined.csv.gz'))
pop_size_tree_span_df <- read.csv(here('figures', 'sim_material', 'output', 'pop_size_tree_span_combined.csv.gz'))

pop_size_tree_count_df <- pop_size_tree_height_df %>% 
  mutate(#pop_size_character = as.character(pop_size),
         #sim_index_character = as.character(sim_index),
         pop_size_ind = paste0(pop_size, '_', sim_index)) %>% 
  group_by(pop_size_ind) %>% 
  summarize(pop_size = first(pop_size),
                             tree_count = n(), .groups = 'drop')

pop_size_tree_span_count_combined <- pop_size_tree_span_df %>% 
  mutate(span = right - left,
         pop_size_factor = factor(pop_size),
         pop_size_ind = paste0(pop_size, '_', sim_index)) %>% 
  group_by(pop_size_ind) %>% 
  summarize(mean_span = mean(span)) %>% 
  left_join(., pop_size_tree_count_df, 
            by = 'pop_size_ind')

pop_size_tree_height_plot <- pop_size_tree_height_df %>% 
  mutate(pop_size_character = as.character(pop_size),
         sim_index_character = as.character(sim_index)) %>%
  #group_by(pop_size_character, sim_index_character) %>% 
  filter(tree_height < quantile(tree_height, prob = 0.95)) %>% 
  mutate(pop_size = factor(pop_size, levels = sort(unique(pop_size_tree_height_df$pop_size)))) %>% 
  ggplot(aes(x = tree_height, y = pop_size)) +
  geom_density_ridges(aes(height = after_stat(density),
                          fill = sim_index_character), scale = 2,
                      size = 0.2, 
                      rel_min_height = 0.01,
                      alpha = 0.4,
                      color = "#8f84e6",
                      #fill = '#989ea4'#, 
                      #panel_scaling = FALSE
                      ) + 
  scale_fill_manual(values = rep('#989ea4', 30)) +
  geom_point(data = . %>% 
               group_by(pop_size_character, sim_index_character) %>% 
               summarize(pop_size = first(pop_size),
                         tree_height = mean(tree_height)),
             aes(x = tree_height, y = pop_size, fill = sim_index_character),
             position=position_jitter(width = 0, height = 0.05),
             size = 1.2, shape = 21,
             color = '#363a3e', stroke = 0.1,
             fill = "#8f84e6") +
  theme_ridges() +
  theme(legend.position = 'none',
        axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 0, 5, 5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4")) +
  xlab('Tree height (TMRCA of local trees)') +
  ylab('Population size')


pop_size_tree_count_plot <- pop_size_tree_span_count_combined %>% 
  ggplot() +
  geom_sina(aes(y = factor(pop_size), x = tree_count, fill = mean_span),
            shape = 21, color = '#363a3e', size = 2, stroke = 0.1) +
  scale_fill_gradientn(name = 'Mean non-recomb.\nregion size',
                       colors = c('#e1e3e5','#989ea4', '#616970', '#363a3e', '#2d2d2d', '#000000')) +
  #theme_classic() +
  theme_ridges() +
  theme(axis.title.x = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.7, color = 'black'),
        plot.margin = margin(5, 5, 5, 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.35,
                                          linetype = 'dashed',
                                          color = "#cfcaf4"),
        legend.key.height= unit(0.35, 'cm'),
        legend.key.width= unit(0.35, 'cm'),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9)) +
  xlab('Local tree count')

popsize_sim_multipanel <- egg::ggarrange(pop_size_tree_height_plot, pop_size_tree_count_plot, ncol = 2,
          widths = c(0.7, 0.35))

popsize_multipanel_fig <- annotate_figure(popsize_sim_multipanel, top = text_grob("Population size", 
                                      color = "Black", face = "plain", size = 18))

ggsave(filename = here('figures', 'pdf', 'popsize_multipanel_fig.png'),
       plot = popsize_multipanel_fig,
       bg = 'transparent',
       width = 8,
       height = 5,
       units = 'in',
       device = "png")



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

### VERSION OF THE NODE SHARING PLOT THAT IS A SINA PLOT ###
# node_sharing_plot <- node_composition_df_decompose_update1 %>%
#   mutate(group = 'Nodes') %>% 
#   arrange(prop_total) %>% #sort so that the highest proportion are plotted on top (https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2)
#   ggplot() +
#   geom_sina(aes(x = group, y = node_time, color = prop_total), 
#             size = 2.5, alpha = 0.5, jitter_y = FALSE) +
#   scale_color_gradientn(name="Proportion\nshared",
#                         colours = purple_vec, #c('#e1e3e5', '#989ea4', '#616970', '#363a3e', '#202325'),
#                         limits = c(0,1)) +
#   theme_classic() +
#   theme(plot.margin = margin(5.5, 5.5, 5.5, 1),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         #axis.ticks.y = element_blank(),
#         #axis.line.y = element_blank(),
#         axis.text.x = element_text(color = 'white'),
#         panel.grid.major.y = element_line(linewidth = 0.45, linetype = 'dashed', color = "#cfcaf4"),
#         legend.key.height= unit(0.35, 'cm'),
#         legend.key.width= unit(0.35, 'cm'),
#         legend.text = element_text(size = 9),
#         legend.title = element_text(size = 9)
#   ) +
#   ylim(0, 1225) +
#   xlab('Nodes')


# tree_1 <- ggtree(main_treeseq_list[[tree_indices[1]]]) + 
#   geom_tippoint(color = "#FDAC4F", shape = 19, size = 3) #+ 
#   #layout_dendrogram()
# 
# 
# tree_2 <- ggtree(main_treeseq_list[[tree_indices[2]]]) + 
#   geom_tippoint(color = "#FDAC4F", shape = 19, size = 3) #+ 
#   #layout_dendrogram()
# 
# 
# tree_3 <- ggtree(main_treeseq_list[[tree_indices[3]]]) + 
#   geom_tippoint(color = "#FDAC4F", shape = 19, size = 3) #+ 
#   #layout_dendrogram()
# 
# tree_4 <- ggtree(main_treeseq_list[[tree_indices[4]]]) + 
#   geom_tippoint(color = "#FDAC4F", shape = 19, size = 3) #+ 
#   #layout_dendrogram()
# 
# color_df <- data.frame(id = main_treeseq_list[[tree_indices[4]]]$tip.label,
#                        index = as.character(1:length(main_treeseq_list[[tree_indices[4]]]$tip.label)) )





#ggtree(main_treeseq_list[[1]]) +
#  geom_tiplab()



# tree_scale_x <- 650
# 
# ggtree(main_treeseq_list[[tree_indices[1]]]) %<+%
#   color_df + 
#   geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
#   scale_fill_manual(values = color_vec) +
#   layout_dendrogram() +
#   theme(legend.position = 'none') +
#   geom_treescale(x = tree_scale_x)
# 
# ggtree(main_treeseq_list[[tree_indices[2]]]) %<+%
#   color_df + 
#   geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
#   scale_fill_manual(values = color_vec) +
#   layout_dendrogram() +
#   theme(legend.position = 'none') +
#   geom_treescale(x = tree_scale_x)
# 
# ggtree(main_treeseq_list[[tree_indices[3]]]) %<+%
#   color_df + 
#   geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
#   scale_fill_manual(values = color_vec) +
#   layout_dendrogram() +
#   theme(legend.position = 'none') +
#   geom_treescale() +
#   geom_treescale(x = tree_scale_x)
# 
# ggtree(main_treeseq_list[[tree_indices[4]]]) %<+%
#   color_df + 
#   geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
#   scale_fill_manual(values = color_vec) +
#   layout_dendrogram() +
#   theme(legend.position = 'none') +
#   geom_treescale() +
#   geom_treescale(x = tree_scale_x)
# 
# 
# 
# 
# 
# 
# ggtree(main_treeseq_list[tree_indices], size = 1.15) +
#   facet_wrap(~.id, ncol=2) %<+%
#   color_df + 
#   geom_tippoint(aes(fill = as.character(index)), color = 'black', shape = 21, size = 3) +
#   scale_fill_manual(values = color_vec)
#   #layout_dendrogram() +
#   
# 
# 
# 


# tip_vec <- paste0('n', c(4, 2, 5, 7, 1, 13, 11, 6, 18, 12, 10, 23, 22, 16, 0, 21, 8, 14, 19, 9, 17, 15, 3, 20))
# #tip_vec <- paste0('n', c(20, 17, 11, 18, 23, 3, 12, 10, 21, 8, 1, 0, 13, 19, 9, 16, 5, 4, 2, 7, 22, 15, 14, 6))
# 
# #paste0('n', c(6, 14, 15, 22, 9, 19, 7, 2, 4, 5, 16, 18, 11, 17, 20, 13, 0, 1, 8, 21, 10, 12, 3, 23))
# color_df <- data.frame(id =tip_vec,
#            index = 1:length(main_treeseq_list[[tree_indices[4]]]$tip.label))
# 
# 
# 
# 
# plot_list(tree_1, tree_2, tree_3, tree_4, tag_levels = 'A')
# 
# 
# ggtree(main_treeseq_list[tree_indices], size = 1.15) + 
#   geom_tiplab(size = 2, align = TRUE) +
#   #geom_tippoint(aes(color = ), shape = 19, size = 3) +
#   layout_dendrogram() +
#   facet_wrap(~.id, ncol=2)





#################################
### TREE HEIGHT VS. SPAN SIZE ###
#################################

# cor_output <- cor.test(tree_info_merged$tree_height, 
#                        tree_info_merged$tree_span, 
#                        method = 'spearman', exact = FALSE)
# 
# treespan_vs_age_plot <- tree_info_merged %>% 
#   ggplot(aes(x = tree_span, y = tree_height)) +
#   geom_point(color = '#616970', alpha = 0.7, size = 3) +
#   geom_smooth(method = 'loess', 
#               se = FALSE,
#               color = "#8f84e6",
#               size = 2) +
#   theme_classic() +
#   theme(panel.grid.major = element_line(size = 0.35, 
#                                         linetype = 'dashed', 
#                                         color = "#cfcaf4")) +
#   xlab('Tree span (# base pairs)') +
#   ylab('Tree height (TMRCA)') +
#   annotate("text",
#            x = 88, y = 1130,
#              color = '#616970',
#            label = paste0('rho = ', round(cor_output$estimate, 3), 
#                           '\nP = ', signif(cor_output$p.value, digits=3)),
#            size = 5)

#correlation in tree heights
# as.data.frame(gtools::permutations(n = 200, r = 2)) %>% 
#   rename(x = V1, y = V2) %>% 
#   mutate(similarity = runif(nrow(gtools::permutations(n = 200, r = 2)), 0, 1)) %>% 
#   ggplot() +
#   geom_tile(aes(x = x, y = y, fill = similarity)) +
#   theme_minimal()
# 
# tree_height_example <- data.frame(position = 1:200,
#            value = runif(200, 0, 2)) 
# 

#beeswarm/violin plot of tree count vs. sample size?
#average segment size?


# data.frame(pop_size = rep(0:10, each = 200),
#            segment_size = rnorm(length(rep(0:10, each = 200)), 50, 2)) %>% 
#   ggplot() +
#   geom_violin(aes(x = as.character(pop_size), y = segment_size), color = 'gray') +
#   geom_sina(aes(x = as.character(pop_size), y = segment_size), 
#                    size = 0.6, color = 'gray', alpha = 0.5) +
#   theme_classic()


#beeswarm/violin plot of tree height vs. population size?


#line plot of position in the genome vs. time
#show that ancestral segments get smaller further back in time


# edge_info <- data.frame(
#   time_start = c(0, 0, 0, 0, 3, 3, 5, 7),
#   time_end = c(2, 2, 3, 4, 4, 5, 10, 12),
#   genome_start = c(0, 0, 0, 0, 0, 6, 0, 3),
#   genome_end = c(9, 9, 9, 9, 6, 9, 6, 9)
# )
# 
# 
# coord_offset <- 0
# edge_vis_list <- list()
# for (i in 1:nrow(edge_info)) {
#   edge_vis_list[[i]] <- data.frame(xstart = seq(edge_info$genome_start[i], 
#                                                 edge_info$genome_end[i] - 0.01, 
#                                                 by = 0.01),
#                                    xstartcoord = seq(edge_info$genome_start[i] - edge_info$genome_start[i], 
#                                                      edge_info$genome_end[i] - edge_info$genome_start[i] - 0.01, 
#                                                      by = 0.01) + coord_offset,
#                                    xend = seq(edge_info$genome_start[i] - edge_info$genome_start[i], 
#                                               edge_info$genome_end[i] - edge_info$genome_start[i] - 0.01, 
#                                               by = 0.01) + 0.01 + coord_offset,
#                                    ystart = edge_info$time_start[i],
#                                    yend = edge_info$time_end[i])
#   
#   #coord_offset <- coord_offset + edge_info$genome_end[i] - edge_info$genome_start[i] 
#   coord_offset <- sum(edge_info$genome_end[1:i] - edge_info$genome_start[1:i])
# }

# edge_vis_list %>% 
#   bind_rows() %>% 
#   ggplot() +
#   geom_rect(aes(xmin = xstartcoord, xmax = xend,
#                 ymin = ystart, ymax = yend, fill = xstart)) +
#   scale_fill_gradient(low = "#efedfb", high = "#393084", 
#                       limits = c(0,10)) +
#   theme_minimal()
# 


# segment_df <- data.frame(genome_start = runif(1000, 0, 100),
#                          time = sample(1:100, 1000, replace = TRUE))
# segment_df$genome_end <- NA

# for (i in 1:nrow(segment_df)) {
#   genome_end_init <- segment_df$genome_start[i] + runif(1, 0, 50)
#   
#   if (genome_end_init > 100) genome_end_init <- 100
#   
#   segment_df$genome_end[i] <- genome_end_init
# }


# segment_summary <- segment_df %>% 
#   mutate(size = genome_end - genome_start) %>% 
#   mutate(time_interval = cut(node_time,
#                              seq(0, 2000, by = 10), 
#                              include.lowest = TRUE, 
#                              right = FALSE)
#   )
#   group_by(time) %>% 
#   summarize(mean_size = mean(size),
#             lower0 = quantile(size, prob = 0),
#             lower5 = quantile(size, prob = 0.05),
#             lower25 = quantile(size, prob = 0.25),
#             upper95 = quantile(size, prob = 0.95),
#             upper75 = quantile(size, prob = 0.75),
#             upper100 = quantile(size, prob = 1))
# 
# 

# segment_summary %>% 
#   ggplot() +
#   geom_segment(aes(x = lower0, xend = upper100, y = time, yend = time), color = 'gray') +
#   geom_point(aes(x = mean_size, y = time), shape = 21, color = 'black', fill = 'white', stroke = 1) +
#   theme_classic()


# segment_summary <- merged_tract_df %>% 
#   mutate(size = right - left) %>% 
#   group_by(time) %>% 
#   summarize(mean_size = mean(size),
#             lower0 = quantile(size, prob = 0),
#             lower5 = quantile(size, prob = 0.05),
#             lower25 = quantile(size, prob = 0.25),
#             upper95 = quantile(size, prob = 0.95),
#             upper75 = quantile(size, prob = 0.75),
#             upper100 = quantile(size, prob = 1))

# segment_df %>% 
#   ggplot() +
#   geom_segment(aes(x = genome_start, xend = genome_end, y = time, yend = time)) +
#   theme_classic()
# merged_tract_df %>% 
#   ggplot() +
#   geom_segment(aes(x = left, xend = right, 
#                    y = node_time, yend = node_time,
#                    color = as.character(sample)),
#                size = 1.25, alpha = 0.75) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# 

#c("#efedfb", "#dfdcf8", "#cfcaf4", "#afa7ed", "#8f84e6", "#6f61df", "#5648c6", "#43389a", "#262058")
# node_composition_df_decompose_update1 %>%
#   mutate(group = 'a') %>% 
#   ggplot() +
#   geom_sina(aes(x = group, y = node_time, color = prop_available), 
#             size = 2.5, alpha = 0.5, jitter_y = FALSE) +
#   scale_color_gradientn(colours = c("#dfdcf8", "#cfcaf4", "#8f84e6", "#6f61df", "#5648c6", "#43389a"),
#                       limits = c(0,1)) +
#   theme_classic()
#   

#for (i in unique(node_composition_df$node_index))) {
#  node_height <- unique(node_composition_df$node_time[node_composition_df$node_index == i])
#}


# node_info_full_arg_df <- read.csv(here('figures', 'sim_material', 'output', 'node_info_full_arg_df.csv'))
# 
# node_info_full_arg_df_sumstep <- node_info_full_arg_df %>% 
#   filter(flags == 131072) %>% 
#   group_by(time) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   mutate(time_interval_one = cut(time,
#                                  seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1), 
#                                  include.lowest = TRUE, 
#                                  right = FALSE)
#   ) %>% 
#   group_by(time_interval_one) %>% 
#   summarize(time_interval_sum = n(),
#             .groups = 'drop')
#   
# 
# recombo_df_cumsum <- data.frame(time_interval_one = cut(seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
#                                    seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1), 
#                                    include.lowest = TRUE, 
#                                    right = FALSE),
#            time_interval_sum = 0) %>% 
#   rbind(., node_info_full_arg_df_sumstep) %>% 
#   mutate(time_interval_one = factor(time_interval_one, levels = levels(cut(seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1),
#                                                                            seq(0, ceiling(max(node_info_full_arg_df$time)), by = 1), 
#                                                                            include.lowest = TRUE, 
#                                                                            right = FALSE))
#   )) %>% 
#   group_by(time_interval_one) %>% 
#   summarize(time_interval_sum = sum(time_interval_sum), .groups = 'drop') %>%  
#   mutate(recombo_cumsum = cumsum(time_interval_sum),
#          lower_bound = as.numeric(gsub('\\[([0-9e+\\.]*),([0-9e+\\.]*).*$', '\\1', time_interval_one)),
#          upper_bound = as.numeric(gsub('\\[([0-9e+\\.]*),([0-9e+\\.]*).*$', '\\2', time_interval_one)))
# 
# 
# recombo_df_cumsum %>% 
#   ggplot() +
#   geom_rect(aes(xmin = 0, xmax = recombo_cumsum,
#                 ymin = lower_bound, ymax = upper_bound, fill = recombo_cumsum),
#             alpha = 1) +
#   scale_fill_gradientn(name = 'Genome\nposition',
#                        colors = c('#e7e8ea', '#b4b5b7', '#5a5a5b'),
#                        limits = c(0, max(recombo_df_cumsum$recombo_cumsum)),
#                        breaks = c(0, max(recombo_df_cumsum$recombo_cumsum)),
#                        labels = c('Start', 'End')) +
#   theme_bw() +
#   xlab('Cumulative recombo. sum')

