root   <- "/home/ailtonpcf/vast/proj/02-compost-microbes/cache"
task   <- "116-hpf-transposases-curatedHMM-funannotate-v2"

apereira_save_plot <- function(plot, plot_dir, plot_name, width = 9, height = 7, dpi = 150, ...) {
  
  dir.create(path = plot_dir, recursive = T)
  
  ggsave(filename = paste0(paste(plot_dir, plot_name, sep = "/"), ".pdf"), 
         plot = plot, width = width, height = height, dpi = dpi, ...)
  
  ggsave(filename = paste0(paste(plot_dir, plot_name, sep = "/"), ".png"), 
         plot = plot, width = width, height = height, dpi = dpi, ...)
  
}
