
# Function to save a ggplot object in multiple formats
save_plot <- function(plot, out_path, base_filename, width = 12, height = 8, dpi = 300) {
    # Save as PDF
    ggsave(glue("{out_path}{base_filename}.pdf"), plot = plot, device = "pdf", width = width, height = height)
    
    # Save as TIFF
    ggsave(glue("{out_path}{base_filename}.tiff"), plot = plot, device = "tiff", width = width, height = height, dpi = dpi, bg = "white")
    
    # Save as SVG
    ggsave(glue("{out_path}{base_filename}.svg"), plot = plot, device = "svg", width = width, height = height)
    
    # Save as JPEG
    ggsave(glue("{out_path}{base_filename}.jpeg"), plot = plot, device = "jpeg", width = width, height = height, dpi = dpi)
}

