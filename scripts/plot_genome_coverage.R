
# script to plot coverage over the JEV genomes

library(ggplot2)
library(scales)

plot_coverage <- function(name, depth, min_depth, primer_bed, pdf_file, png_file) {
 
   # get depth of coverage information
  depth_df = read.csv(depth, sep="\t", header=F)
  colnames(depth_df) <- c("ref", "Position", "Depth")

  # get bed file to draw in primer binding regions
  # haven't put these on the plot yet
  bed_df = read.csv(primer_bed, sep="\t", header=F)
  colnames(bed_df) <- c("ref", "start", "end", "name", "Depth", "strand")

  p <- ggplot(depth_df, aes(x=Position, y=Depth, fill=ref)) +
	geom_area(aes(alpha=0.5)) +
	geom_hline(yintercept=strtoi(min_depth), linetype="dotted", colour="red") +
	# add segment to show where the JEV genome starts and ends
	geom_segment(x=0, xend=10941, y=0, yend=0, size=2) +
	geom_point(x=0, y=0, size=3) +
	geom_point(x=10941, y=0, size=3) +
	scale_fill_brewer(palette = "Dark2") +
	#scale_fill_manual(values=cols) +
	#theme(axis.title.y = element_blank(), legend.title = element_blank()) +
	scale_y_continuous(trans = "pseudo_log", labels = comma, expand = c(0,0), breaks=c(1, 10, 100, 1000,10000)) +
	scale_x_continuous(breaks = seq(0, 20000, by = 1000)) +
	#coord_flip() +
	#facet_grid(rows=vars(Reference_Name), scales="free") +
	labs(
		title = paste0(name, sep=" ", "genome coverage")
	) +
	#theme_minimal(base_family = "Roboto Condensed") +
	theme(
		plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
		plot.title = element_text(size = 12, face = "bold"),
		axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
		axis.text = element_text(size = 10),
		legend.position = "none",
		panel.grid.major.y = element_blank(),
		)
   
  ggsave(pdf_file, p, width=12, height=4)
  ggsave(png_file, width=12, height=4, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
name = args[1]
depth = args[2]
min_depth = args[3]
primer_bed = args[4]
pdf_file = args[5]
png_file = args[6]

plot_coverage(name, depth, min_depth, primer_bed, pdf_file, png_file)


