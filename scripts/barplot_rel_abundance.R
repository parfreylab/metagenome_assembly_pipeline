# A simple barplot for relative abundances.
library(optparse)

# read in normalized count data
options.list <- list(make_option(c("-f",
                                   "--file_to_plot"),
                                   help="Path to input count data to plot, string."),
                     make_option(c("-t",
                                    "--title_for_graph"),
                                 default="TITLE",
                                 help="Title for barplot graph, string."),
                     make_option(c("-x",
                                    "--x_graph_label"),
                                 default="X AXIS",
                                 help="X axis label for barplot, string."),
                     make_option(c("-y",
                                   "--y_graph_label"),
                                 default="Y AXIS",
                                 help="Y axis label for barplot, string."),
                     make_option(c("-o",
                                   "--output_path"),
                                 default="./",
                                 help="Path to output folder, string. Default: current working directory."))

opt.parser <- OptionParser(option_list = options.list)
opt <- parse_args(opt.parser)

path <- opt$file_to_plot
graph.title <- opt$title_for_graph
graph.x.axis <- opt$x_graph_label
graph.y.axis <- opt$y_graph_label
out <- opt$output_path

if (is.null(path)) {
  stop("Path to data for plotting must be specified. See script usage (--help)")
}
if (!endsWith(out, "/")) {
  out <- paste(out, "/", sep = "")
}

to.plot <- read.delim(path, sep = "\t", header = T)

# reformat the sample names (collapse a little bit) so they look nicer on the graph
names(to.plot) <- gsub("quadrat", "qd", names(to.plot))
names(to.plot) <- gsub("water", "w", names(to.plot))
names(to.plot) <- gsub("pre_treatment", "pt", names(to.plot))
# double check the names before proceeding

# ALL THIS STUFF IS ACCORDING TO MELISSA'S BARPLOT TUTORIAL
rownames(to.plot) <- to.plot[,1]
to.plot <- to.plot[, -1]
to.plot <- data.matrix(to.plot)

numcolor <- nrow(to.plot)

randocolors <- sample(colors()[-grep("gr(e|a)y|white", colors())], numcolor)

quartz(graph.title, 15,10)
par(fig = c(0, 0.65, 0, 1))
barplot(to.plot,
        col = randocolors[factor(rownames(to.plot))],
        las = 2,
        ylab = graph.y.axis,
        xlab = graph.x.axis,
        main = graph.title,
        cex.names = 0.5
)
par(fig = c(0.55,1,0,1),
    new = TRUE)
plot(0,0 # MAKE A BLANK PLOT
     , pch = '' # No points
     , xaxt = 'n' # Supress x axis
     , yaxt = 'n' # Supress y axis
     , bty = 'n' # Supress box around plot
     , xlab = '' # no xlable
     , ylab = '' # no ylable
)
legend('center' # puts in center of space
       , legend = levels(factor(rownames(to.plot))) # OTUs - levels will force legend to use same ordering at barplot()
       , pch = 21
       , col = 'black' # outline is black for pch 19
       , pt.bg = randocolors # fills in pch 19 with color
       , pt.cex = 2
       , bty = 'n' # I personally don't like a box around my legend
       #, cex = 0.8
       #, ncol = 2
)
quartz.save(file = paste(out, graph.title, ".png", sep = ""), type = "png")
