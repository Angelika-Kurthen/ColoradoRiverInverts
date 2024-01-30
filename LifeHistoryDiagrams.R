# Ensure tikzDevice is installed and loaded
if (!require(tikzDevice)) install.packages("tikzDevice")
library(tikzDevice)

# Create a new TikZ picture
tikz("circles_arrows.tex", width = 5, height = 5, standAlone = TRUE)

# Begin the TikZ picture
cat("\\begin{tikzpicture}\n")

# Draw the circles
cat("\\node [circle, draw] (S1) at (0,0) {S1};\n")
cat("\\node [circle, draw] (S2) at (3,0) {S2};\n")
cat("\\node [circle, draw] (S3) at (1.5,3) {S3};\n")

# Close the TikZ device
dev.off()
cat("\\draw [->, bend right] (A) to node[midway, below] {$G_i$} (B);\n")
cat("\\draw [->, bend right] (B) to (C);\n")
cat("\\draw [->, bend right] (C) to (A);\n")
cat("\\draw [->, bend right] (A) to (A);\n")
cat("\\draw [->, bend right] (C) to (B);\n")
cat("\\end{tikzpicture}\n")

dev.off()


