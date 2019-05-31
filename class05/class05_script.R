help("read.table")


#practice line plot
weight_table = read.table('bimm143_05_rstats/weight_chart.txt', header = TRUE)

plot(weight_table$Age, weight_table$Weight, type = "o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby weight vs age")


# practice bar plot 

feature_count = read.table('bimm143_05_rstats/feature_counts.txt', sep="\t", header = TRUE)

par(mar=c(5,8,5,5))
barplot(feature_count$Count, horiz=TRUE, names.arg = feature_count$Feature, main="Number of Features in Mouse Genome", las=1)


# practice histogram
random <- c(rnorm(10000),rnorm(10000)+4)
hist(random, breaks=10)

# color exercise

female_male <- read.table('bimm143_05_rstats/male_female_counts.txt', sep="\t", header = TRUE)
colors = rainbow(nrow(female_male))
barplot(female_male$Count, names.arg= female_male$Sample, col=c("blue2","red2"))


# color scatter practice 

genes <- read.table('bimm143_05_rstats/up_down_expression.txt', header = TRUE)
nrow(genes)
table(genes$State)
plot(genes$Condition1, genes$Condition2, col = genes$State)
palette()
levels(genes$State)
palette(c("blue","gray","red"))


# methylation

meth <- read.table('bimm143_05_rstats/expression_methylation.txt', header=TRUE)
inds <- meth$expression > 0
dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                                      colramp = colorRampPalette(c("blue2",
                                                                   "green2",
                                                                   "red2",
                                                                   "yellow")) )
plot(meth$gene.meth[inds], meth$expression[inds], col = dcols.custom,pch=20)
