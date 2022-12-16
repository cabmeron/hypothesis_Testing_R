

# Question 1

# creating complete pathway for data by concatenating two strings before loading with function load()


  x <- "/Users/Girlmccoy/Desktop/Desktop/stinkyMATICS/bioinformatics/"
  y <- "GSE9782(2).RData"
  
full.data.path.file <- paste0(x,y)

load(full.data.path.file)




# Question 2

# using log2() to acquire new matrix from original data
# finding indices of treated/control groups by using which()
# consolidate function to acquire mean difference that can be used within apply() to acquire values from entire matrix
# return values into a variable named LogFC


dataGSE9782.log2 <- log2(dataGSE9782)

indices.Treated <- which(group == "treated")
indices.Control <- which(group == "control")

myMeanDiff <- function(row.value)
{
  treated.value <- row.value[indices.Treated]
  control.value <- row.value[indices.Control]
  
  treated.mean <- mean(treated.value)
  control.mean <- mean(control.value)
  mean.diff <- treated.mean - control.mean
  
  return(mean.diff)
  
  
}

LogFC <- apply(dataGSE9782.log2,1,myMeanDiff)


# Question 3

# function to extract t.scores from matrix with t.test

my.t.score <- function(row.value)
{
  
  
  treated.value <- row.value[indices.Treated]
  control.value <- row.value[indices.Control]
  
  
  
  
  expression.t.score <-(t.test(treated.value,control.value)$statistic[["t"]])
  
  
  
  return(expression.t.score)
  
  
}

# function to extract p.values from matrix with t.test

my.p.value <- function(row.value)
{
  
  
  treated.value <- row.value[indices.Treated]
  control.value <- row.value[indices.Control]
  
  
  
  expression.p.values <-(t.test(treated.value,control.value)$p.value)
  
  
  
  return(expression.p.values)
  
  
}

# using functions inside of apply() to get statistics from matrix rows

# extracting individual values from t.test into variables
# putting gene id's into variable with rownames()
# organizing all acquired values to data frame

data.p.values <- apply(dataGSE9782,1,my.p.value)
data.t.values <- apply(dataGSE9782,1,my.t.score)
gene.ids <- rownames(dataGSE9782)

my.expression.data.frame <- data.frame(gene.ids = gene.ids,
                                       p.values = data.p.values,
                                       t.scores = data.t.values,
                                       logFC = LogFC
                                       
                                       
                                       
                                       
)

save(my.expression.data.frame, file = "DE-results.rds")

# Question 4

# gathering proper variable values for proper volcano plot
# plotting / saving volcano plot
#

negLog <- -log10(my.expression.data.frame$p.values)
colors <- ifelse(abs(LogFC) > 1 & data.p.values < 0.05, "red", "black" )
plot(my.expression.data.frame$logFC,negLog, col = colors)


pdf(file = "/Users/Girlmccoy/Desktop/Desktop/stinkyMATICS/bioinformatics/volcano.pdf" )
plot(my.expression.data.frame$logFC,negLog, col = colors)
dev.off()

# Question 5

my.wilcoxon <- function(row.value)
{
  treated.value <- row.value[indices.Treated]
  control.value <- row.value[indices.Control]
  
  wilcoxon.result <- wilcox.test(treated.value,control.value)$p.value
 
  
}

wilcox.test.values <- apply(dataGSE9782,1,my.wilcoxon)

gene.ids <- rownames(dataGSE9782)[data.p.values < 0.05 & wilcox.test.values < 0.05]


saveRDS(gene.ids, "wilcoxon-results.")

# Question 6

# plotting histogram of log2(data file)
# saving histogram as pdf

log2.expression.levels <- log2(dataGSE9782)

pdf(file = "/Users/Girlmccoy/Desktop/Desktop/stinkyMATICS/bioinformatics/hist.pdf")
hist(log2.expression.levels)
dev.off()

# Question 7

# acquiring index containing maximum value in LogFC variable
# using index to acquire all values in the row containing our maximum value
# splitting row values into two sets before boxplotting
# saving boxplot as pdf
gene.index <- which(abs(LogFC) == max(abs(LogFC)))

print(gene.index)
View(gene.index)
expression.values <- log2.expression.levels[gene.index,]


treated <- expression.values[which(group == "treated")]
control <- expression.values[which(group == "control")]


my.boxplot.input <- list(control,treated)
pdf(file = "/Users/Girlmccoy/Desktop/Desktop/stinkyMATICS/bioinformatics/boxplot.pdf")
boxplot(my.boxplot.input, col = "orange", names = c("Control", "Treated"), ylab = "LogFC")
dev.off()


