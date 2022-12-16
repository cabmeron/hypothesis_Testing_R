#PROJECT 1 master

#1) Retrieve data matrix and patient group information of the dataset GSE19804 from NIH GEO

library(GEOquery)

eset <- getGEO("GSE19804", filename = "GSE19804_series_matrix.txt.gz")

myData <- exprs(eset[1:1000,]) #take 1000 for running faster


my.p.Data <- pData(eset)

control <- rownames(my.p.Data[grep("Lung Normal",my.p.Data$title),])

cancer <- rownames(my.p.Data[grep("Lung Cancer",my.p.Data$title),])


 #Function for T.score
my.t.score <- function(row.value)
{
  
  
  cancer.value <- row.value[cancer]
  control.value <- row.value[control]
  
  
  
  expression.t.score <-(t.test(cancer.value,control.value)$statistic[["t"]])
  
  
  
  return(expression.t.score)
  
  
}


#Function for P-value

my.p.value <- function(row.value)
{
  
  
  cancer.value <- row.value[cancer]
  control.value <- row.value[control]
  
  
  
  expression.p.values <-(t.test(cancer.value,control.value)$p.value)
  
  
  
  return(expression.p.values)
  
  
}

#Function for Mean Difference

myMeanDiff <- function(row.value)
{
  cancer.value <- row.value[cancer]
  control.value <- row.value[control]
  
  cancer.mean <- mean(cancer.value)
  control.mean <- mean(control.value)
  mean.diff <- cancer.mean - control.mean
  
  return(mean.diff)
  
  
}

#Acquiring P-values
my.p.values <- apply(myData,1,my.p.value)

#Acquiring T-scores
my.t.scores <- apply(myData,1,my.t.score)

gene.ids <- rownames(myData)

logFC <- apply(myData,1,myMeanDiff)

#Creating of data frame containing acquired statistics

my.expression.data.frame <- data.frame(gene_ids = gene.ids,
                                       p.values = my.p.values,
                                       t.scores = my.t.scores,
                                       LogFC = logFC
)

saveRDS(my.expression.data.frame, file = "DE-results.rds")


#Organization and plotting of volcano plot
pdf("volcano.pdf")
negLog <- -log10(my.expression.data.frame$p.values)
colors <- ifelse(abs(logFC) > 1 & my.p.values < 0.05, "red", "black")
plot(my.expression.data.frame$LogFC,negLog, col = colors)
dev.off()

print(my.expression.data.frame$LogFC)
#Storing T-scores in new data frame

t.OBSERVED <- data.frame(t.SCORE = my.t.scores)
saveRDS(t.OBSERVED, file = "t-observed.rds")


#Start of calculating empirical P-Values for T-scores
t.NULL.DISTRIBUTION = NULL
for(i in 1:100)
{
  indices <- sample(1:ncol(myData),ncol(myData))
  control.indices1 <- indices[1:length(control)]
  cancer.indices2 <- indices[-control.indices1]
  
  t.test.results <- apply(myData,
                          MARGIN = 1,
                          FUN = function(x)
                          {
                            test.result <- t.test(x[cancer.indices2], x[control.indices1])
                            c(test.result$statistic)
                          }
                          
                          
                          
  )
  
  t.NULL.DISTRIBUTION = cbind(t.NULL.DISTRIBUTION, t.test.results)
  
  
}


#Permutation analysis
p.T <- NULL

for(i in 1:nrow(myData))
{
  
  
  
  
  t.Observed.i <- t.OBSERVED[i,]
  t.NULL.DIST.i <- t.NULL.DISTRIBUTION[i,]
  count <- sum((abs(t.Observed.i) < abs(t.NULL.DIST.i)))
  
  
  
  p.T[i] <- (count / 100) 
}

saveRDS(p.T, file = "p-empirical-tscore.rds")

#Start of permutation analysis using Euclidian distance

myMeanDiff2 <- function(row.value)
{
  cancer.value <- row.value[cancer]
  control.value <- row.value[control]
  
  cancer.mean <- mean(cancer.value)
  control.mean <- mean(control.value)
  mean.diff <- abs(cancer.mean - control.mean)
  
  return(mean.diff)
  
  
}


e.OBSERVED.test <- apply(myData, 1, myMeanDiff2)

e.OBSERVED <- data.frame(e.OBSERVED.test)


e.NULL.DISTRIBUTION = NULL
for(i in 1:100)
{
  indices <- sample(1:ncol(myData),ncol(myData))
  control.indices1 <- indices[1:length(control)]
  cancer.indices2 <- indices[-control.indices1]
  
  e.score.results <- apply(myData,
                           MARGIN = 1, 
                           FUN = function(x)
                           {
                             cancer.value <- x[cancer.indices2]
                             control.value <- x[control.indices1]
                             
                             cancer.mean <- mean(cancer.value)
                             control.mean <- mean(control.value)
                             mean.diff <- abs(cancer.mean - control.mean)
                             
                             return(mean.diff)
                           }
                           
  )
  
  e.NULL.DISTRIBUTION = cbind(e.NULL.DISTRIBUTION, e.score.results)
  
  
}

View(e.NULL.DISTRIBUTION)

p.E <- NULL

for(i in 1:nrow(myData))
{
  
  
  
  
  e.Observed.i <- e.OBSERVED[i,]
  e.NULL.DIST.i <- e.NULL.DISTRIBUTION[i,]
  count2 <- sum((abs(e.Observed.i) < abs(e.NULL.DIST.i)))
  
  
  
  p.E[i] <- (count2 / 100) 
}

saveRDS(p.E, file = "p-empirical-euclidian.rds")
vector.correlation <- (cor(p.T,p.E))

pdf("hist-pT.pdf")
hist(p.T)
dev.off()

pdf("hist-pE.pdf")
hist(p.E)
dev.off()


print(vector.correlation)


#thanks!

