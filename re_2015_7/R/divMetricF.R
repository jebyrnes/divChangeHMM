# RE 150619
# diversity metric function

# x = dataframe
# first1 = the first column of the species matrix within the data frame
# last1 = the last column of the species matrix within the data frame
# returns the dataframe with richness, Shannon, and Pielou's evenness appended 
# at the end of the dataframe.

divMetF <- function(x, first1, last1) {
  specMat <- x[, first1:last1]
  x$rich <- specnumber(specMat)
  x$div <- diversity(specMat)
  x$even <- x$div/(log(x$rich))
  x$abund <- rowSums(specMat)
  return(x)
}