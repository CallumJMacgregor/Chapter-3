# creates a new data frame from an existing one containing repeated
# measures as columns, e.g. a data frame named exp1.df:
# subject  age	repeat1	repeat2
# 001		34	1.45	1.67
# 002		38	1.20	1.54
# ...
# called like:
# make.rm(1:2,3:4,exp1.df)
# would multiply the "constant" variables by the number of
# repeats and reformat the repeats to a single column.
# subject	age	repdat	contrasts
# 001		34	1.45	T1
# 002		38	1.20	T1
# ...
# 001		34	1.67	T2
# 002		38	1.54	T2
# ...
# this allows a "univariate" repeated measures analysis of the data.

make.rm<-function(constant,repeated,data,contrasts) {
  if(!missing(constant) && is.vector(constant)) {
    if(!missing(repeated) && is.vector(repeated)) {
      if(!missing(data)) {
        dd<-dim(data)
        replen<-length(repeated)
        if(missing(contrasts))
          contrasts<-
          ordered(sapply(paste("T",1:length(repeated),sep=""),rep,dd[1]))
        else
          contrasts<-matrix(sapply(contrasts,rep,dd[1]),ncol=dim(contrasts)[2])
        if(length(constant) == 1) cons.col<-rep(data[,constant],replen)
        else cons.col<-lapply(data[,constant],rep,replen)
        new.df<-data.frame(cons.col,
                           repdat=as.vector(data.matrix(data[,repeated])),
                           contrasts)
        return(new.df)
      }
    }
  }
  cat("Usage: make.rm(constant, repeated, data [, contrasts])\n")
  cat("\tWhere 'constant' is a vector of indices of non-repeated data and\n")
  cat("\t'repeated' is a vector of indices of the repeated measures data.\n")
}