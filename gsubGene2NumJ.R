gsubGene2NumJ <- function(x) {
  x	<-	gsub(	"IGKJ1",	"0.4",	x,	fixed=TRUE)
  x	<-	gsub(	"IGKJ2",	"0.3",	x,	fixed=TRUE)
  x	<-	gsub(	"IGKJ4",	"0.2",	x,	fixed=TRUE)
  x	<-	gsub(	"IGKJ5",	"0.1",	x,	fixed=TRUE)
  x	<-	gsub(	"IGLJ2",	"0.04",	x,	fixed=TRUE)
  x	<-	gsub(	"IGLJ3",	"0.02",	x,	fixed=TRUE)
  x	<-	gsub(	"IGLJ1",	"0.01",	x,	fixed=TRUE)}
