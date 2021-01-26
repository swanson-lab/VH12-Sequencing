gsubNum2GeneJ <- function(x) {
  x	<-	gsub(	"0.4",	"IGKJ1",	x,	fixed=TRUE)
  x	<-	gsub(	"0.3",	"IGKJ2",	x,	fixed=TRUE)
  x	<-	gsub(	"0.2",	"IGKJ4",	x,	fixed=TRUE)
  x	<-	gsub(	"0.1",	"IGKJ5",	x,	fixed=TRUE)
  x	<-	gsub(	"0.04",	"IGLJ2",	x,	fixed=TRUE)
  x	<-	gsub(	"0.02",	"IGLJ3",	x,	fixed=TRUE)
  x	<-	gsub(	"0.01",	"IGLJ1",	x,	fixed=TRUE)
}