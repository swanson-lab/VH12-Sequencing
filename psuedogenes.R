psuedogenes<- function(x) {
  x <- gsub( "IGKV4-56", "\U03A8 IGKV4-56", x, fixed=TRUE)
  x <- gsub( "IGKV4-60", "\U03A8 IGKV4-60", x, fixed=TRUE)
  x <- gsub( "IGKV4-75", "\U03A8 IGKV4-75", x, fixed=TRUE)
  x <- gsub( "IGKV4-77", "\U03A8 IGKV4-77", x, fixed=TRUE)
}