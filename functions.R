
annotationMatrixFromLD = function( ldframe ) {
  require(magrittr)
  aid = unique( ldframe[,1] ) %>% setNames(1:length(.), . )
  gid = unique( ldframe[,2] ) %>% setNames(1:length(.), . )

  amat = Matrix( nrow = length(aid), 
                 ncol= length(gid), 
                 0, 
                 dimnames = list(names(aid), names(gid)) )
  
  amat[ cbind( aid[ as.character(ldframe[,1]) ], gid[ as.character(ldframe[,2]) ]) ] = 1
  
  return(amat)
}




rs.signifDataframe = function(x, digits=2) {
  for (i in 1:ncol(x)) {
    if (class(x[,i]) == "numeric") {
      x[,i] = signif(x[,i], digits)
    }
  }
  return(x)
}