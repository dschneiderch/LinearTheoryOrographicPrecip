#' Computes k and l for FFT
#'
#' @param dx resolution in x direction
#' @param dy resolution in y direction
#' @param nxpad number of columns in dem
#' @param nypad number of rows in in dem
#' @details as the name suggests, this is Kim's approach. Don't use this because it doesn't seem to replicate the guassian hill example.

fftfreq_kim <- function(dx,dy,nxpad,nypad){
  ## ***** Kim's approach
  Fsx=1/dx
  Fsy=1/dy
  # my notation is different than Kims and Smith&Barstad2004.
  # Kim: nx is number of columns (or # of x coordinates), ny is number of rows (or # of y coordinates)
  # me: ny is the number of columns (or "# of lines of y dimension"), nx is number of rows (or "# of lines in x dimension")
  kxv = Fsx*seq(0,1,length.out=nxpad)
  lyv = Fsy*seq(0,1,length.out=nypad)

  k=matrix(rep(kxv,each=nxpad),ncol=nypad,byrow=F) #since kxv has length ny, repeat each nx times.
  l=matrix(rep(lyv,each=nypad),ncol=nypad,byrow=T) #since lyv has length nx, repeat each ny times.
  # The dimensions of k and l should be the same as original h. Since there are only ny unique values in kxv,  each row must be the same. since there are only nx unique values in lyv, each column must be the same.
  # With byrow=T for lyv, output matches matlab::meshgrid from Kim's code.
  # matlab::meshgrid(kxv,lyv) # matlab meshgrid function creates a grid for both dimensions such that all combinations are represented.

  return(list('k'=k,'l'=l))
}
