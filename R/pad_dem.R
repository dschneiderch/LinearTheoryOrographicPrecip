#' Pad DEM with 0 to remove edge effects from FFT ----
#'
#' @export
#' @param h matrix of topo
#' @param pad (optional) 0 for no padding. Leave blank for automatic padding.
#' @details This will also remove any NA and convert to 0
#'

pad_dem <- function(h,pad=NULL){

  pad_max=250

  nx=dim(h)[1]
  ny=dim(h)[2]

  if(is.null(pad)){
    pad=ceiling(((nx+ny)/2)/100)*100
    if ( pad > pad_max){
      pad=pad_max
    }
  }

  hpad=matrix(0,nrow=nx+2*pad,ncol=ny+2*pad)
  h[is.na(h)] <- 0
  hpad[(pad+1):(nx+pad),(pad+1):(ny+pad)] <- h


  return(list('dem'=hpad,'pad'=pad))

  }
