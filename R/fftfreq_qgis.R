#' Computes k and l matrixes for FFT
#'
#' @param nxpad number of columns
#' @param nypad number of rows
#' @param dx distance per grid cell in x direction
#' @param dy distance per grid cell in y direction


fftfreq_qgis <- function(dx,dy,nxpad,nypad){

  ### This function replicates np.fftfreq from python
  np_fftfreq <- function(n,d=1){
    # n is length of vector desired
    # d is
    if(n!=floor(n)) stop()
    val=1/(n*d)
    results=rep(NA,n)
    N=floor((n-1)/2)+1
    p1=seq(0,N)
    results[1:(N+1)]=p1
    p2=seq(-(floor(n/2)),-1)
    results[(N+1):length(results)]=p2
    results*val
  }

  x_n_value <- np_fftfreq(nypad,1/nxpad)
  y_n_value <- np_fftfreq(nxpad,1/nypad)

  x_len <- nxpad*dx
  y_len <- nypad*dy

  kx_line <- 2*pi*x_n_value/x_len
  ky_line <- 2*pi*y_n_value/y_len


  # # # https://blogs.mathworks.com/steve/2017/06/13/frequency-samples-for-the-output-of-fft2/
  # # # for even N: fa = ((-N/2):(N/2-1))/N but you can generalize for even or odd N

  # Fsx = 1/dx; # sampling frequency (m^-1). smallest change we could detact is intuitively 1/dx which is linked to the resolution of the dem
  # Fsy = 1/dy;
  #
  # fax = seq(-ceiling((nypad-1)/2),floor((nypad-1)/2))/nypad
  # fay = seq(-ceiling((nxpad-1)/2),floor((nxpad-1)/2))/nxpad
  # kx_line = fax*Fsx;
  # ky_line = fay*Fsy;

  k=matrix(rep(kx_line,each=nxpad),ncol=nypad,byrow=F)
  l=matrix(rep(ky_line,each=nypad),ncol=nypad,byrow=T)

  return(list('k'=k,'l'=l))
}
