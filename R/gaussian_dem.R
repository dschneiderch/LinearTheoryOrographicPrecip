#' Replicates the topography used in Figures 3 and 4 in Smith and Barstad 2004
#'
#' @export
#' @param dx resolution along x axis
#' @param dy resolution along y axis

gaussian_dem <- function(dx, dy){
  xv = seq(-100e3,200e3,dx)
  yv = seq(-150e3,150e3,dy)
  nx = length(yv)
  ny = length(xv)
  h_max = 500
  x0 = -25e3
  y0 = 0
  sigma_x = 15e3
  sigma_y = 15e3
  X=matrix(rep(xv,each=ny),ncol=ny,byrow=F) #since kxv has length ny, repeat each nx times.
  Y=matrix(rep(yv,each=nx),ncol=ny,byrow=T) #since lyv has length nx, repeat each ny times.
  # ll=matlab::meshgrid(xv, yv) #use the matlab command meshgrid if you want.
  h = h_max * exp(-(((X - x0)**2 / (2 * sigma_x**2)) + ((Y - y0)**2 / (2 * sigma_y**2))))
#
  # image(h)

  return(list('xcoords'=xv,'ycoords'=yv,'zcoords'=h))
  }


