#' Computes the vertical wavenumber
#'
#' @export
#' @param ny if ny=1 wavenumber is for a 1D profile
#' @param k sampling freq for x
#' @param l sampling freq for y
#' @param method default is kim's approach. 'qgis' for qgis approach
#' @details the default method is Kim's
#'
#'
#'

vertical_wavenumber <- function(ny,k,l,method=NULL,coriollis=FALSE){

  if(ny!=1){
    sigma = U*k + V*l; # intrinsic frequency (s^-1) - pp. 1379

    if(is.null(method)){

      m = ( ((Nm^2-sigma^2)/sigma^2) * (k^2+l^2) )^0.5 # vertical wavenumber (m^-1)

      #if sigma^2 is > Nm^2 then m = NaN for these cells. Find these and apply a different eqn (in text betwen (12) and (13))
      m_ind <- matrix(seq(1,nxpad*nypad),ncol=nypad)
      m_ind <- as.vector(m_ind*is.na(m))
      m_ind <- m_ind[m_ind>0]

      if(length(m_ind)>0){
        # With weak stratification (Nm^2 << sigma^2), airflow is irrotational and m
        # reduces to
        m[m_ind] = 1i*(k[m_ind]^2+l[m_ind]^2) # vertical wavenumber (m^-1)
      }

    } else if(method=='qgis'){

      sigma_sqr_reg = sigma^2

      if(coriollis){
        m_denom = sigma^2 - 2 * 7.2921e-5 * sin(lat_coriolis * pi / 180)
      } else {
        m_denom = sigma^2
      }

      # sigma_sqr_reg doesn't seem to get used anywhere....
      sigma_sqr_reg[abs(sigma_sqr_reg) < eps & abs(sigma_sqr_reg) >= 0] <- eps
      # sigma_sqr_reg[abs(sigma_sqr_reg) < eps & sigma_sqr_reg < 0] <- -eps #why would this ever happen if sigma was squared?

      # The vertical wave number
      # Eqn. 12
      # Regularization
      m_denom[abs(m_denom) < eps & abs(m_denom) >= 0] = eps
      m_denom[abs(m_denom) < eps & abs(m_denom) < 0] = -eps

      m1 = (Nm^2 - sigma^2)/m_denom
      m2 = k^2 + l^2
      m_sqr = m1*m2
      m = sqrt(-1*m_sqr)

      # Regularization
      m[m_sqr >= 0 & sigma == 0] <- sqrt(m_sqr[m_sqr >= 0 & sigma == 0])
      m[m_sqr >= 0 & sigma != 0] <- sqrt(m_sqr[m_sqr >= 0 & sigma != 0]) * sign(sigma[m_sqr >= 0 & sigma != 0])

    }

  } else {# In 2D hydrostatic case (for a profile instead of 2D DEM) where  # of columns (ny) = 1 and therefore l=0
    sigma = U*k # intrinsic frequency (s^-1)
    m = (Nm/U)*sign(k) # vertical wavenumber (m^-1)
  }

  return(list('m'=m,'sigma'=sigma))
}
