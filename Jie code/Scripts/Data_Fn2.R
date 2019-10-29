
#' Build data input for SSST
#'
#' \code{Data_Fn} builds a tagged list of data inputs used by TMB for running the model
#'
#' @param Version a version number 
#' @param obsmodel_p The observation model for each size class p
#' \describe{
#'   \item{ObsModel=0}{Poisson}
#'   \item{ObsModel=1}{Lognormal}
#'   \item{ObsModel=2}{Zero-inflated lognormal}
#'   \item{ObsModel=3}{lognormal-Poisson}
#'   \item{ObsModel=4}{Normal}
#' }
#' @param b_i Sampled N per unit area for each observation i
#' @param s_i Spatial knot (e.g., grid cell) for each observation i
#' @param t_i Time interval (e.g., year) for each observation i
#' @param p_i size class for each observation i
#' @param a_x Area associated with each knot
#' @param n_factors Rank of covariance matrix for process error
#' 
#' @param c_pkt harvest (N per area, in same units as \code{b_i} and \code{a_x}), where \code{c_pkt=0} involves no harvest (the default)
#' @param catchlogsd log-standard deviation of observed harvest specified via \code{c_pkt}
#' @param MeshList, tagged list representing location information for the SPDE mesh hyperdistribution, i.e., from \code{SpatialDeltaGLMM::Spatial_Information_Fn}
#' @param spatial_method DEPRECATED, always uses "Mesh" approximation
#' @param CheckForErrors Boolean, whether to check for errors in data inputs

#' @return Tagged list containing inputs to function \code{MIST::Build_TMB_Fn()}

#' @export
Data_Fn = function(Version, obsmodel_p=NULL, b_i, s_i, t_i, p_i, a_x, n_r, R_sex, R_size, mature, size_midpoints, M_pkt=NULL, growthmale_tran,  MeshList, n_factors=1, c_pkt=NULL, catchlogsd, CheckForErrors=TRUE ){

  # Assemble options vector
  options_vec = c(0,0)

  # Expand a_x for auxiliary knots
  a_k = c(a_x, rep(0,MeshList$anisotropic_mesh$n-nrow(MeshList$loc_x)))

  # Check for bad data entry
  data_frame = cbind( "sitenum"=s_i, "year"=t_i, "size_class"=p_i, "catch"=b_i )
  if( CheckForErrors==TRUE ){
    #if( !all(length(b_i)!=n_i | length(a_i)!=n_i | length(v_i)!=n_i | length(s_i)!=n_i | length(t_i)!=n_i ) stop("b_i, a_i, v_i, s_i, or t_i doesn't have length n_i")
  }

  # Fill in defaults
  if( is.null(c_pkt)) c_pkt = array(0, dim=c(length(unique(data_frame[,'size_class'])),MeshList$anisotropic_mesh$n,diff(range(data_frame[,'year']))+1) )
  if( is.null(M_pkt)) M_pkt = array(0.25, dim=c(length(unique(data_frame[,'size_class'])),MeshList$anisotropic_mesh$n,diff(range(data_frame[,'year']))+1) )
 
   # Spatial
  if( Version%in%c("spatial_size_8","spatial_size_9")){
    require( Matrix )
    # Sparse matrices for 2D AR1 process
    # Generate sparse matrices for precision matrix of 2D AR1 process
    M0 = Spatial_List$GridList[["M0"]]
    M1 = Spatial_List$GridList[["M1"]]
    M2 = Spatial_List$GridList[["M2"]]

    # Data
    # Necessary columns: sitenum, year, catch, spp
    if(Version%in%c("spatial_size_8","spatial_size_9")) Return = list("Options_vec"=options_vec, "ObsModel_p"=obsmodel_p, "n_i"=nrow(data_frame), "n_s"=length(unique(data_frame[,'sitenum'])), "n_t"=diff(range(data_frame[,'year']))+1, "n_k"=MeshList$anisotropic_mesh$n, "n_p"=length(unique(data_frame[,'size_class'])), "n_j"=n_factors, "n_r"=n_r, "b_i"=data_frame[,'catch'], "p_i"=as.numeric(data_frame[,'size_class'])-1, "s_i"=data_frame[,'sitenum']-1, "t_i"=data_frame[,'year']-min(data_frame[,'year']), "R_size"=R_size, "mature"=mature, "size_midpoints"=size_midpoints, "M_pkt"=M_pkt, "growthmale_tran"=growthmale_tran, "c_pkt"=c_pkt, "catchlogsd"=catchlogsd, "a_k"=a_k, "spde"=NULL, "M0"=M0, "M1"=M1, "M2"=M2)
    
    if("spde" %in% names(Return)) Return[['spde']] = list("n_s"=MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)) )

  } # End spatial

  return(Return)
}

