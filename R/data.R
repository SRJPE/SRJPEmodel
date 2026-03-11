
#' @title BT-SPAS-X BUGS parameters
#' @name bt_spas_x_bayes_params
#' @description Parameters for running BUGS on the BT-SPAS-X model
#' @format A named list with 4 elements:
#' \itemize{
#'   \item \code{number_mcmc}: TODO
#'   \item \code{number_burnin}:
#'   \item \code{number_thin}:
#'   \item \code{number_chains}:
#'   }
'bt_spas_x_bayes_params'

#' @title BT-SPAS-X STAN model code
#' @name bt_spas_x_model_code
#' @description A nested named list containing model code for versions of the bt-spas-x models.
#' Model code written in WinBUGS has two versions: "with" and "without" cut functions. The
#' STAN
#' @format A named list with 4 elements:
#' \itemize{
#'   \item \code{pCap_all}: pCap model code
#'   \item \code{abundance}: abundance model code (STAN)
#'   \item \code{abundance_BUGS}: abundance model code (BUGS)
#'   }
'bt_spas_x_model_code'

#' @title Passage to Spawner model code
#' @name p2s_model_code
#' @description A character string containing the STAN model code for the Passage to Spawner submodel.
'p2s_model_code'

#' @title Survival model code
#' @name survival_model_code
#' @description A nested named list containing model code for versions of the survival models.
#' @format A named list with 2 elements:
#' \itemize{
#'   \item \code{survival_CovWY}:
#'   \item \code{survival_NoCov}:
#'   }
'survival_model_code'

