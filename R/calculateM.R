#' Calculation of transition probabilities for age error readings.
#'
#' Default values on parameters are from fitted model to ...
#' @param beta0 beta0
#' @param beta1 beta1
#' @param alpha0 alpha0
#' @param alpha1 alpha1
#' @param phi phi
#' @return Transition matrix
#' @export
calculateM = function(beta0=4.699026,beta1=-0.3507445,alpha0=0.7443011,alpha1=-0.0381708,phi=0.325192)
{
  fitted_matrix = matrix(0, nrow = 20, ncol = 20)

  data_exp = data_slong[data_slong$reader%in%names(data_s)[7:12] & !is.na(data_slong$age), ]
  true_age = data_exp$`Modal age (exp)`
  age = data_exp$age

  beta0 = fit_m$par[1]
  beta1 = fit_m$par[2]
  alpha0 = fit_m$par[3]
  alpha1 = fit_m$par[4]
  phi = fit_m$par[5]

  for(i in 3:15){
    paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1*i))
    qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))

    for(j in 3:15){
      if(i==j){
        fitted_matrix[i, j] = paa
      }
      if(i == 15 & j == 15){
        fitted_matrix[i, j] = qa + paa - qa * paa
      }
      if(i == 3 & j == 3){
        fitted_matrix[i, j] = 1 - qa + paa * qa
      }
      if(j > i){
        if(j == 15){
          #fitted_matrix[i, j] = (1 - paa) * qa * (1+phi^15-phi^i)
          fitted_matrix[i, j] = (1 - paa) * qa * phi ^ (15 - i - 1)
        }else{
          fitted_matrix[i, j] = (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
        }

      }
      if(i > j){
        d = i - 3
        fitted_matrix[i, j] = (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
      }
    }
  }
  fitted_matrix = fitted_matrix[rowSums(fitted_matrix) != 0, rowSums(fitted_matrix) != 0]
  rownames(fitted_matrix) = 3:15
  colnames(fitted_matrix) = 3:15
  fitted_matrix
}
