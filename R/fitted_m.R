#' Calculates matrix for reading errors given a set of parameters
#'
#' @param para A vector of length 5 with (beta_0,beta_1,alpha_0,alpha_1,log(phi))
#' @param Amin Mininum age
#' @param Amax Maximum age
#' @return Matrix containing probabilities for error readings
#' @export
fitted_m = function(para)
{
  beta0 = para[1]
  beta1 = para[2]
  alpha0 = para[3]
  alpha1 = para[4]
  phi = exp(para[5])
  fitted_matrix = matrix(0, nrow = 20, ncol = 20)

  data_exp = data_slong[data_slong$reader%in%names(data_s)[7:12] & !is.na(data_slong$age), ]
  true_age = data_exp$`Modal age (exp)`
  age = data_exp$age


  for(i in Amin:Amax){
    paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1*i))
    qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))

    for(j in Amin:Amax){
      if(i==j){
        fitted_matrix[i, j] = paa
      }
      if(i == Amax & j == Amax){
        fitted_matrix[i, j] = qa + paa - qa * paa
      }
      if(i == Amin & j == Amin){
        fitted_matrix[i, j] = 1 - qa + paa * qa
      }
      if(j > i){
        if(j == Amax){
          fitted_matrix[i, j] = (1 - paa) * qa * phi ^ (Amax - i - 1)
        }else{
          fitted_matrix[i, j] = (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
        }

      }
      if(i > j){
        d = i - Amin
        fitted_matrix[i, j] = (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
      }
    }
  }

  fitted_matrix = fitted_matrix[rowSums(fitted_matrix) != 0, rowSums(fitted_matrix) != 0]
  rownames(fitted_matrix) = Amin:Amax
  colnames(fitted_matrix) = Amin:Amax
  fitted_matrix
}
