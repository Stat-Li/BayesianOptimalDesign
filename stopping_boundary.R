f_futility = function(gamma_a, gamma_b, nevents, tot_trans, p_lambda, p_gamma, t_j){
  pgamma(1, shape = (gamma_a+nevents), rate = (gamma_b+tot_trans), lower.tail = T) - 
    p_lambda*t_j^p_gamma
}

tot_stop_futility = function(gamma_a, gamma_b, p_lambda, p_gamma, t_interim, event_tot){
  result = rep(NA, length(t_interim))
  for (i in 1:length(t_interim)) {
    nevents = ceiling(event_tot*t_interim[i])
    t_i = t_interim[i]
    res_opt = uniroot(f_futility, interval = c(0, 100), 
                      gamma_a = gamma_a, gamma_b = gamma_b, nevents = nevents, 
                      p_lambda = p_lambda, 
                      p_gamma = p_gamma, t_j = t_i)
    result[i] = res_opt$root
  }
  return(cbind(result))
}

f_superiority = function(gamma_a, gamma_b, nevents, tot_trans, p_lambda, t_j){
  pgamma(1, shape = (gamma_a+nevents), rate = (gamma_b+tot_trans), lower.tail = T) - 
    (2*pnorm(qnorm((1 + p_lambda)/2)/sqrt(t_j)) - 1)
}

tot_stop_superiority = function(gamma_a, gamma_b, p_lambda, t_interim, event_tot){
  result = rep(NA, length(t_interim))
  for (i in 1:length(t_interim)) {
    nevents = ceiling(event_tot*t_interim[i])
    t_i = t_interim[i]
    res_opt = uniroot(f_superiority, interval = c(0, 100), 
                      gamma_a = gamma_a, gamma_b = gamma_b, nevents = nevents, 
                      p_lambda = p_lambda, t_j = t_i)
    result[i] = res_opt$root
  }
  return(cbind(result))
}