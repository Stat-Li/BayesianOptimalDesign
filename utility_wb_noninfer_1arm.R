# Gamma family
cut_futilit1 = function(t_i, p_lambda, p_gamma){
    p_l = p_lambda*(1 - exp(-p_gamma*t_i))/(1 - exp(-p_gamma))
    return(p_l)
}
# power function
cut_futility = function(t_i, p_lambda, p_gamma){
    p_l = p_lambda*t_i^p_gamma
    return(p_l)
}
cut_super = function(t_i, p_lambda){
    p_l = 2*pnorm(qnorm((1 + p_lambda)/2)/sqrt(t_i)) - 1
    return(p_l)
}
prob = function(delta, a, b, m, k){
    pgamma(delta, shape = (a + m), rate = (b + k), lower.tail = TRUE)
}
U_total = function(t_obs, S0, t, p_shape){
    w_i = -log(S0)*(t_obs/t)^p_shape
    U = sum(w_i)
    return(U)
}
pmin_stop = function(matrix_stop){
    matrix_stop = t(t(matrix_stop)*(1:dim(matrix_stop)[2]))
    matrix_stop[matrix_stop == 0] = NA
    res=do.call(pmin, c(as.data.frame(matrix_stop), list(na.rm=TRUE)))
    return(res)
}
stage_i = function(mi, S0, t, dataSurv_sort1_s1, dataSurv_s1, dataSurv_sort_s1,
                   p_shape, ta, tf){
    if(dim(dataSurv_sort1_s1)[1] < mi) {
        dataSurv_s1$surv_obs[dataSurv_s1$cnsr==0] = 
            (ta + tf) - dataSurv_s1$enrollT[dataSurv_s1$cnsr==0]
        mi_event = dim(dataSurv_sort1_s1)[1]
        size_i = dim(dataSurv_s1)[1]; duration_i = ta + tf
        data_obs_surv = dataSurv_s1$surv_obs
    } else {
        duration_i = dataSurv_sort1_s1$studytime[mi]
        enroll_cut = dataSurv_sort_s1$enrollT[which(dataSurv_sort_s1$studytime==
                                                        dataSurv_sort1_s1$studytime[mi])]
        data_n = dataSurv_sort_s1[dataSurv_sort_s1$enrollT<=dataSurv_sort1_s1$studytime[mi] & 
                                      dataSurv_sort_s1$studytime>=dataSurv_sort1_s1$studytime[mi], ]
        data_n$surv_obs = dataSurv_sort1_s1$studytime[mi] - data_n$enrollT
        data_obs_surv = c(dataSurv_sort1_s1$surv_obs[1:(mi-1)], data_n$surv_obs)
        size_i = dim(data_n)[1] + mi - 1
        mi_event = mi
    }
    tot_obs = U_total(data_obs_surv, S0=S0, t=t, p_shape=p_shape)
    return(list(n.events = mi_event, Uobs = tot_obs, 
                duration_i = duration_i, size_i = size_i))
}

optim_tune = function(delta_0, p_shape=1 , m, n, ti=c(0.5, 1), pram_gamma_a, 
                      pram_gamma_b, ta = 4, tf = 2, S0 = 0.53, t0 = 3, 
                      nsim = 10000, delta_noninferior){
    p_scale0 = (-log(0.53))^(1/p_shape)/3
    p_scale1 = p_scale0*delta_0^(1/p_shape)
    m_interim = ceiling(m*ti)
    results_post = matrix(NA, nrow = nsim, ncol = length(ti))
    results_dur = matrix(NA, nrow = nsim, ncol = length(ti))
    results_size = matrix(NA, nrow = nsim, ncol = length(ti))
    results_event = matrix(NA, nrow = nsim, ncol = length(ti))
    for(i in 1:nsim){
        set.seed(i)
        surv0_s1 = (-log(runif(n)))^(1/p_shape)/p_scale1
        enrollT_s1 = runif(n, 0, ta)
        studyT_s1 = surv0_s1 + enrollT_s1
        studyT_obs_s1 = pmin(studyT_s1, (ta + tf))
        ind_cnsr_s1 = 1*(studyT_s1 < (ta + tf))
        dataSurv_s1 = data.frame(studytime = studyT_obs_s1, cnsr = ind_cnsr_s1,
                                 surv_obs = surv0_s1, enrollT = enrollT_s1)
        dataSurv_sort_s1 = dataSurv_s1[order(dataSurv_s1$studytime),]
        dataSurv_sort1_s1 = dataSurv_sort_s1[which(dataSurv_sort_s1$cnsr==1),]
        for (j in 1:length(ti)) {
            mi = m_interim[j]
            result_j = stage_i(mi, S0, t0, dataSurv_sort1_s1 = dataSurv_sort1_s1, 
                               dataSurv_s1 = dataSurv_s1, 
                               dataSurv_sort_s1 = dataSurv_sort_s1, p_shape = p_shape,
                               ta = ta, tf = tf)
            results_post[i,j] = prob(delta = delta_noninferior, a = pram_gamma_a, 
                                     b = pram_gamma_b, m = result_j$n.events, 
                                     k = result_j$Uobs)
            results_dur[i,j] = result_j$duration_i
            results_size[i,j] = result_j$size_i
            results_event[i,j] = result_j$n.events
        }
    }
    return(list(prob = results_post, dur = results_dur, 
                size = results_size, event = results_event))
}

grid_func_size = function(param_tune, results, ti = c(0.5, 1)){
    p_lambda = param_tune[1]; p_gamma = param_tune[2]
    probs = results$prob; size = results$size; 
    dur = results$dur; event = results$event
    # -- superiority stop -----#
    cutoffs_super = cut_super(ti, p_lambda = p_lambda)
    prob_check_super = t(t(probs) >= cutoffs_super)
    # PRN = 1 - mean(rowSums(prob_check_super) == 0) # power
    prob_check1 = prob_check_super
    prob_check1[, dim(prob_check1)[2]] = 1
    stop_i_super = pmin_stop(prob_check1) # superiority stop at interim i
    # -- futility stop -----#
    cutoffs_fut = cut_futility(ti, p_lambda = p_lambda, p_gamma = p_gamma)
    prob_check_fut = t(t(probs) < cutoffs_fut)
    prob_check2 = prob_check_fut
    prob_check2[, dim(prob_check2)[2]] = 1
    stop_i_fut = pmin_stop(prob_check2) # futility stop at interim i
    
    PRN = mean(stop_i_super<stop_i_fut) + 
        mean((stop_i_super==stop_i_fut)*
                 (rowSums(prob_check_fut)==0))
    PET_EFS = mean(stop_i_fut<stop_i_super)
    PET_ESS = mean(stop_i_fut>stop_i_super)
    PET = PET_EFS + PET_ESS
    stop_i = pmin(stop_i_fut, stop_i_super)
    size_i = apply(cbind(1:dim(prob_check1)[1],stop_i),1, function(x) size[x[1],x[2]])
    dur_i = apply(cbind(1:dim(prob_check1)[1],stop_i),1, function(x) dur[x[1],x[2]])
    event_i = apply(cbind(1:dim(prob_check1)[1],stop_i),1, function(x) event[x[1],x[2]])
    size_avg = mean(size_i)
    dur_avg = mean(dur_i)
    event_avg = mean(event_i)
    return(c(PRN, PET, PET_EFS, PET_ESS, size_avg, dur_avg, event_avg))
}
grid_func = function(param_tune, results, ti = c(0.5, 1)){
    # x=p_lambda; y=p_gamma
    p_lambda = param_tune[1]; p_gamma = param_tune[2]
    probs = results$prob; size = results$size; dur = results$dur
    # -- superiority stop -----#
    cutoffs_super = cut_super(ti, p_lambda = p_lambda)
    prob_check_super = t(t(probs) >= cutoffs_super)
    prob_check1 = prob_check_super
    prob_check1[, dim(prob_check1)[2]] = 1
    # stop_i_super=apply(prob_check1,1, function(x) which(x!=0)[1]) # superiority stop at interim i
    stop_i_super = pmin_stop(prob_check1)
    # -- futility stop -----#
    cutoffs_fut = cut_futility(ti, p_lambda = p_lambda, p_gamma = p_gamma)
    prob_check_fut = t(t(probs) < cutoffs_fut)
    prob_check2 = prob_check_fut
    prob_check2[, dim(prob_check2)[2]] = 1
    stop_i_fut = pmin_stop(prob_check2)
    # stop_i_fut=apply(prob_check2,1, function(x) which(x!=0)[1]) # futility stop at interim i
    
    p_type1 = mean(stop_i_super<stop_i_fut) + 
        mean((stop_i_super==stop_i_fut)*
                 (rowSums(prob_check_fut)==0))
    return(p_type1)
}
grid_func_fut = function(param_tune, results, ti = c(0.5, 1)){
    # x=p_lambda; y=p_gamma
    p_lambda = param_tune[1]; p_gamma = param_tune[2]
    probs = results$prob; size = results$size; dur = results$dur
    
    # -- futility stop -----#
    cutoffs_fut = cut_futility(ti, p_lambda = p_lambda, p_gamma = p_gamma)
    prob_check_fut = t(t(probs) < cutoffs_fut)
    # prob_check2 = prob_check_fut
    # prob_check2[, dim(prob_check2)[2]] = 1
    # stop_i_fut = pmin_stop(prob_check2)

    p_type1 = mean(rowSums(prob_check_fut)==0)
    return(p_type1)
}

grid_func_size_fut = function(param_tune, results, ti = c(0.5, 1)){
    p_lambda = param_tune[1]; p_gamma = param_tune[2]
    probs = results$prob; size = results$size; 
    dur = results$dur; event = results$event
    
    # -- futility stop -----#
    cutoffs_fut = cut_futility(ti, p_lambda = p_lambda, p_gamma = p_gamma)
    prob_check_fut = t(t(probs) < cutoffs_fut)
    prob_check2 = prob_check_fut
    prob_check2[, dim(prob_check2)[2]] = 1
    stop_i_fut = pmin_stop(prob_check2) # futility stop at interim i
    
    PRN = mean(rowSums(prob_check_fut)==0) 
    PET = mean(stop_i_fut < length(ti))
    stop_i = stop_i_fut
    size_i = apply(cbind(1:dim(prob_check2)[1],stop_i),1, function(x) size[x[1],x[2]])
    dur_i = apply(cbind(1:dim(prob_check2)[1],stop_i),1, function(x) dur[x[1],x[2]])
    event_i = apply(cbind(1:dim(prob_check2)[1],stop_i),1, function(x) event[x[1],x[2]])
    size_avg = mean(size_i)
    dur_avg = mean(dur_i)
    event_avg = mean(event_i)
    return(c(PRN, PET, size_avg, dur_avg, event_avg))
}
