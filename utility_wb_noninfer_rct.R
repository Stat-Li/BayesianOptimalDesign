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
prob_post = function(theta, theta_hat, theta_0, sig2_0, nevent, clt_r){
    sig2 = 1/nevent/clt_r/(1-clt_r)
    sig2_tilde = 1/(1/sig2 + 1/sig2_0)
    theta_tilde = sig2_tilde*(theta_hat/sig2 + theta_0/sig2_0)
    prob = pnorm(theta, mean = theta_tilde, sd = sqrt(sig2_tilde), lower.tail = F)
    return(prob)
}
# prob = function(delta, a, b, m, k){
#     pgamma(delta, shape = (a + m), rate = (b + k), lower.tail = TRUE)
# }
# U_total = function(t_obs, S0, t, p_shape){
#     w_i = -log(S0)*(t_obs/t)^p_shape
#     U = sum(w_i)
#     return(U)
# }
pmin_stop = function(matrix_stop){
    matrix_stop = t(t(matrix_stop)*(1:dim(matrix_stop)[2]))
    matrix_stop[matrix_stop == 0] = NA
    res=do.call(pmin, c(as.data.frame(matrix_stop), list(na.rm=TRUE)))
    return(res)
}
# mi = 15; df = df; df_sort = df_sort; df_sort1 = df_sort1
# ta = 4; tf = 2
stage_i = function(mi, df, df_sort, df_sort1, ta, tf){
    if(dim(df_sort1)[1] < mi) {
        df$surv[df$cnsr==0] = 
            (ta + tf) - df$enrollT[df$cnsr==0]
        mi_event = dim(df_sort1)[1]
        size_i = dim(df)[1]; duration_i = ta + tf
        # estimate log hazard ratio using cox model
        fit.cox = coxph(Surv(df$surv, df$cnsr==1) ~ as.factor(df$arm))
        theta_hat = -as.numeric(fit.cox$coefficients)
        clt_r = sum(df$arm=='clt')/size_i
    } else {
        duration_i = df_sort1$studyT_obs[mi]
        enroll_cut = df_sort1$enrollT[mi]
        data_n = df_sort[df_sort$enrollT<=duration_i & df_sort$studyT_obs>duration_i, ]
        data_n$surv = duration_i - data_n$enrollT
        data_n$cnsr = 0; data_n$studyT_obs = duration_i
        df_all = rbind(df_sort1[1:mi, ], data_n)
        size_i = dim(df_all)[1]
        mi_event = mi
        fit.cox = coxph(Surv(df_all$surv, df_all$cnsr==1) ~ as.factor(df_all$arm))
        theta_hat = -as.numeric(fit.cox$coefficients)
        clt_r = sum(df_all$arm=='clt')/size_i
    }
    return(list(n.events = mi_event, theta_hat = theta_hat, clt_ratio = clt_r,
                duration_i = duration_i, size_i = size_i))
}

# delta_0=0.6; p_shape=1; m=20; nclt=20; ntrt=30
# S0 = 0.53; t0 = 3; ti=c(0.5, 1); ta = 4; tf = 2
optim_tune = function(delta_0, p_shape=1, m, nclt, ntrt, ti=c(0.5, 1), 
                      ta, tf, S0, t0, theta_0, sig2_0, nsim = 10000,
                      theta_noninferior){
    p_scale0 = (-log(S0))^(1/p_shape)/t0
    p_scale1 = p_scale0*delta_0^(1/p_shape)
    m_interim = ceiling(m*ti)
    results_post = matrix(NA, nrow = nsim, ncol = length(ti))
    results_dur = matrix(NA, nrow = nsim, ncol = length(ti))
    results_size = matrix(NA, nrow = nsim, ncol = length(ti))
    results_event = matrix(NA, nrow = nsim, ncol = length(ti))
    for(i in 1:nsim){
        set.seed(i)
        surv0_clt = (-log(runif(nclt)))^(1/p_shape)/p_scale0
        surv0_trt = (-log(runif(ntrt)))^(1/p_shape)/p_scale1
        enrollT_clt = runif(nclt, 0, ta)
        enrollT_trt = runif(ntrt, 0, ta)
        df_clt = data.frame(surv = surv0_clt, enrollT = enrollT_clt, arm = "clt")
        df_trt = data.frame(surv = surv0_trt, enrollT = enrollT_trt, arm = "trt")
        df = rbind(df_clt, df_trt)
        df$studyT = df$surv + df$enrollT
        df$studyT_obs = pmin(df$studyT, (ta + tf))
        df$cnsr = 1*(df$studyT < (ta + tf))
        df_sort = df[order(df$studyT), ]
        df_sort1 = df_sort[which(df_sort$cnsr==1), ]
        
        for (j in 1:length(ti)) {
            mi = m_interim[j]
            result_j = stage_i(mi=mi, df = df, df_sort = df_sort, df_sort1 = df_sort1,
                               ta = ta, tf = tf)
            results_post[i,j] = prob_post(theta = theta_noninferior,
                                          theta_hat=result_j$theta_hat, 
                                          theta_0 = theta_0, sig2_0 = sig2_0, 
                                          nevent = result_j$n.events, 
                                          clt_r = result_j$clt_ratio)
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
    probs = results$prob
    
    # -- futility stop -----#
    cutoffs_fut = cut_futility(ti, p_lambda = p_lambda, p_gamma = p_gamma)
    prob_check_fut = t(t(probs) < cutoffs_fut)

    p_type1 = mean(rowSums(prob_check_fut)==0)
    return(p_type1)
}

grid_func_size_fut = function(param_tune, results, ti = c(0.5, 1)){
    p_lambda = param_tune[1]; p_gamma = param_tune[2]
    probs = results$prob; size = results$size
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
