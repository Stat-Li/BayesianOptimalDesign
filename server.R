library(survival)
server <- function(input, output, session) {
  

  df_out <- eventReactive(input$go, {
    
    design_stop <- input$design_stop
    delta0 <- input$delta0
    delta_1 <- input$delta1
    type1 <- input$type1
    type2 <- input$type2
    weibull_shape <- input$weibull_shape
    ta <- input$ta
    tf <- input$tf
    surv_3year <- input$surv_3year
    x_year <- input$x_year
    n_sims <- input$n_sims
    timing = c(as.numeric(unlist(strsplit(input$interim_timing," "))), 1)
    results = NULL
    results_tune = NULL
    results_stop = NULL
    if(input$design_arm=="Single-arm"){
      source("utility_wb_noninfer_1arm.R")
      pram_gamma_a = as.numeric(unlist(strsplit(input$gamma_a," ")))
      m_target = ceiling((qnorm(1-type1)+qnorm(1-type2))^2/(log(delta_1/delta0))^2)
      ff = function(t){
        exp(log(surv_3year)*(t/x_year)^weibull_shape)^delta_1
      }
      n_max=ceiling(m_target/(1-1/ta*integrate(ff,tf,ta+tf)$value))
      for (h in 1:length(pram_gamma_a)) {
        
        pram_gamma_b = 2*pram_gamma_a[h]/(1 + delta_1)
        
        results_h0 = optim_tune(delta_0 = delta0, m = m_target, n=n_max, ti=timing, nsim = n_sims,
                                pram_gamma_a = pram_gamma_a[h], pram_gamma_b = pram_gamma_b,
                                p_shape = weibull_shape, ta = ta, tf = tf, delta_noninferior = delta0,
                                S0 = surv_3year, t0 = x_year)
        
        
        p1 = seq (0.01, 0.99, by = 0.01)
        p2 = seq (0.01, 3, by = 0.01)
        
        p.tune = expand.grid(x = p1, y = p2)
        if(input$design_stop == "Futility and superiority stopping"){
          grid_search = apply(p.tune, 1, grid_func, results = results_h0, ti = timing)
          grid_search = round(grid_search, digits = 3)
          sum_0=p.tune[which(abs(grid_search-type1) < 0.01 & grid_search < type1), ]
          # dim(sum_0)
          
          results_h1 = optim_tune(delta_0 = delta_1, m = m_target, n=n_max, ti=timing, nsim = n_sims,
                                  pram_gamma_a = pram_gamma_a[h], pram_gamma_b = pram_gamma_b,
                                  p_shape = weibull_shape, ta = ta, tf = tf, delta_noninferior = delta0,
                                  S0 = surv_3year, t0 = x_year)
          
          grid_search_power=t(apply(sum_0, 1, grid_func_size, results = results_h1,
                                    ti = timing))
          grid_search_power = round(grid_search_power, digits = 3)
          sum_1=sum_0[which(grid_search_power[,1] >= (1-type2)), ]
          # dim(sum_1)
          if(dim(sum_1)[1]==0){
            sum_1 = sum_0[which(grid_search_power[,1] == max(grid_search_power[,1])),]
          }
          #--------------- minimize expected sample size under H0 --------------------
          
          grid_search_all = t(apply(sum_1, 1, grid_func_size, results = results_h0, ti = timing))
          grid_search_all = round(grid_search_all, digits = 3)
          p_opt = sum_1[which.min(grid_search_all[,5]),]
          p_opt = as.numeric(p_opt)
          p_cut_fut = cut_futility(timing, p_lambda = p_opt[1], p_gamma = p_opt[2])
          p_cut_sup = cut_super(timing, p_lambda = p_opt[1])
          results_stop = rbind(results_stop, cbind(p_cut_fut, p_cut_sup))
          
          PRN0 = grid_search_all[which.min(grid_search_all[,5]),1]
          PET0 = grid_search_all[which.min(grid_search_all[,5]),2]
          PET_EFS0 = grid_search_all[which.min(grid_search_all[,5]),3]
          PET_ESS0 = grid_search_all[which.min(grid_search_all[,5]),4]
          size0 = grid_search_all[which.min(grid_search_all[,5]),5]
          dur0 = grid_search_all[which.min(grid_search_all[,5]),6]
          event0 = grid_search_all[which.min(grid_search_all[,5]),7]
          grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          if(sum(grid_search_power[,1] >= (1-type2))<=1) {
            grid_search_power_min_size = matrix(grid_search_power[which.max(grid_search_power[,1]),], nrow = 1)
          }else{
            grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          }
          
          PRN1 = grid_search_power_min_size[which.min(grid_search_all[,5]),1]
          PET1 = grid_search_power_min_size[which.min(grid_search_all[,5]),2]
          PET_EFS1 = grid_search_power_min_size[which.min(grid_search_all[,5]),3]
          PET_ESS1 = grid_search_power_min_size[which.min(grid_search_all[,5]),4]
          size1 = grid_search_power_min_size[which.min(grid_search_all[,5]),5]
          dur1 = grid_search_power_min_size[which.min(grid_search_all[,5]),6]
          event1 = grid_search_power_min_size[which.min(grid_search_all[,5]),7]
          xx1=c(PRN0*100, PET0*100, PET_EFS0*100, PET_ESS0*100, size0, event0, dur0)
          xx2=c(PRN1*100, PET1*100, PET_EFS1*100, PET_ESS1*100, size1, event1, dur1)
          
          results = rbind(results, xx1, xx2)
        }
        if(input$design_stop == "Futility stopping only"){
          grid_search = apply(p.tune, 1, grid_func_fut, results = results_h0, ti = timing)
          grid_search = round(grid_search, digits = 3)
          sum_0=p.tune[which(abs(grid_search-type1) < 0.01 & grid_search < type1), ]
          
          results_h1 = optim_tune(delta_0 = delta_1, m = m_target, n=n_max, ti=timing, nsim = n_sims,
                                  pram_gamma_a = pram_gamma_a[h], pram_gamma_b = pram_gamma_b,
                                  p_shape = weibull_shape, ta = ta, tf = tf, delta_noninferior = delta0,
                                  S0 = surv_3year, t0 = x_year)
          
          grid_search_power=t(apply(sum_0, 1, grid_func_size_fut, results = results_h1, 
                                    ti = timing))
          grid_search_power = round(grid_search_power, digits = 3)
          sum_1=sum_0[which(grid_search_power[,1] >= (1-type2)), ]
          
          if(dim(sum_1)[1]==0){
            sum_1 = sum_0[which(grid_search_power[,1] == max(grid_search_power[,1])),]
          }
          #--------------- minimize expected sample size under H0 --------------------
          
          grid_search_all = t(apply(sum_1, 1, grid_func_size_fut, results = results_h0, ti = timing))
          grid_search_all = round(grid_search_all, digits = 3)
          p_opt = sum_1[which.min(grid_search_all[,3]),]
          p_opt = as.numeric(p_opt)
          p_cut_fut = cut_futility(timing, p_lambda = p_opt[1], p_gamma = p_opt[2])
          
          results_stop = rbind(results_stop, cbind(p_cut_fut))
          
          PRN0 = grid_search_all[which.min(grid_search_all[,3]),1]
          PET0 = grid_search_all[which.min(grid_search_all[,3]),2]
          size0 = grid_search_all[which.min(grid_search_all[,3]),3]
          dur0 = grid_search_all[which.min(grid_search_all[,3]),4]
          event0 = grid_search_all[which.min(grid_search_all[,3]),5]
          grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          if(sum(grid_search_power[,1] >= (1-type2))<=1) {
            grid_search_power_min_size = matrix(grid_search_power[which.max(grid_search_power[,1]),], nrow = 1)
          }else{
            grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          }
          
          PRN1 = grid_search_power_min_size[which.min(grid_search_all[,3]),1]
          PET1 = grid_search_power_min_size[which.min(grid_search_all[,3]),2]
          size1 = grid_search_power_min_size[which.min(grid_search_all[,3]),3]
          dur1 = grid_search_power_min_size[which.min(grid_search_all[,3]),4]
          event1 = grid_search_power_min_size[which.min(grid_search_all[,3]),5]
          xx1=c(PRN0*100, PET0*100, size0, event0, dur0)
          xx2=c(PRN1*100, PET1*100, size1, event1, dur1)
          results = rbind(results, xx1, xx2)
        }
        
      }
      if(input$design_stop == "Futility stopping only"){
        df_final_ops = data.frame(cbind(rep(1:length(pram_gamma_a), each = 2), weibull_shape,
                                    rep(pram_gamma_a, each = 2), rep(c(1, delta_1), length(pram_gamma_a)),
                                    results))
        names(df_final_ops) = c('Design', 'Weibull Shape', 'Gamma prior shape',
                            'Delta', 'PRN (%)', 'PET (%)', 
                            'ES', 'EE', 'Duration')
        results_stop = round(results_stop, digits = 2)
        df_final_stop = data.frame(cbind(rep(1:length(pram_gamma_a), each = length(timing)),
                                         rep(pram_gamma_a, each = length(timing)), 
                                         rep(timing, length(pram_gamma_a)),
                                         results_stop))
        names(df_final_stop) = c('Design', 'Gamma prior shape',
                                 'Information time', 'Futility stopping if posterior prob. <')
      }
      if(input$design_stop == "Futility and superiority stopping"){
        df_final_ops = data.frame(cbind(rep(1:length(pram_gamma_a), each = 2), weibull_shape,
                                    rep(pram_gamma_a, each = 2), rep(c(1, delta_1), length(pram_gamma_a)),
                                    results))
        names(df_final_ops) = c('Design', 'Weibull Shape', 'Gamma prior shape',
                            'Delta', 'PRN (%)', 'PET (%)', 'PET futility (%)',
                            'PET superiority (%)', 'ES', 'EE', 'Duration')
        results_stop = round(results_stop, digits = 2)
        df_final_stop = data.frame(cbind(rep(1:length(pram_gamma_a), each = length(timing)),
                                         rep(pram_gamma_a, each = length(timing)), 
                                         rep(timing, length(pram_gamma_a)),
                                         results_stop))
        names(df_final_stop) = c('Design', 'Gamma prior shape',
                                 'Information time', 'Futility stopping if posterior prob. <', 
                                 'Superiority stopping if posterior prob. >=')
      }
    }
    if(input$design_arm=="Two-arm RCT"){
      source("utility_wb_noninfer_rct.R")
      eta = -log(delta0)
      theta0 = -0.5*log(delta_1)
      sigma2_0 = as.numeric(unlist(strsplit(input$sigma2_0," ")))
      clt_ratio = input$ss_ratio/(input$ss_ratio + 1)
      m_target = ceiling((qnorm(1-type1)+qnorm(1-type2))^2/(log(delta_1)+eta)^2/clt_ratio/(1 - clt_ratio))
      S1 = function(t){
        exp(log(surv_3year)*(t/x_year)^weibull_shape*delta_1)
      }
      n_max=ceiling(m_target/(1-1/ta*integrate(S1,tf,ta+tf)$value))
      n_trt = ceiling(n_max*(1 - clt_ratio)); n_clt = ceiling(n_max*clt_ratio)
      for (h in 1:length(sigma2_0)){
        results_h0 = optim_tune(delta_0 = exp(-eta), p_shape=weibull_shape, m = m_target,
                                nclt = n_clt, ntrt = n_trt, ti = timing, S0 = surv_3year, t0 = x_year,
                                ta = ta, tf = tf, theta_0 = theta0, sig2_0 = sigma2_0[h], nsim = n_sims,
                                theta_noninferior = eta)
        p1 = seq (0.01, 0.99, by = 0.01)
        p2 = seq (0.01, 3, by = 0.01)
        p.tune = expand.grid(x = p1, y = p2)
        if(input$design_stop == "Futility and superiority stopping"){
          grid_search = apply(p.tune, 1, grid_func, results = results_h0, ti = timing)
          grid_search = round(grid_search, digits = 3)
          sum_0=p.tune[which(abs(grid_search-type1) < 0.01 & grid_search <= type1), ]
          results_h1 = optim_tune(delta_0 = delta_1, p_shape=weibull_shape, m = m_target,
                                  nclt = n_clt, ntrt = n_trt, ti = timing, S0 = surv_3year, t0 = x_year,
                                  ta = ta, tf = tf, theta_0 = theta0, sig2_0 = sigma2_0[h], nsim = n_sims,
                                  theta_noninferior = eta)
          grid_search_power=t(apply(sum_0, 1, grid_func_size, results = results_h1,
                                    ti = timing))
          grid_search_power = round(grid_search_power, digits = 3)
          sum_1=sum_0[which(grid_search_power[,1] >= (1-type2)), ]
          dim(sum_1)
          if(dim(sum_1)[1]==0){
            sum_1 = sum_0[which(grid_search_power[,1] == max(grid_search_power[,1])),]
          }
          grid_search_all = t(apply(sum_1, 1, grid_func_size, results = results_h0, ti = timing))
          grid_search_all = round(grid_search_all, digits = 3)
          p_opt = sum_1[which.min(grid_search_all[,5]),]
          p_opt = as.numeric(p_opt)
          p_cut_fut = cut_futility(timing, p_lambda = p_opt[1], p_gamma = p_opt[2])
          p_cut_sup = cut_super(timing, p_lambda = p_opt[1])
          results_stop = rbind(results_stop, cbind(p_cut_fut, p_cut_sup))
          
          PRN0 = grid_search_all[which.min(grid_search_all[,5]),1]
          PET0 = grid_search_all[which.min(grid_search_all[,5]),2]
          PET_EFS0 = grid_search_all[which.min(grid_search_all[,5]),3]
          PET_ESS0 = grid_search_all[which.min(grid_search_all[,5]),4]
          size0 = grid_search_all[which.min(grid_search_all[,5]),5]
          dur0 = grid_search_all[which.min(grid_search_all[,5]),6]
          event0 = grid_search_all[which.min(grid_search_all[,5]),7]
          grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          if(sum(grid_search_power[,1] >= (1-type2))<=1) {
            grid_search_power_min_size = matrix(grid_search_power[which.max(grid_search_power[,1]),], nrow = 1)
          }else{
            grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          }

          PRN1 = grid_search_power_min_size[which.min(grid_search_all[,5]),1]
          PET1 = grid_search_power_min_size[which.min(grid_search_all[,5]),2]
          PET_EFS1 = grid_search_power_min_size[which.min(grid_search_all[,5]),3]
          PET_ESS1 = grid_search_power_min_size[which.min(grid_search_all[,5]),4]
          size1 = grid_search_power_min_size[which.min(grid_search_all[,5]),5]
          dur1 = grid_search_power_min_size[which.min(grid_search_all[,5]),6]
          event1 = grid_search_power_min_size[which.min(grid_search_all[,5]),7]
          xx1=c(PRN0*100, PET0*100, PET_EFS0*100, PET_ESS0*100, size0, event0, dur0)
          xx2=c(PRN1*100, PET1*100, PET_EFS1*100, PET_ESS1*100, size1, event1, dur1)
          results = rbind(results, xx1, xx2)
        }
        if(input$design_stop == "Futility stopping only"){
          
          grid_search = apply(p.tune, 1, grid_func_fut, results = results_h0, ti = timing)
          grid_search = round(grid_search, digits = 3)
          sum_0=p.tune[which(abs(grid_search-type1) < 0.01 & grid_search <= type1), ]
          # dim(sum_0)

          results_h1 = optim_tune(delta_0 = delta_1, p_shape=weibull_shape, m = m_target,
                                  nclt = n_clt, ntrt = n_trt, ti = timing, S0 = surv_3year, t0 = x_year,
                                  ta = ta, tf = tf, theta_0 = theta0, sig2_0 = sigma2_0[h], nsim = n_sims,
                                  theta_noninferior = eta)

          grid_search_power=t(apply(sum_0, 1, grid_func_size_fut, results = results_h1,
                                    ti = timing))
          grid_search_power = round(grid_search_power, digits = 3)
          sum_1=sum_0[which(grid_search_power[,1] >= (1-type2)), ]

          if(dim(sum_1)[1]==0){
            sum_1 = sum_0[which(grid_search_power[,1] == max(grid_search_power[,1])),]
          }

          grid_search_all = t(apply(sum_1, 1, grid_func_size_fut, results = results_h0, ti = timing))
          grid_search_all = round(grid_search_all, digits = 3)
          p_opt = sum_1[which.min(grid_search_all[,3]),]
          p_opt = as.numeric(p_opt)
          p_cut_fut = cut_futility(timing, p_lambda = p_opt[1], p_gamma = p_opt[2])
          
          results_stop = rbind(results_stop, cbind(p_cut_fut))
          
          PRN0 = grid_search_all[which.min(grid_search_all[,3]),1]
          PET0 = grid_search_all[which.min(grid_search_all[,3]),2]
          size0 = grid_search_all[which.min(grid_search_all[,3]),3]
          dur0 = grid_search_all[which.min(grid_search_all[,3]),4]
          event0 = grid_search_all[which.min(grid_search_all[,3]),5]
          grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          if(sum(grid_search_power[,1] >= (1-type2))<=1) {
            grid_search_power_min_size = matrix(grid_search_power[which.max(grid_search_power[,1]),], nrow = 1)
          }else{
            grid_search_power_min_size = grid_search_power[which(grid_search_power[,1] >= (1-type2)), ]
          }

          PRN1 = grid_search_power_min_size[which.min(grid_search_all[,3]),1]
          PET1 = grid_search_power_min_size[which.min(grid_search_all[,3]),2]
          size1 = grid_search_power_min_size[which.min(grid_search_all[,3]),3]
          dur1 = grid_search_power_min_size[which.min(grid_search_all[,3]),4]
          event1 = grid_search_power_min_size[which.min(grid_search_all[,3]),5]
          xx1=c(PRN0*100, PET0*100, size0, event0, dur0)
          xx2=c(PRN1*100, PET1*100, size1, event1, dur1)
          results = rbind(results, xx1, xx2)
        }
      }
      if(input$design_stop == "Futility stopping only"){
        df_final_ops = data.frame(cbind(rep(1:length(sigma2_0), each = 2), weibull_shape,
                                    rep(sigma2_0, each = 2), rep(c(1, delta_1), length(sigma2_0)),
                                    results))
        names(df_final_ops) = c('Design', 'Weibull Shape', 'Normal prior variance',
                            'Delta', 'PRN (%)', 'PET (%)',
                            'ES', 'EE', 'Duration')
        results_stop = round(results_stop, digits = 2)
        df_final_stop = data.frame(cbind(rep(1:length(sigma2_0), each = length(timing)),
                                         rep(sigma2_0, each = length(timing)), 
                                         rep(timing, length(sigma2_0)),
                                         results_stop))
        names(df_final_stop) = c('Design', 'Normal prior variance',
                                 'Information time', 'Futility stopping if posterior prob. <')
      }
      if(input$design_stop == "Futility and superiority stopping"){
        df_final_ops = data.frame(cbind(rep(1:length(sigma2_0), each = 2), weibull_shape,
                                    rep(sigma2_0, each = 2), rep(c(1, delta_1), length(sigma2_0)),
                                    results))
        names(df_final_ops) = c('Design', 'Weibull Shape', 'Normal prior variance',
                            'Delta', 'PRN (%)', 'PET (%)', 'PET futility (%)',
                            'PET superiority (%)', 'ES', 'EE', 'Duration')
        results_stop = round(results_stop, digits = 2)
        df_final_stop = data.frame(cbind(rep(1:length(sigma2_0), each = length(timing)),
                                         rep(sigma2_0, each = length(timing)), 
                                         rep(timing, length(sigma2_0)),
                                         results_stop))
        names(df_final_stop) = c('Design', 'Normal prior variance',
                                'Information time', 'Futility stopping if posterior prob. <', 
                                'Superiority stopping if posterior prob. >=')
      }
    }
    
    return(list(df_final_ops, df_final_stop))

  })
  
  output$out_table_ops <- renderDataTable({
    df_out()[[1]]
  })
  output$out_table_stop <- renderDataTable({
    df_out()[[2]]
  })
}