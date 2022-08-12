# #1a)
# estim_gama_params_A=function(entry_data){
#   entry_data=rgamma(100,shape=2,scale=1)
#   sum_data=sum(entry_data)
#   sum_ln_data=sum(log(entry_data))
#   n=length(entry_data)
#   ll_gamma=(n*alpha*ln(alpha))-(beta*sum_data)+((alpha-1)sum_ln_data)-n*ln(gamma(alpha))
#   optim(2,ll_gamma,method = 'L-BFGS-B',
#         lower=0.00001)
#   
# }
# 
# estim_gama_params_A(c(sample(1:50,20)))
# 
# #b)
# estim_gama_params_B=function(entry_data){
#   n=length(entry_data)
#   beta_hat=alpha_hat/mean(entry_data)
# }

ll_func = function(density_func, sample, params_names){
  
  # Parêmtros: density_func - função da densidade
  #            sample       - vetor da amostra
  #            params_names - vetor com os nomes dos parâmetros a serem 
  #                           otimizados da função da densidade
  
  # Retorno: função de máxima verossimilhança onde o parâmetro será um vetor 
  #          passado na mesma ordem do vetor params_names
  
  function(params){
    
    # Definindo a lista dos parâmetros
    params_list = as.list(params)
    # Definindo os nomes dos parâmetros
    names(params_list) = params_names
    # Definindo a amostra na lista de parâmetros
    params_list$x = sample
    # Calculando a funçãp de verossimilança
    sum(log(call(density_func, params_list)))
  }
}

sample_example=c(22.0,23.9,20.9,23.8,25,24,21.7,23.8,22.8,23.1,23.1,23.5,23,23)
sample_example=rgamma(100,shape=2,rate = 2)


ll_gamma = likelihood_function(dgamma, sample_example,params = c("shape", "rate"))

x_bar = mean(sample_example)
# Definindo chute inicial para beta
(beta0 = x_bar/var(sample_example))
alpha0 = x_bar*beta0


estimated_params = optim(c(0,0), ll_gamma, method = c("L-BFGS-B"),
                     control=list(fnscale=-1), lower = 0.01)







##### parte 2
iterative_cauchy=function(data,error_choice=1e-6){
  
  teta_estim=mean(data)
  s_teta=function(teta_hat) sum(2*(data-teta_hat)/(1+((data - teta_hat)**2)))
  h_teta=function(teta_hat) 2*sum((1-((data-teta_hat)**2))/((1+((data - teta_hat)**2))**2))
  
  error_eval=1
  i=1
  while(error_eval>error_choice)
    {
      teta_estim[i+1]=teta_estim[i]+(s_teta(teta_estim[i])/h_teta(teta_estim[i]))
      error_eval=teta_estim[i+1]-teta_estim[i]
      print(teta_estim)
        i=i+1
    }
  
   return(teta_estim)
  }
teste=rcauchy(100,5,1)  
mean(teste)
iterative_cauchy(teste)



gera_bolfarine=function(teta){
  
  ##Armazenando a função alva
  func_alva=function(x) (1+teta*x)/2
  
  ##Armazenando a função candidata
  func_candidata=function(x)dunif(x,min=-1,max=1)
  
  alva=vector()
  M <- optimize(f = function(x)func_alva(x)/func_candidata(x), interval = c(-1,1), maximum = T)$objective
  
  #Utilizando o metodo da aceitação rejeicao para a criacao dos valores aleatorios
  j=1
  for(i in 1:1000){
    u=runif(1)
    y=runif(1,-1,1)
    if(u*M<=func_alva(y)/func_candidata(y)){
      alva[j]=  y
      j=j+1
    }
  }
  ##calculando a proporcao de pontos aceitos
  proporcao=(j/1000)*100
  hist(alva,freq=F,col='red',xlim=c(-1,1))
  curve(func_alva, add = T,col='black',lwd=2)
  return(alva)
}

sample_bolfarine=gera_bolfarine(0.5)

iterative_bolfarine=function(data,error_choice=1e-6){
  
  teta_estim=mean(data)
  s_teta=function(teta_hat) sum(data/(1+(teta_hat*data)))
  h_teta=function(teta_hat) sum((data**2)/((1+teta_hat*data)**2))
  
  error_eval=1
  i=1
  while(error_eval>error_choice)
  {
    teta_estim[i+1]=teta_estim[i]+(s_teta(teta_estim[i])/h_teta(teta_estim[i]))
    error_eval=teta_estim[i+1]-teta_estim[i]
    print(teta_estim)
    i=i+1
  }
  
  return(teta_estim)
}

iterative_bolfarine(sample_bolfarine)


