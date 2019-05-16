nStrata <- function(data,stratanames,alpha,moe,S2=NULL,pq=NULL,N=Inf,
                    method=c("none","prop","optimum")){
  if(missing(method)) {warning("the method is not specified; by default, the method is none.")
    method="none"
  }
  if(!(method %in% c("none","prop","optimum"))){ 
    stop("the method name is not in the list")}
  if(method=="optimum" && N==Inf){
    stop("You have to input N if method is optimum.")
  }
  dados=as.data.frame(data)
  m=match(stratanames,colnames(dados))
  dados$stratum=do.call(paste, c(dados[m], sep = "_"))
  n_strata=length(unique(dados$stratum))
  if(!is.null(S2)){var=S2
  if(length(S2)==1){var=rep(S2,length(n_strata))}
  }else{var=pq
  if(length(pq)==1){var=rep(pq,length(n_strata))}}
  Wh=table(dados$stratum)/sum(table(dados$stratum))
  Wh2var2=Wh^2*var
  wh=Wh *sqrt(var) /sum(Wh *sqrt(var))
  z_alpha=qnorm(1 - ((1 - alpha) / 2))
  if(method=="none"){
    n_linha=(z_alpha^2/moe^2)*sum(Wh2var2/wh)
    if(is.infinite(N)){
      return(n_linha)
    }else{
      return(ceiling(n_linha/(1+(z_alpha^2/(moe^2 * N))*sum(Wh2var2))))
    }
  }
  if(method=="prop"){
    n_linha=(z_alpha^2 *sum(Wh2var2))/moe^2
    if(is.infinite(N)){
      return(n_linha)
    }else{
      return(ceiling(n_linha/(1+(n_linha/N))))
    }
    
  }
  if(method=="optimum"){
    WhSh <- Wh*sqrt(var)
    return(ceiling(sum(WhSh)^2/(1+(z_alpha^2/(N*moe^2))*sum(Wh2var2))))
  }
}
