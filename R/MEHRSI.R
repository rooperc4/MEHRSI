#' A function to compile the right side of a habitat equation
#'
#' This function compiles text of equation for model. It takes the habitat variables
#' and makes them into a text version of the right side of an equation that is then read into the
#' habitat model function.
#' @param variables Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  starting values for each of the parameters 
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' hab_equation()

hab_equation<-function(variables,form,par)
  
  
{
  n<-0
  part<-form
  
  
  for(i in 1:length(variables))
  {n<-sum(form[1:i])
  if(form[i]==1)
    part[i]<-paste(par[n],"*",variables[i])
  
  if(form[i]==2)
    part[i]<-paste(par[n-1],"*",variables[i],"*exp(-",par[n],"*",variables[i],")")
  
  
  if(form[i]==3)
    part[i]<-paste(par[n-2],"+",par[n-1],"*",variables[i],"+",par[n],"*",variables[i],"^2")
  if(form[i]==0)
    part[i]<-0
  }
  
  
  partial_eq<-paste(part,collapse="+")
  
  return(partial_eq)
}


#' A function to predict CPUE from a habitat equation
#'
#' This function uses the habitat model equation and parameter values from
#' the hab_equation function to calculate a predicted CPUE. A dummy variable
#' for a year effect is added from the yearp function using the year of the 
#' survey for each of the stations. The function also #' can incorporate a 
#' presence-absence prediction (defaults to 1 for present) #' from an external 
#' step (usually from the presence_absence function)
#' @param variables Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  values for each of the parameters 
#' @param pred_pa the predicted presence or absence at each transect
#' @param year year of the data collection
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' HabModel()

HabModel<-function(variables, form, par,pred_pa=1,year)
{
  equation1<-hab_equation(variables,form,par)
  
  PCPUE<-(eval(parse(text = equation1))+yearp(year,par,form))*pred_pa
  
  return(PCPUE)
  
}


#' A function to create a dummy year variable for a habitat equation
#'
#' This function creates a dummy variable for a year effect that can be
#' added to the habitat model equation. 
#' @param year Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  values for each of the parameters 
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' yearp()

yearp<-function(year,par,form){
  
  rowpy<-as.integer(by(year,year,nrow))
  p1<-length(par)
  p2<-sum(form)+1
  pary<-par[p2:p1]
  yearps<-rep(pary,rowpy)
  
  return(yearps)}

#' This function calculates negative log-liklihood for a habitat model
#'
#' This function takes observed and predicted CPUE and calculates a negative log-likelhood for the data
#' it uses either a log-normal or gamma distribution. The model form, parameters, habitat
#' variables, presence or absence and year are passed through this function to the HabModel function which computes the
#' predicted values. It accepts "normal" or "gamma" as arguments for distribution.
#' @param year Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  values for each of the parameters 
#' @param variables Set of habitat variables to use in the modeling
#' @param ob_CPUE Observed CPUE at the transect
#' @param pred_pa Predicted presence or absence at the transect
#' @param distribution Either "normal" (the default) or "gamma"
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' modlike()

modlike<-function(par, form, variables,ob_CPUE,pred_pa=1,year,distribution="normal")
{
 if(distribution=="normal"){
   PCPUE<-HabModel(variables,form,par,pred_pa,year)
  sigma<-sd(ob_CPUE)
  likj<-ob_CPUE
  likj<-log(sigma)+0.5*log(2*3.141593)+(ob_CPUE-PCPUE)^2/(2*sigma^2)
  lik_sum=sum(likj)
  return(lik_sum)}
 if(distribution=="gamma"){
   PCPUE<-HabModel(variables,form,par,pred_pa,year)
	gresid<-(ob_CPUE1-unlist(PCPUE))^2
lik_sum<--sum(dgamma(unlist(gresid),1, 1, log = TRUE))
return(lik_sum)}
  
}


#' A wrapper function to find the best fitting MEHRSI habitat model
#'
#' This function takes the habitat variables, observed CPUE, the equation form, starting 
#' parameters, the predicted presence or absence, the year, and the distribution. The
#' function fits an intital full model and then iteratively reduces the number of
#' parameters for each habitat variable. It then chooses the best-fitting reduced model
#' and repeats the process until there is no improvement in AIC with the removal of 
#' additional parameters or variables. The function then returns a list with the fitting
#' history and the best model results.
#' @param year Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  values for each of the parameters 
#' @param variables Set of habitat variables to use in the modeling
#' @param ob_CPUE Observed CPUE at the transect
#' @param pred_pa Predicted presence or absence at the transect
#' @param distribution Either "normal" (the default) or "gamma"
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' iteration()

iteration<-function(variables,form,ob_CPUE,pred_pa=1,par,year,distribution="normal"){

  #for testing
#  variables<-AI_data[,(var1)]
#  form<-forms1
#  ob_CPUE<-ob_CPUE1
#  pred_pa<-pred_pa1
#  par<-par1
#  year<-year1
#CPUE data plotted against habitat variables 
  png(filename="Raw_data.png",width=7,height=7,res=300,units="in") 
  par(mfrow=c(ceiling(length(variables)/3),3)) 
  for(i in 1:length(variables)){
    xdata<-cbind(eval(parse(text=variables[i])),ob_CPUE)
    plot(xdata,ylab="ob_CPUE",xlab=colnames(variables[i]))
  } 
  dev.off()
  
   iter<-1
  fit<-rep(NA,sum(form)*length(par)+30)
  AIC<-rep(NA,sum(form)*length(par)+30)
  form_next<-array(,dim=c(length(AIC)+30,length(variables)))
  AICm<-1
  AIC[AICm]<-0
  AICf<-1
  
  while(AICf>AIC[AICm]){
    n<-1
    yearn<-length(unlist(unique(year)))
    
    #run full model at beginning of each round
   par<-rep(0,(sum(form)+yearn))
#    fit1<-nlm(modlike,par,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal",ndigit=12, gradtol=.0000015, stepmax=3, steptol=.000001, iterlim=1000,print.level=2,hessian=TRUE)
    fit1<-nlminb(par,modlike,gradient=NULL,hessian=NULL,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal",control=list(trace=1))
    fit[iter]<-fit1$objective
#print(fit[iter])    
    AIC[iter]<-2*fit[iter]+2*sum(form)
    AICf<-2*fit[iter]+2*sum(form)
    form_next[iter,]<-form 
    
    #repeat fitting for every variable
    for(k in 1:dim(variables)[2]){
      nadd<-form[n]
      #fit model for each removal of parameter
      if(form[n]>0) {
        for(i in 1:form[n]){
          
          form[n]<-form[n]-1
          par<-rep(0,(sum(form)+yearn))
          iter<-iter+1
        fit1<-nlminb(par,modlike,gradient=NULL,hessian=NULL,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal")
#          print(i) 
        #fit1<-nlm(modlike,par,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal",ndigit=100, gradtol=.0000015, stepmax=3, steptol=.000001, iterlim=1000,print.level=2,hessian=TRUE)
        #  fit1<-optim(par,modlike,gr=NULL,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal",method="BFGS",control=list(trace=1))
          fit[iter]<-fit1$objective
          
          AIC[iter]<-2*fit[iter]+2*sum(form)
          form_next[iter,]<-form 
 print(AIC[1:20]) 
 print(form_next[1:20,])
          
        }}
      form[n]<-form[n]+nadd
      n<-n+1
      
      print(k)
   #   browser()
    }
    
    #estimate best AIC for each round
    AICm<-which.min(AIC)
    
    #use best AIC to determine the full model for the next round
    form<-form_next[AICm,]
    iter<-iter+1
    print(form)
  #  print(date())
    flush.console()
  }
  
  #need to fit model at end and extract residuals, do a plot, calculate r2, etc
  par<-rep(0,(sum(form)+yearn))
  form_hist<-form_next
  best_fit<-fit[AICm]
  best_AIC<-AIC[AICm]
  best_model<-form_next[AICm,]
#  fit2<-nlm(modlike,par,best_model,variables,ob_CPUE,pred_pa=1,year,distribution="normal",ndigit=100, gradtol=.0000015, stepmax=3, steptol=.000001, iterlim=1000)
  fit2<-nlminb(par,modlike,gradient=NULL,hessian=NULL,best_model,variables,ob_CPUE,pred_pa=1,year,distribution="normal",control=list(trace=1))
  best_loglik<-fit2$objective
  best_parameters<-fit2$par
  best_PCPUE<-HabModel(variables,best_model,best_parameters,pred_pa,year)
  resids<-ob_CPUE-best_PCPUE
  r_squared<-cor(ob_CPUE,best_PCPUE)^2
  
 #######################################################################
  #FIGURES
  
   #Variable relationships  
  png(filename="VariableRelationships.png",width=7,height=7,res=300,units="in")
  par(mfrow=c(ceiling(length(variables)/3),3))
  n<-1
  for(i in 1:length(variables)){
    maxdata<-max(eval(parse(text=variables[i])))
    mindata<-min(eval(parse(text=variables[i])))
    xdata<-seq(mindata,maxdata,by=(maxdata-mindata)/50)
    pars<-c(unlist(fit2$par))
    form1<-best_model
    if(form1[i]==1) {
      part<-paste(round(pars[n],4),"*x")
      ydata<-pars[n]*xdata  }
    if(form1[i]==2) {
      part<-paste(round(pars[n],4),"*x*exp(-",round(pars[n+1],4),"*x)")
      ydata<-pars[n]*xdata*exp(-pars[n+1]*xdata) }
    if(form1[i]==3)   {
      part<-paste(round(pars[n],4),"+",round(pars[n+1],4),"*x+",round(pars[n+2],4),"*x^2")
      ydata<-pars[n]+pars[n+1]*xdata+pars[n+2]*xdata^2 }
    if(form1[i]==0)  {
      part<-"Not significant"
      ydata<-1         }
    n=n+form1[i]
    zdata<-cbind(xdata,ydata)
    plot(zdata,xlab=colnames(variables[i]),main=part,ylab="Partial Effect on CPUE")}  
  dev.off()   
  
  #Observed v. predicted
  png(filename="ObservedPredicted.png",width=7,height=7,res=300,units="in")
  par(mfrow=c(2,2))
  plot(ob_CPUE,best_PCPUE, xlab="Observed", ylab="Predicted");
  regr<-lm(best_PCPUE~ob_CPUE);
  abline(regr);
  title("observed and predicted values")
  mtext(paste("r_squared = ",round(unlist(r_squared),3)),side=3)
  dev.off()
  
#QQPLOT
  png(filename="QQplot.png",width=7,height=7,res=300,units="in")
  qqdata1<-cbind(pred_pa,resids)
  qqdata1<-subset(qqdata1,pred_pa>0)
  qqnorm((qqdata1[,"resids"]), main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  qqline((qqdata1[,"resids"]))
    dev.off()       
    
#TABLES
    tabledata<-array(dim=c(2,5))
    colnames(tabledata)<-c("Model","Number of parameters","AIC","R2","Contribution")
    par<-rep(0,(sum(form)+yearn))
    fitf<-nlminb(par,modlike,gradient=NULL,hessian=NULL,form,variables,ob_CPUE,pred_pa=1,year,distribution="normal",control=list(trace=1))
    tabledata[1,1]<-"Full model"
    tabledata[1,2]<-sum(form)
    tabledata[1,3]<-AIC[1]
    full_PCPUE<-HabModel(variables,form,fitf$par,pred_pa,year)
    tabledata[1,4]<-cor(ob_CPUE,full_PCPUE)^2
    tabledata[1,5]<-"--"
    
    tabledata[2,1]<-"Best model"
    tabledata[2,2]<-sum(best_model)
    tabledata[2,3]<-best_AIC
    tabledata[2,4]<-r_squared
    tabledata[2,5]<-"--"
    
    contribf<-array(dim=c(dim(variables)[2],5))
    colnames(contribf)<-c("Model","Number of parameters","AIC","R2","Contribution")
        contribf[,1]<-colnames(variables)
    for(i in 1:dim(contribf)[1]){
    formt<-form1
    formt[i]<-0
    par<-rep(0,(sum(formt)+yearn))
    fitf<-nlminb(par,modlike,gradient=NULL,hessian=NULL,formt,variables,ob_CPUE,pred_pa=1,year,distribution="normal",control=list(trace=1))
    contribf[i,2]<-sum(formt)
    contribf[i,3]<-2*fitf$objective+2*sum(formt)
    part_PCPUE<-HabModel(variables,formt,fitf$par,pred_pa,year)
    contribf[i,4]<-cor(ob_CPUE,part_PCPUE)^2
    contribf[i,5]<-1-best_loglik/fitf$objective}
    contribf<-contribf[form1!=0,]
    contribf<-contribf[order(contribf[,5],decreasing=TRUE),]
    contribf[,5]<-as.numeric(contribf[,5])/max(as.numeric(contribf[,5]))
    tabledata<-rbind(tabledata,contribf)
    pander::pandoc.table(tabledata)
    
  return(list(form_hist=form_hist,best_fit=best_fit,best_AIC=best_AIC,best_loglik=best_loglik,r_squared=r_squared,best_parameters=best_parameters,best_model=best_model,best_PCPUE=best_PCPUE,resids=resids))
}

#' This function calculates cdf of CPUE for predicting presence or absence in a habitat model
#'
#' This function calculates the cumulative distribution of the CPUE of a species weighted across a variable.
#' It computes the upper 95% and lower 5% boundaries (default) for the distribution of a species over a variable.
#' @param variable Set of habitat variables used to predict presence or absence
#' @param CPUE Observed CPUE at the transect
#' @param cutoff level at which you want to set the limits, defaults to 0.05
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' cdf_bound()

cdf_bound<-function(variable, CPUE,cutoff){
  
  cdfu<-0
  ord<-0
  cdfl<-0
  
  vars<-c(unlist(variable))
  weights<-rep(0,length(variable))
  
  sumwt<- sum(CPUE)
  
  for(i in 1:length(vars)){
    ord<-which.min(vars)
    
    weights[i] = CPUE[ord]
    
    if(cdfu<=1-cutoff){ 
      cdfu<-cdfu + weights[i]/sumwt
      cdf_95<-variable[ord]
    }
    
    if(cdfl<=cutoff){
      cdfl<-cdfl + weights[i]/sumwt
      cdf_05<-variable[ord]
    }
    
    
    
    vars[ord]<-9999
    
  }
  return(list(upper=cdf_95,lower=cdf_05))
}


#' This function calculates error for cdf of CPUE for predicting presence or absence in a habitat model
#'
#' This function estimates the error around the upper 95% and lower 5% boundaries f
#' or the distribution of a species over a variable by bootstraping the upper and lower 
#' boundaries. The bootstrap errors are calculated by resampling the data 1000 times and 
#' then computing the limits. This is done for each of the variables used to determine presence or absence.
#' @param variables Haibtat variables used for estimating presence or absence
#' @param CPUE Observed CPUE at the transect
#' @param breps Number of bootstraps to calculate error estimate
#' @param cutoff level at which you want to set the limits, defaults to 0.05
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' cdf_limits()

cdf_limits<-function(CPUE,variables,breps=1000,cutoff=.05){

  pajpopdata<-cbind(CPUE,variables)
  dimt<-dim(variables)
  cols<-dimt[2]
  
  depth_limitu<-array(,dim=c(breps,cols))
  depth_limitl<-array(,dim=c(breps,cols))
  cdf_upper<-c(rep(0,cols))
  cdf_lower<-c(rep(0,cols))
  cdf_upper_SD<-c(rep(0,cols))
  cdf_lower_SD<-c(rep(0,cols))
  
  for(i in 1:breps){
    
    pabootdata <- pajpopdata[sample(1:nrow(pajpopdata), replace=TRUE),]
    
    for(j in 1:cols){
      depth_limit<-cdf_bound(pabootdata[,j+1],pabootdata[,1],cutoff)
      depth_limitu[i,j]<-depth_limit$upper
      depth_limitl[i,j]<-depth_limit$lower
    }
  }
  for(i in 1:cols){
    cdf_upper[i]<-mean(depth_limitu[,i])
    cdf_upper_SD[i]<-sd(depth_limitu[,i])
    cdf_lower[i]<-mean(depth_limitl[,i])
    cdf_lower_SD[i]<-sd(depth_limitl[,i])
  }
  nam<-names(variables)
  out<-cbind(cdf_lower,cdf_upper,cdf_lower_SD,cdf_upper_SD)
  row.names(out)<-colnames(variables)
  return(limits=out)
}


#' This function estimates presence or absence at each transect
#'
#' This function uses the cdf_limits function to calculate the limits of presence or absence 
#' for a species based on a variable of interest. For example, if the cdf_limits found that 
#' a species was only present at depths between 100 m and 200 m, this function would input 
#' the depths from a trawl survey haul and determine presence or absence or absence of the 
#' species at that site.
#' @param limits Limits on the distribution of a species from cdf_limits function
#' @param variables Variables used for determining presence or absence
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' presence_absence()

presence_absence<-function(variables,limits,ob_CPUE){
library(PresenceAbsence)
library(pander)
    dime<-dim(variables)
  cols2<-dime[2]
  rows2<-dime[1]
  predpa<-c(rep(1,rows2))
  
  
  
  for(i in 1:rows2){	
    for(j in 1:cols2){
      if(variables[i,j]<limits[j,1]){predpa[i]=0}
      if(variables[i,j]>limits[j,2]){predpa[i]=0}
    }
  }
  out<-c(predpa)
  
  #Table 2
  obspres<-ifelse(ob_CPUE>0,1,0)
  aucdata<-data.frame(cbind(seq(1,length(ob_CPUE),1),obspres,predpa))
  pander::pandoc.table(cmx(aucdata), caption="Presence and Absence",row.names=c("Predicted Present","Predicted Absent"),col.names=c("Present","Absent"))
  
  return(out)
}

#' This function estimates distance of one set of points to another
#'
#' This function takes two lat-long pairs and calculates the distance between
#' them in m or km
#' @param lat1 Latitude of point 1
#' @param lat2 Latitude of point 2
#' @param long1 Longitude of point 1
#' @param long2 Longitude of point 2
#' @param unit either "km" or "m"
#' @keywords distance
#' @export
#' @examples
#' dist_xy()

dist_xy<-function(lat1,long1,lat2,long2,unit){
  #lat1<-56
  #lat2<-57
  #long1<-(-174)
  #long2<-(-171)
  #unit="m"
  
  
  r<-ifelse(unit == "m", 6371200,6371.2)
  
  la1 = lat1 * pi / 180
  la2 = lat2 *pi / 180
  lo1 = long1 * pi / 180
  lo2 = long2 * pi/ 180
  
  
  rads <- acos((sin(la1) * sin(la2)) + (cos(la1) * cos(la2) * cos(lo2 - lo1)))
  return(r * rads)
}

#' Function to estimate boot-strapped errors for annual index
#'
#' This function takes the best fitting model and resamples the data refitting and recomputing
#' the indices to come up with a bootstrapped error estimate for the annual index value.
#' @param year Set of habitat variables to use in the modeling
#' @param form for each variable, a form for the model equation (1, 2 or 3) 
#' @param par  values for each of the parameters 
#' @param variables Set of habitat variables to use in the modeling
#' @param ob_CPUE Observed CPUE at the transect
#' @param pred_pa Predicted presence or absence at the transect
#' @param distribution Either "normal" (the default) or "gamma"
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' dist_xy()

boot_survey_error<-function(best_model,boot_reps=500,year,variables,pred_pa=1,ob_CPUE,distribution="normal"){

yearn<-length(unlist(unique(year)))
parameter_ests<-array(0,dim=c(boot_reps,(sum(best_model$best_model)+yearn)))
likes_boot<-array(0,boot_reps)

for(i in 1:boot_reps){	
  bootdata<- sample(1:length(pred_pa), replace=TRUE)
  par1<-rep(0,(sum(best_model$best_model)+yearn))
  yearu<-data.frame(year=year[bootdata,])
  variablesu<-variables[bootdata,]
  ob_CPUEu<-ob_CPUE[bootdata]
  formu<-best_model$best_model
  pred_pau<-pred_pa[bootdata]
#  fit3<-optim(par1,modlike,gr=NULL,formu,variablesu,ob_CPUEu,pred_pau,yearu,distribution="normal",method="BFGS",control=list(trace=1))
  fit3<-nlminb(par1,modlike,gradient=NULL,formu,variablesu,ob_CPUEu,pred_pau,yearu,distribution="normal",control=list(trace=1))
  best_parameters<-fit3$par
  parameter_ests[i,]<-fit3$par
  likes_boot[i]<-fit3$objective
  print(i)
#  print(date())
  flush.console()
}

return(list(likes_boot,parameter_ests))
}


#' Function to estimate annual index
#'
#' This function takes the best fitting model and resamples the data refitting and recomputing
#' the indices to come up with a bootstrapped error estimate for the annual index value.
#' @param parameter_ests Bootstrap estimates of parameters
#' @param best_model for each variable, a form for the model equation (1, 2 or 3) 
#' @param best_parameters  fitted values for each of the parameters 
#' @param variables Set of habitat variables to use in the modeling
#' @param ob_CPUE Observed CPUE at the transect
#' @param pred_pa Predicted presence or absence at the transect
#' @param year years of the survey
#' @param minob 1/2 of the minimum CPUE observation for backtransformation of the data
#' @boot_reps number of bootstrap replicates
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' index_calc()

index_calc<-function(parameter_ests,best_model,best_parameters,variables,year,minob,boot_reps){
library(ggplot2)
parameter_error<-parameter_ests
year_cols<-sum(best_model)+1
yearn<-length(unlist(unique(year)))
var2<-array(dim=c(0,dim(variables)[2]))
#var2<-best_model$best_model
for(i in 1:dim(variables)[2]){
var2[i]<- median(variables[[i]])}
model_equation<-hab_equation(var2,best_model,best_parameters)
overall_mean<-eval(parse(text = model_equation))
index_est<-exp(best_parameters[year_cols:(year_cols+yearn-1)]+overall_mean)-minob
index_est_sd<-apply(exp(parameter_error[,year_cols:(year_cols+yearn-1)]-minob),2,FUN=sd)
index_est_cv<-100*index_est_sd/index_est
index_table<-data.frame(Year=unlist(unique(year),use.names=FALSE),Index=index_est,SD=index_est_sd,CV=index_est_cv)
index_table<-index_table[order(index_table$Year),]
pander::pandoc.table(index_table,row.names=FALSE,digits=4)
#PLOT OF INDEX

p<-ggplot(index_table,aes(x=Year,y=Index))+geom_line()+geom_point()+
  geom_ribbon(aes(ymin=Index-SD*1.96/sqrt(boot_reps), ymax=Index+SD*1.96/sqrt(boot_reps)),
              alpha=0.2)
png(filename="IndexPlot.png",width=7,height=7,units="in",res=300)
    print(p)
dev.off()


return(index_table)

}

#' Function to spatial autocorrelation in residuals
#'
#' This function takes the best fitting model and resamples the data refitting and recomputing
#' the indices to come up with a bootstrapped error estimate for the annual index value.
#' @param longitude Longitude of haul
#' @param latitude latitude of haul
#' @param residuals  residuals of the best fitting model 
#' @param pred_CPUE predicted CPUE of the best fitting model
#' @param ob_CPUE Observed CPUE at the transect
#' @param region either "GOA" (default) or "AI"
#' @keywords habitat model, survey abundance index
#' @export
#' @examples
#' spatial_resids()

spatial_resids<-function(longitude,latitude,residuals,pred_CPUE,ob_CPUE,region="GOA"){					
#latitude<-Juvenile_POP_Data$lat
#longitude<-Juvenile_POP_Data$long
#residuals<-jpopdata_out$resids
#pred_CPUE<-jpopdata_out$best_PCPUE
#ob_CPUE<-ob_CPUE1
  library(spatial)
  library(MASS)	
 # spatial_data<-data.frame(x=longitude,y=latitude,z=residuals)
  png(filename="spatial_resids.png",width=7,height=7,units="in",res=300)
  par(mfrow=c(3,2))
  dist1<-seq(0,10,.1)
  dist2<-seq(0,40,1)
  r_sq_spatial.kr<-array(0,dim=c(100,2))
  range1<-.01
  for(i in 1:100){
    jpop.kr<-surf.gls(2,sphercov,x=longitude,y=latitude,z=residuals,r=dist1, d=range1)
    jpop.kr.pred<-predict.trls(jpop.kr,longitude,latitude)
    
    jpop.pred<-jpop.kr.pred+pred_CPUE
    r_sq_spatial.kr[i,1]<-range1
    r_sq_spatial.kr[i,2]<-(cor(jpop.pred,ob_CPUE))^2
    range1<-range1+.01
  }
  
  range2i<-which.max(r_sq_spatial.kr[,2])
  range2<-r_sq_spatial.kr[range2i,1]
  jpop.kr<-surf.gls(2,sphercov,x=longitude,y=latitude,z=residuals,r=dist1, d=range2)
  jpop.kr.pred<-predict.trls(jpop.kr,longitude,latitude)
  jpop.pred<-jpop.kr.pred+pred_CPUE
  r_sq_spatial<-(cor(jpop.pred,ob_CPUE))^2
  trsurf <- trmat(jpop.kr,  -169, -133, 52, 61,5)
  eqscplot(trsurf, type = "n")
  contour(trsurf, add = TRUE)
if(region=="GOA"){  prsurf <- prmat(jpop.kr, -169, -133, 52, 61, 5)}
if(region=="AI"){  prsurf <- prmat(jpop.kr, -190, -165, 52, 61, 5)}
  contour(prsurf, levels=seq(-1, 1, .1))
  variogram(jpop.kr,100)
  correlogram(jpop.kr,100)
  lines(dist2,sphercov(dist2,range2))
  plot(r_sq_spatial.kr[,1],r_sq_spatial.kr[,2])
  plot(ob_CPUE,jpop.pred)
  title(round(r_sq_spatial,3))
  
  dev.off()
}
