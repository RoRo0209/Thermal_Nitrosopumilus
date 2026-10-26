# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(.release_root)

args<-commandArgs(trailingOnly=TRUE)
root <- .release_root
out<-file.path(root,'60weeks','dense_temperature_reanalysis_20261015')
d<-read.csv(file.path(root,'60weeks','lag_reanalysis_20261015','source_observations.csv'))
env<-read.csv(file.path(out,'environment_temperature.csv'))
d$station<-factor(d$station,levels=c('E','W','Z'))
stopifnot(nrow(d)==39,all(table(d$station)==13),nrow(env)==165,
 !anyDuplicated(paste(env$station,env$week)),all(is.finite(env$temperature)))
envkey<-paste(env$station,env$week)
stopifnot(max(abs(env$temperature[match(paste(d$station,d$week),envkey)]-d$temperature))<1e-12)
d$annual_sin<-sin(2*pi*d$week/52);d$annual_cos<-cos(2*pi*d$week/52)
d$log_ratio<-log10(d$HighT/d$LowT)
d$log_ratio_eps1e6<-log10((d$HighT+1e-6)/(d$LowT+1e-6))
lags<--8:8; joined<-list()
for(L in lags){
 z<-d;z$lag_weeks<-L;z$temperature_week<-z$week-L
 ix<-match(paste(z$station,z$temperature_week),envkey)
 z$temp_lag<-env$temperature[ix];z$temperature_sample<-env$sample_id[ix]
 joined[[as.character(L)]]<-z
}
full<-do.call(rbind,joined);rownames(full)<-NULL
write.csv(full,file.path(out,'all_explicit_matches.csv'),row.names=FALSE)
common<-function(ls) Reduce(intersect,lapply(as.character(ls),function(k)joined[[k]]$sample[is.finite(joined[[k]]$temp_lag)]))
ids04<-common(0:4);ids08<-common(0:8)
write.csv(data.frame(sample=ids04),file.path(out,'common_0_to_4_samples.csv'),row.names=FALSE)
write.csv(data.frame(sample=ids08),file.path(out,'common_0_to_8_samples.csv'),row.names=FALSE)
coverage<-aggregate(is.finite(temp_lag)~station+lag_weeks,full,sum)
names(coverage)[3]<-'n_matched';write.csv(coverage,file.path(out,'pair_counts.csv'),row.names=FALSE)
pear<-function(x,y)if(length(x)>=4&&sd(x)>1e-12&&sd(y)>1e-12)cor(x,y)else NA_real_
describe<-function(z,response,cohort,adjustment){
 ans<-list()
 for(s in c('E','W','Z','pooled_station_adjusted')){
  g<-if(s=='pooled_station_adjusted')z else z[z$station==s,]
  g<-g[is.finite(g$temp_lag),];g$y<-g[[response]]
  covars<-c(if(s=='pooled_station_adjusted')'station',if(adjustment=='annual_harmonic')c('annual_sin','annual_cos'))
  if(nrow(g)<max(4,length(covars)+3))next
  xf<-if(length(covars))reformulate(covars,response='temp_lag')else temp_lag~1
  yf<-if(length(covars))reformulate(covars,response='y')else y~1
  x<-resid(lm(xf,g));y<-resid(lm(yf,g))
  ans[[length(ans)+1]]<-data.frame(cohort=cohort,station=s,response=response,adjustment=adjustment,
    lag_weeks=unique(g$lag_weeks),n=nrow(g),pearson_r=pear(x,y),
    spearman_rho=if(s!='pooled_station_adjusted'&&adjustment=='unadjusted')cor(g$temp_lag,g$y,method='spearman')else NA_real_)
 }
 do.call(rbind,ans)
}
assoc<-list()
for(cohort in c('lag_specific','common_0_to_4','common_0_to_8')){
 ls<-switch(cohort,lag_specific=lags,common_0_to_4=0:4,common_0_to_8=0:8)
 ids<-switch(cohort,lag_specific=d$sample,common_0_to_4=ids04,common_0_to_8=ids08)
 for(L in ls)for(response in c('log_ratio','log_ratio_eps1e6'))for(adj in c('unadjusted','annual_harmonic')){
  z<-joined[[as.character(L)]];z<-z[z$sample %in% ids,]
  assoc[[length(assoc)+1]]<-describe(z,response,cohort,adj)
 }
}
assoc<-do.call(rbind,assoc);write.csv(assoc,file.path(out,'correlations.csv'),row.names=FALSE)
del<-list()
for(w in sort(unique(d$week[d$sample%in%ids04])))for(L in 0:4)for(adj in c('unadjusted','annual_harmonic')){
 z<-joined[[as.character(L)]];z<-z[z$sample%in%ids04 & z$week!=w,]
 a<-describe(z,'log_ratio','common_0_to_4',adj);a$omitted_week<-w
 del[[length(del)+1]]<-a
}
write.csv(do.call(rbind,del),file.path(out,'leave_one_response_week_out.csv'),row.names=FALSE)
X<-sapply(0:4,function(L){z<-joined[[as.character(L)]];z$temp_lag[match(ids04,z$sample)]})
colnames(X)<-paste0('lag',0:4);write.csv(cor(X),file.path(out,'candidate_temperature_correlations.csv'))

suppressPackageStartupMessages(library(mgcv))
gam_results<-list();gls_results<-list()
for(L in 0:4){
 z<-joined[[as.character(L)]];z<-z[z$sample%in%ids04,];z<-z[order(z$station,z$week),]
 for(response in c('log_ratio','log_ratio_eps1e6'))for(station_adjusted in c(FALSE,TRUE)){
  z$y<-z[[response]]
  fm<-if(station_adjusted)y~station+s(temp_lag,k=4) else y~s(temp_lag,k=4)
  m<-gam(fm,data=z,method='REML');sm<-summary(m)
  gam_results[[length(gam_results)+1]]<-data.frame(lag_weeks=L,n=nrow(z),response=response,
    station_adjusted=station_adjusted,edf=sm$s.table[1,'edf'],nominal_p_unadjusted_for_dependence=sm$s.table[1,'p-value'],
    total_deviance_explained_pct=sm$dev.expl*100,AIC=AIC(m))
 }
 z$y<-z$log_ratio
 for(season in c(FALSE,TRUE)){
  controls<-c('station',if(season)c('annual_sin','annual_cos'))
  fit_gls<-function(addtemp){
   fm<-reformulate(c(controls,if(addtemp)'temp_lag'),response='y')
   fits<-lapply(c(.2,.6,.9),function(start)tryCatch(nlme::gls(fm,data=z,method='ML',
      correlation=nlme::corCAR1(value=start,form=~week|station),
      control=nlme::glsControl(msMaxIter=200)),error=function(e)e))
   ok<-vapply(fits,inherits,logical(1),'gls')
   if(!any(ok))return(list(error=paste(vapply(fits,conditionMessage,character(1)),collapse='; ')))
   fits<-fits[ok];m<-fits[[which.max(vapply(fits,function(f)as.numeric(logLik(f)),numeric(1)))]]
   list(model=m)
  }
  base<-fit_gls(FALSE);fit<-fit_gls(TRUE)
  if(is.null(fit$model)||is.null(base$model)){
   gls_results[[length(gls_results)+1]]<-data.frame(lag_weeks=L,n=nrow(z),season_adjusted=season,
     temperature_beta=NA_real_,model_based_se=NA_real_,rho_per_week=NA_real_,AIC=NA_real_,
     delta_AIC_vs_baseline=NA_real_,error=paste(base$error,fit$error))
  }else{
   m<-fit$model;tab<-summary(m)$tTable
   gls_results[[length(gls_results)+1]]<-data.frame(lag_weeks=L,n=nrow(z),season_adjusted=season,
    temperature_beta=tab['temp_lag','Value'],model_based_se=tab['temp_lag','Std.Error'],
    rho_per_week=as.numeric(coef(m$modelStruct$corStruct,unconstrained=FALSE)),AIC=AIC(m),
    delta_AIC_vs_baseline=AIC(m)-AIC(base$model),error='')
  }
 }
}
gam_results<-do.call(rbind,gam_results);gls_results<-do.call(rbind,gls_results)
write.csv(gam_results,file.path(out,'common_cohort_GAM_models.csv'),row.names=FALSE)
write.csv(gls_results,file.path(out,'common_cohort_CAR1_models.csv'),row.names=FALSE)
cat('COMMON 0..4 n=',length(ids04),' 0..8 n=',length(ids08),'\n')
print(table(d$station[d$sample%in%ids04]))
print(assoc[assoc$cohort=='lag_specific'&assoc$lag_weeks%in%c(0,2,3,4)&assoc$response=='log_ratio'&assoc$adjustment=='unadjusted',],row.names=FALSE)
print(assoc[assoc$cohort=='common_0_to_4'&assoc$response=='log_ratio'&assoc$station=='pooled_station_adjusted',],row.names=FALSE)
print(gam_results[gam_results$response=='log_ratio_eps1e6',],row.names=FALSE)
print(gls_results,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
cat('DENSE_TEMPERATURE_ANALYSIS_COMPLETE\n')
