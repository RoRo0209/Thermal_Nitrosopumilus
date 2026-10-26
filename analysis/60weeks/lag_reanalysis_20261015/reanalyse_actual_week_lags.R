#!/usr/bin/env Rscript
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

# Set THERMAL_PROJECT_ROOT; pass --figures to create the optional figures.
args <- commandArgs(trailingOnly=TRUE)
make_figures <- '--figures' %in% args
root_args <- args[args!='--figures']
root <- .release_root
out <- file.path(root,'60weeks','lag_reanalysis_20261015')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Use the active R library environment.
required <- if(make_figures)c('ggplot2','patchwork','svglite','ragg') else character()
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Required R packages unavailable: ',paste(missing,collapse=', '))
if(make_figures) suppressPackageStartupMessages({library(ggplot2);library(patchwork)})

book <- .release_input("THERMAL_TIME_SERIES_SOURCE_WORKBOOK")
# Frozen workbook extraction, independently reconciled to the workbook using openpyxl.
# Stop rather than silently use a stale extraction if the source workbook changes.
stopifnot(identical(unname(tools::md5sum(book)),
                   readLines(file.path(out,'input_workbook_md5.txt'),warn=FALSE)))
d <- read.csv(file.path(out,'input_submitted_Data3.csv'),stringsAsFactors=FALSE)
d <- d[order(d$station,d$week),]; rownames(d)<-NULL
old_path <- file.path(root,'60weeks','out_60weeks_by_station_mgonly','60weeks_station_week_table_mgonly.csv')
old <- read.csv(old_path); old <- old[match(d$sample,old$sample),]
stopifnot(nrow(d)==39,all(table(d$station)==13),!anyDuplicated(d$sample),
          identical(sort(d$sample),sort(old$sample)),all(is.finite(as.matrix(d[,c('week','temperature','HighT','LowT')]))),
          all(d$HighT>0),all(d$LowT>0),all(d$week==old$week_num),
          max(abs(d$temperature-old$Temperature))<1e-12,
          max(abs(d$HighT-old$HighT))<1e-12,max(abs(d$LowT-old$LowT))<1e-12,
          max(abs(d$submitted_ratio-d$HighT/d$LowT))<1e-12)
d$log_ratio <- log10(d$HighT/d$LowT)
d$log_ratio_eps1e6 <- log10((d$HighT+1e-6)/(d$LowT+1e-6))
d$log_ratio_station_eps <- old$ratio
write.csv(d,file.path(out,'source_observations.csv'),row.names=FALSE)
write.csv(data.frame(source=c('Submitted Supplementary Data 3','Original analysis table'),
 md5=unname(tools::md5sum(c(book,old_path)))),file.path(out,'source_checksums.csv'),row.names=FALSE)
stations<-c('E','W','Z'); lags <- -16:16
common_weeks <- Reduce(intersect,lapply(stations,function(s)d$week[d$station==s]))
write.csv(data.frame(week=sort(common_weeks)),file.path(out,'common_sampling_weeks.csv'),row.names=FALSE)

residualize <- function(g, response, adjustment) {
  z<-data.frame(z=g[[response]],week=g$week,
               annual_sin=sin(2*pi*g$week/52),annual_cos=cos(2*pi*g$week/52))
  if(adjustment=='unadjusted') return(z$z)
  formula<-if(adjustment=='annual_harmonic') z~annual_sin+annual_cos else z~annual_sin+annual_cos+week
  model<-lm(formula,data=z)
  stopifnot(model$rank==length(coef(model)))
  as.numeric(residuals(model))
}
pair_indices <- function(g,k) {
  j<-match(g$week+k,g$week); i<-which(!is.na(j)); j<-j[i]
  stopifnot(all(g$week[j]-g$week[i]==k),!anyDuplicated(i),!anyDuplicated(j))
  list(i=i,j=j)
}
estimate <- function(data,response,adjustment,dataset='all_observed',min_n=3) {
  ans<-list()
  for(s in stations) {
    g<-data[data$station==s,];g<-g[order(g$week),]
    x<-residualize(g,'temperature',adjustment); y<-residualize(g,response,adjustment)
    for(k in lags) {
      ij<-pair_indices(g,k); i<-ij$i;j<-ij$j;n<-length(i)
      ok<-n>=min_n && sd(x[i])>1e-12 && sd(y[j])>1e-12
      ans[[length(ans)+1]]<-data.frame(dataset=dataset,station=s,response=response,
       adjustment=adjustment,lag_weeks=k,n_pairs=n,
       pearson_r=if(ok)cor(x[i],y[j]) else NA_real_,
       spearman_rho=if(ok)cor(x[i],y[j],method='spearman') else NA_real_)
    }
  }
  do.call(rbind,ans)
}
responses<-c('log_ratio','log_ratio_eps1e6','log_ratio_station_eps')
adjustments<-c('unadjusted','annual_harmonic','annual_harmonic_linear_trend')
all_estimates<-list()
for(ds in c('all_observed','common_weeks')) {
 sub<-if(ds=='all_observed')d else d[d$week %in% common_weeks,]
 for(resp in responses) for(adj in adjustments)
   all_estimates[[length(all_estimates)+1]]<-estimate(sub,resp,adj,ds)
}
est<-do.call(rbind,all_estimates)
write.csv(est,file.path(out,'exact_week_lag_estimates.csv'),row.names=FALSE)

# Each source pairing is retained; weeks, not row positions, define the lag.
pairs<-list();support<-list()
for(s in stations) {
 g<-d[d$station==s,]
 for(k in lags) {
  ij<-pair_indices(g,k);i<-ij$i;j<-ij$j
  support[[length(support)+1]]<-data.frame(station=s,lag_weeks=k,n_pairs=length(i),
     n_observed_weeks=length(unique(c(g$week[i],g$week[j]))))
  if(length(i)) pairs[[length(pairs)+1]]<-data.frame(station=s,lag_weeks=k,
    temperature_sample=g$sample[i],temperature_week=g$week[i],temperature=g$temperature[i],
    ratio_sample=g$sample[j],ratio_week=g$week[j],log_ratio=g$log_ratio[j])
 }
}
support<-do.call(rbind,support);pairs<-do.call(rbind,pairs)
stopifnot(all(support$n_pairs[support$lag_weeks==2]==0),
          all(support$n_pairs[support$lag_weeks==4]==7))
write.csv(support,file.path(out,'exact_week_pair_counts.csv'),row.names=FALSE)
write.csv(pairs,file.path(out,'explicit_actual_week_pairs.csv'),row.names=FALSE)

# Influence analysis is not a bootstrap or a confidence-interval calculation.
deletions<-list()
for(w in sort(unique(d$week))) {
 sub<-d[d$week!=w,]
 for(adj in adjustments) {
  z<-estimate(sub,'log_ratio',adj,min_n=4)
  z$omitted_week<-w
  deletions[[length(deletions)+1]]<-z
 }
}
del<-do.call(rbind,deletions)
write.csv(del,file.path(out,'leave_one_week_out_estimates.csv'),row.names=FALSE)
influence<-do.call(rbind,lapply(split(del,list(del$station,del$adjustment,del$lag_weeks),drop=TRUE),function(z){
 valid<-is.finite(z$pearson_r)
 data.frame(station=z$station[1],adjustment=z$adjustment[1],lag_weeks=z$lag_weeks[1],
  influence_min=if(any(valid))min(z$pearson_r[valid]) else NA_real_,
  influence_max=if(any(valid))max(z$pearson_r[valid]) else NA_real_,
  n_eligible_deletions=sum(valid),min_remaining_pairs=if(any(valid))min(z$n_pairs[valid]) else NA_integer_)
}))
write.csv(influence,file.path(out,'leave_one_week_out_ranges.csv'),row.names=FALSE)

# Export observed/residual data and model summaries for inspectable season adjustment.
season<-list();fits<-list()
for(s in stations) {
 g<-d[d$station==s,]
 for(adj in adjustments) {
  season[[length(season)+1]]<-data.frame(station=s,sample=g$sample,week=g$week,adjustment=adj,
    temperature_value=residualize(g,'temperature',adj),ratio_value=residualize(g,'log_ratio',adj))
 }
 for(v in c('temperature','log_ratio')) {
  z<-data.frame(value=g[[v]],annual_sin=sin(2*pi*g$week/52),annual_cos=cos(2*pi*g$week/52),week=g$week)
  for(adj in adjustments[-1]) {
   fm<-if(adj=='annual_harmonic') value~annual_sin+annual_cos else value~annual_sin+annual_cos+week
   m<-lm(fm,data=z)
   fits[[length(fits)+1]]<-data.frame(station=s,variable=v,adjustment=adj,
     n=nrow(g),residual_df=df.residual(m),r_squared=summary(m)$r.squared,
     residual_sd=sqrt(sum(resid(m)^2)/df.residual(m)))
  }
 }
}
write.csv(do.call(rbind,season),file.path(out,'seasonal_residuals.csv'),row.names=FALSE)
write.csv(do.call(rbind,fits),file.path(out,'seasonal_fit_diagnostics.csv'),row.names=FALSE)
primary<-est[est$dataset=='all_observed' & est$response=='log_ratio',]
primary<-merge(primary,influence,by=c('station','adjustment','lag_weeks'),all.x=TRUE)
write.csv(primary[primary$lag_weeks %in% 2:4,],file.path(out,'focal_2_to_4_week_results.csv'),row.names=FALSE)
write.csv(primary[primary$n_pairs>=6,],file.path(out,'displayed_estimates_and_influence.csv'),row.names=FALSE)

# R-only figure exports: independent points, no interpolation across unsupported lags.
if(make_figures) {
pal<-c(E='#0072B2',W='#D55E00',Z='#009E73');shapes<-c(E=16,W=17,Z=15)
theme_set(theme_classic(base_size=8,base_family='Arial')+
 theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7,colour='#444444'),
       axis.title=element_text(size=8),axis.text=element_text(size=7),
       legend.position='bottom',legend.title=element_text(size=8),legend.text=element_text(size=7),
       axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.3),
       plot.margin=margin(4,7,4,4)))
d$station<-factor(d$station,levels=rev(stations))
p_a<-ggplot(d,aes(week,station,colour=station,shape=station))+
 geom_segment(aes(x=4,xend=60,yend=station),colour='#DDDDDD',linewidth=.4)+
 geom_point(size=2.1)+scale_colour_manual(values=pal)+scale_shape_manual(values=shapes)+
 scale_x_continuous(breaks=c(seq(0,56,8),60),limits=c(0,62))+
 labs(title='a  Observed sampling calendar',subtitle='13 metagenomes per station; no weekly interpolation',
      x='Sampling week',y=NULL)+theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank())
sp<-support[support$lag_weeks>=0,];sp$station<-factor(sp$station,levels=rev(stations))
p_b<-ggplot(sp,aes(lag_weeks,station,fill=n_pairs))+
 geom_tile(width=.92,height=.85)+geom_text(aes(label=n_pairs,colour=n_pairs>=9),size=2.35)+
 scale_colour_manual(values=c('FALSE'='#222222','TRUE'='white'),guide='none')+
 scale_fill_gradient(low='#F5F5F5',high='#366B89',limits=c(0,13),guide='none')+
 scale_x_continuous(breaks=0:16,expand=expansion(add=.5))+
 labs(title='b  Number of observed pairs at each exact time difference',
      subtitle='Counts are per station and are identical for positive and negative lags',
      x='Absolute time difference (weeks)',y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank())
plot_cor<-function(adj,title) {
 z<-primary[primary$adjustment==adj & primary$n_pairs>=6,]
 z$station<-factor(z$station,levels=stations)
 # Small within-lag visual offsets separate station glyphs; data retain exact integer lags.
 z$x_display<-z$lag_weeks+c(E=-.45,W=0,Z=.45)[as.character(z$station)]
 ggplot(z,aes(x_display,pearson_r,colour=station,shape=station))+
  geom_hline(yintercept=0,linetype='dashed',linewidth=.35,colour='#777777')+
  geom_vline(xintercept=0,linewidth=.3,colour='#DDDDDD')+
  geom_linerange(aes(ymin=influence_min,ymax=influence_max),linewidth=.45,alpha=.5,show.legend=FALSE)+
  geom_point(size=2)+scale_colour_manual(values=pal,name='Station')+
  scale_shape_manual(values=shapes,name='Station')+
  scale_x_continuous(breaks=seq(-16,16,4),limits=c(-17.5,17.5))+
  scale_y_continuous(breaks=c(-1,-.5,0,.5,1),limits=c(-1.02,1.02))+
  labs(title=title,subtitle='Shown only when at least 6 observed pairs are available',
       x='Lag (weeks; positive = temperature first)',y='Pearson correlation')
}
p_c<-plot_cor('unadjusted','c  Actual-time association')
p_d<-plot_cor('annual_harmonic','d  Annual-harmonic residual association')
fig<-p_a / p_b / (p_c | p_d)+plot_layout(heights=c(.65,.75,1.7),guides='collect')+
 plot_annotation(caption=paste('Positive lag pairs temperature at week t with the ratio at week t + lag.',
   'Vertical ranges show leave-one-week-out sensitivity, not confidence intervals.',
   'Separate station estimates; no significance thresholds or optimal-lag inference.',sep='\n'),
   theme=theme(plot.caption=element_text(size=7,hjust=0),plot.margin=margin(5,5,5,5)))
w<-183/25.4;h<-165/25.4;stem<-file.path(out,'Figure_S8_actual_week_lag_reanalysis')
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='Arial');print(fig);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=w,height=h);print(fig);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=300);print(fig);dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(fig);dev.off()
}

cat('\nFOCAL RESULTS (primary log-ratio)\n')
print(primary[primary$lag_weeks==4,c('station','adjustment','n_pairs','pearson_r','spearman_rho','influence_min','influence_max')],row.names=FALSE)
cat('\nPRIMARY DISPLAYED UNADJUSTED RESULTS\n')
print(primary[primary$adjustment=='unadjusted' & primary$n_pairs>=6,c('station','lag_weeks','n_pairs','pearson_r')],row.names=FALSE)
cat('\nCOMMON WEEKS\n');print(sort(common_weeks))
capture.output(sessionInfo(),file=file.path(out,if(make_figures)'R_session_info.txt' else 'R_session_info_analysis_only.txt'))
writeLines(c('All source samples reconciled; all abundances positive.',
 'Every retained pairing satisfies ratio_week - temperature_week = lag_weeks.',
 'At +2 weeks: zero pairs per station. At +4 weeks: seven pairs per station.',
 'Correlations with fewer than three pairs are missing, not zero.',
 'Plots display n_pairs >= 6. Leave-one-week ranges require >=4 remaining pairs.',
 'No station-pooled sample size, nominal P values or white-noise significance bounds.',
 'All graphics generated in R. Text in SVG/PDF remains text.',
 'PNG must be visually inspected before delivery.'),file.path(out,'implementation_checks.txt'))
cat(if(make_figures)'\nANALYSIS_AND_FIGURE_EXPORT_COMPLETE\n' else '\nANALYSIS_TABLES_COMPLETE (append --figures to export plots)\n')
