# P8108 HW2: supplied computing tools for the two-task assignment.
# Use unchanged or implement the same definitions in other software.
# Source this file; no analysis is run automatically.
weighted_lr <- function(dat, rho=0, gamma=0) {
 stopifnot(all(c("time","status","group") %in% names(dat)), all(is.finite(dat$time)), all(dat$time>0), all(dat$status %in% 0:1), all(dat$group %in% 0:1), rho>=0, gamma>=0)
 tt<-sort(unique(dat$time)); ix<-match(dat$time,tt); k<-length(tt)
 exit0<-tabulate(ix[dat$group==0],nbins=k);exit1<-tabulate(ix[dat$group==1],nbins=k)
 n0<-rev(cumsum(rev(exit0)));n1<-rev(cumsum(rev(exit1)))
 d0<-tabulate(ix[dat$group==0 & dat$status==1],nbins=k);d1<-tabulate(ix[dat$group==1 & dat$status==1],nbins=k)
 take<-(d0+d1)>0;tt<-tt[take];n0<-n0[take];n1<-n1[take];d0<-d0[take];d1<-d1[take]
 n<-n0+n1;d<-d0+d1;after<-cumprod(1-d/n);before<-if(length(after)) c(1,head(after,-1)) else numeric()
 weight<-before^rho*(1-before)^gamma;e1<-n1*d/n;v<-numeric(length(n));ok<-n>1;v[ok]<-n0[ok]*n1[ok]*d[ok]*(n[ok]-d[ok])/(n[ok]^2*(n[ok]-1))
 tab<-data.frame(time=tt,n0,n1,d0,d1,Sminus=before,weight,e1,u=d1-e1,v)
 U<-sum(weight*(d1-e1));V<-sum(weight^2*v)
 list(U=U,V=V,chisq=if(V>0) U^2/V else NA_real_,p=if(V>0)pchisq(U^2/V,1,lower.tail=FALSE) else NA_real_,table=tab)
}

# The exponential pieces are a simulation device, not a model to fit or derive.
# Group 0 hazard = 0.06/month. Group 1 = 0.06 through delay, then 0.06*hr_after.
simulate_trial <- function(delay, hr_after, censor_rate, tau, seed, n_per_group=80L) {
  stopifnot(is.finite(delay), is.finite(tau), delay>0, delay<tau,
            is.finite(hr_after), hr_after>0, hr_after<=1,
            is.finite(censor_rate), censor_rate>=0, n_per_group>=2)
  set.seed(seed)
  group <- rep(0:1, each=n_per_group)
  e <- rexp(length(group), rate=1)
  after <- ifelse(group==0, .06, .06*hr_after)
  T <- ifelse(e<=.06*delay, e/.06, delay+(e-.06*delay)/after)
  C <- if(censor_rate==0) rep(Inf,length(group)) else rexp(length(group),censor_rate)
  data.frame(time=pmin(T,C,tau), status=as.integer(T<=C & T<=tau), group=group)
}

validate_design <- function(plan) {
  need <- c("scenario","delay","hr_after","censor_rate","tau")
  stopifnot(all(need %in% names(plan)), nrow(plan)==4L,
            !anyDuplicated(plan$scenario), all(nzchar(plan$scenario)))
  stopifnot(all(is.finite(as.matrix(plan[,need[-1]]))),
            all(plan$delay>0 & plan$delay<plan$tau),
            all(plan$hr_after>0 & plan$hr_after<1), all(plan$censor_rate>=0))
  counts <- vapply(plan[,need[-1]], function(z)length(unique(z)), integer(1))
  stopifnot(sum(counts==2L)==2L, sum(counts==1L)==2L,
            nrow(unique(plan[,need[-1]]))==4L)
  invisible(TRUE)
}

make_plan <- function(alternatives) {
  validate_design(alternatives)
  a <- alternatives[,c("scenario","delay","hr_after","censor_rate","tau")]
  a$kind <- "alternative"
  # Under hr_after=1 the delay has no effect on the event-time distribution.
  # Reuse the null for alternatives sharing censoring and follow-up.
  null <- a[!duplicated(a[,c("censor_rate","tau")]),,drop=FALSE]
  null$scenario <- paste0("null_",seq_len(nrow(null)))
  null$hr_after <- 1; null$kind <- "null"
  stopifnot(!any(null$scenario %in% a$scenario))
  key <- function(x)paste(x$censor_rate,x$tau,sep="|")
  a$matched_null <- null$scenario[match(key(a),key(null))]
  null$matched_null <- null$scenario
  rbind(a,null)
}

run_plan <- function(alternatives, B=1000L, seed=81082026L, output_dir="results/simulation") {
  stopifnot(B>=2, B==as.integer(B), B<=10000)
  plan <- make_plan(alternatives)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  set.seed(seed)
  seeds <- matrix(sample.int(.Machine$integer.max,nrow(plan)*B,replace=FALSE),nrow=nrow(plan))
  out <- vector("list",nrow(plan)*B*2); pos <- 0L
  for(h in seq_len(nrow(plan))) {
    message("Running ",plan$scenario[h]," (",B," datasets)")
    z <- plan[h,]
    for(b in seq_len(B)) {
      d <- simulate_trial(z$delay,z$hr_after,z$censor_rate,z$tau,seeds[h,b])
      for(m in c("LR","FH01")) {
        error <- ""
        ans <- tryCatch(weighted_lr(d,0,as.integer(m=="FH01")),
                        error=function(e){error<<-conditionMessage(e); list(U=NA,V=NA,chisq=NA,p=NA)})
        pos<-pos+1L
        out[[pos]] <- data.frame(scenario=z$scenario,kind=z$kind,matched_null=z$matched_null,
          replicate=b,seed=seeds[h,b],method=m,U=ans$U,V=ans$V,chisq=ans$chisq,p=ans$p,
          events=sum(d$status),censor_fraction=mean(d$status==0),
          at_risk_at_delay=sum(d$time>=z$delay),
          events_before_delay=sum(d$status==1 & d$time<=z$delay),
          events_after_delay=sum(d$status==1 & d$time>z$delay),error=error)
      }
    }
  }
  out <- do.call(rbind,out)
  summ <- do.call(rbind,lapply(split(out,list(out$scenario,out$method),drop=TRUE),function(x){
    ok<-is.finite(x$p); N<-sum(ok); rate<-if(N)mean(x$p[ok]<.05) else NA_real_
    data.frame(scenario=x$scenario[1],kind=x$kind[1],matched_null=x$matched_null[1],
      method=x$method[1],B=nrow(x),valid=N,failures=sum(!ok),reject_rate=rate,
      MCSE=if(N)sqrt(rate*(1-rate)/N) else NA_real_,
      mean_events=mean(x$events),mean_censor_fraction=mean(x$censor_fraction),
      mean_at_risk_at_delay=mean(x$at_risk_at_delay),
      mean_events_before_delay=mean(x$events_before_delay),mean_events_after_delay=mean(x$events_after_delay))
  }))
  paired <- do.call(rbind,lapply(split(out,out$scenario),function(x){
    a<-x[x$method=="FH01",];b<-x[x$method=="LR",];b<-b[match(a$replicate,b$replicate),]
    ok<-is.finite(a$p)&is.finite(b$p);D<-as.integer(a$p[ok]<.05)-as.integer(b$p[ok]<.05)
    delta<-if(length(D))mean(D) else NA_real_;se<-if(length(D)>1)sd(D)/sqrt(length(D)) else NA_real_
    data.frame(scenario=x$scenario[1],kind=x$kind[1],valid_pairs=length(D),
      difference=delta,MCSE=se,lower=delta-1.96*se,upper=delta+1.96*se)
  }))
  rownames(summ)<-rownames(paired)<-NULL
  write.csv(plan,file.path(output_dir,"design.csv"),row.names=FALSE)
  write.csv(out,file.path(output_dir,"replicates.csv"),row.names=FALSE,na="NA")
  write.csv(summ,file.path(output_dir,"summary.csv"),row.names=FALSE,na="NA")
  write.csv(paired,file.path(output_dir,"paired_differences.csv"),row.names=FALSE,na="NA")
  writeLines(c(paste("B per cell:",B),paste("Master seed:",seed),
               "Two-sided alpha=0.05; n=80 per group; lambda0=0.06/month",capture.output(sessionInfo())),
             file.path(output_dir,"sessionInfo.txt"))
  invisible(list(design=plan,replicates=out,summary=summ,paired=paired))
}
