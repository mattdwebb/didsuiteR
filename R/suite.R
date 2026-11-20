.timed <- function(expr) {
  start <- proc.time()[[3L]]
  out <- .attempt(expr)
  out$elapsed <- unname(proc.time()[[3L]]-start)
  out
}
.native <- function(adapter,p,options,null) {
  f <- adapter(p,options,TRUE)
  .require(length(f$native_se)==1L && is.finite(f$native_se) && f$native_se>0,
    "Native variance unavailable or zero","IDENTIFICATION_FAILURE")
  df <- f$native_df
  .require(length(df)==1L && !is.na(df) && df>0,"Native reference degrees of freedom unavailable","UNSUPPORTED")
  stat <- (f$estimate-null)/f$native_se
  list(se=f$native_se,stat=stat,df=df,p=2*stats::pt(-abs(stat),df),
       method=f$native_method,raw=f)
}
.jack <- function(adapter,p,options,observed,null) {
  ids <- .ids(p$data[[p$cluster]]); G <- length(ids)
  .require(G>=2L,"Jackknife requires at least two inference clusters","IDENTIFICATION_FAILURE")
  values <- setNames(rep(NA_real_,G),ids)
  policy <- if (is.null(observed$jackknife_support)) "changing" else observed$jackknife_support
  .require(policy %in% c("fixed","changing"),"Invalid jackknife support policy","UNSUPPORTED")
  support_changed <- setNames(rep(FALSE,G),ids)
  failures <- setNames(character(G), ids)
  for (i in seq_along(ids)) {
    pp <- p; pp$data <- p$data[as.character(p$data[[p$cluster]])!=ids[i],,drop=FALSE]
    pp$cluster_deletion <- TRUE
    cohort_changed <- !setequal(unique(pp$data$.didsuite_g),unique(p$data$.didsuite_g))
    refit <- .attempt({
      if (policy=="fixed") .require(!cohort_changed,
        paste("Deletion removes required cohort:",ids[i]),"IDENTIFICATION_FAILURE")
      f <- adapter(pp,options,FALSE)
      support_changed[i] <- cohort_changed ||
        (!is.null(observed$support) && !identical(f$support,observed$support))
      if (policy=="fixed" && !is.null(observed$support)) .require(identical(f$support,observed$support),
        paste("Deletion changes estimand support:",ids[i]),"IDENTIFICATION_FAILURE")
      .require(length(f$estimate)==1L && is.finite(f$estimate),"Undefined deletion","IDENTIFICATION_FAILURE")
      f$estimate
    })
    if (refit$status=="PASS") values[i] <- refit$value else failures[i] <- paste(ids[i], refit$reason)
  }
  if (any(!is.finite(values))) return(list(status="IDENTIFICATION_FAILURE", failures=failures,
    deletions=values, df=G-1L, support_policy=policy, support_changed=support_changed, stat=NA_real_))
  se <- sqrt((G-1)/G*sum((values-mean(values))^2))
  tolerance <- 64*.Machine$double.eps*max(1,abs(observed$estimate),abs(values))
  .require(is.finite(se) && se>tolerance,"Zero jackknife variance","IDENTIFICATION_FAILURE")
  stat <- (observed$estimate-null)/se
  list(se=se,stat=stat,df=G-1,p=2*stats::pt(-abs(stat),G-1),
       deletions=values,centering="mean deletion estimate",
       support_policy=policy,support_changed=support_changed)
}
.ri <- function(method,adapter,p,options,draws,observed) {
  statistic <- function(pp) {
    f <- adapter(pp,options,FALSE)
    if (identical(observed$jackknife_support,"fixed") && !is.null(observed$support)) .require(identical(f$support,observed$support),
      "Assignment changes estimand support","IDENTIFICATION_FAILURE")
    if (method=="ri_beta") return(f$estimate)
    if (method=="ri_t_native") return(.native(adapter,pp,options,0)$stat)
    stop("Unsupported randomization statistic",call.=FALSE)
  }
  obs <- statistic(p)
  values <- rep(NA_real_,nrow(draws$matrix)); failure <- rep("",length(values))
  for (i in seq_along(values)) {
    ans <- .attempt(statistic(.reassign(p,draws$matrix[i,])))
    if (ans$status=="PASS") values[i] <- ans$value else failure[i] <- ans$reason
  }
  ok <- is.finite(values)
  hits <- ifelse(ok,abs(values)>=abs(obs)-64*.Machine$double.eps*pmax(1,abs(values),abs(obs)),NA)
  list(p=if (all(ok)) (1+sum(hits))/(length(values)+1) else NA_real_,
    statistic=obs,statistics=values,hit_indicators=hits,hits=sum(hits,na.rm=TRUE),
    draws=length(values),valid_draws=sum(ok),invalid_draws=sum(!ok),failures=failure,
    status=if (all(ok)) "PASS" else "IDENTIFICATION_FAILURE")
}
.blank_row <- function(id,p) {
  nums <- c("estimate","native_se","native_stat","native_df","native_p",
    "jackknife_se","jackknife_t","jackknife_df","jackknife_p","ri_beta_p",
    "ri_t_native_p")
  out <- as.list(setNames(rep(NA_real_,length(nums)),nums))
  c(list(estimator=id,language="R",estimator_version=NA_character_,
    target_label=NA_character_,target_mode=p$target_mode,target_compatible=FALSE,
    native_method=NA_character_),out,
    list(ri_reps=0L,assignment_space=NA_real_,assignment_hash=NA_character_,
      nobs=nrow(p$data),nclusters=length(unique(p$data[[p$cluster]])),
      n_treated_clusters=length(unique(p$data[[p$cluster]][p$data$.didsuite_g>0])),
      cohorts_used=paste(sort(unique(p$cohorts[p$cohorts>0])),collapse=","),
      control_group=p$control_group,data_mode=p$mode,assignment_unit=p$assignment_unit,
      inference_cluster=p$cluster,covariates_requested=paste(p$covariates,collapse=","),
      covariates_used=NA_character_,covariate_mode=NA_character_,weighting_scheme=NA_character_,
      sample_mode=p$sample_mode,event_window=paste(p$event_window,collapse=","),
      cells_used=NA_integer_, cells_eligible=NA_integer_, min_controls=NA_integer_, max_controls=NA_integer_,
      jackknife_support=NA_character_,
      runtime_setup=0,runtime_point=0,runtime_native=0,runtime_jk=0,runtime_ri=0,
      runtime_total=0,status="UNSUPPORTED",failure_reason=""))
}
did_suite <- function(data,outcome,time,cluster,estimators="all",gvar=NULL,treatment=NULL,
                      assignment_unit=cluster,unit=NULL,covariates=character(),weights=NULL,
                      data_mode="auto",inference="standard",assignments=NULL,reps=999L,seed=0L,
                      exact=FALSE,null=0,target_mode="native",sample_mode="native",
                      event_window=NULL,control_group="native",estimator_options=list(),
                      adapters=list()) {
  start <- proc.time()[[3L]]
  target_mode <- match.arg(target_mode,c("native","harmonized"))
  sample_mode <- match.arg(sample_mode,c("native","common"))
  .require(length(null)==1 && is.finite(null),"Finite scalar null required")
  .require(is.list(estimator_options) && is.list(adapters),"Options and adapters must be lists")
  for (config in list(estimator_options,adapters)) .require(!length(config) ||
    (!is.null(names(config)) && all(nzchar(names(config))) && !anyDuplicated(names(config))),
    "Estimator options and adapters must have unique names")
  methods <- .methods(inference)
  catalog <- discover_estimators()
  if (identical(estimators,"all")) estimators <- catalog$estimator
  .require(length(estimators)>0 && !anyDuplicated(estimators) &&
             all(estimators %in% catalog$estimator),"Unknown or duplicate estimator selection")
  .require(all(names(estimator_options) %in% estimators) && all(names(adapters) %in% estimators),
    "Unused estimator options or adapters")
  p <- .prepare(data,outcome,time,cluster,assignment_unit,unit,gvar,treatment,
                covariates,weights,data_mode)
  p$target_mode <- target_mode; p$sample_mode <- sample_mode
  p$event_window <- event_window; p$control_group <- control_group
  setup <- unname(proc.time()[[3L]]-start)
  ri_methods <- intersect(methods,c("ri_beta","ri_t_native"))
  draws <- if (length(ri_methods)) .attempt(assignment_set(p$cohorts,p$ids,reps,seed,exact,assignments)) else NULL
  if (!length(ri_methods) && !is.null(assignments)) .abort("Assignments supplied without an RI method")
  rows <- list(); details <- list(); event_rows <- list()
  prepared <- p
  for (id in estimators) {
    p <- prepared
    options <- estimator_options[[id]]; if (is.null(options)) options <- list()
    if ("event_window" %in% names(options)) {
      p$event_window <- options$event_window; options$event_window <- NULL
    }
    if ("control_group" %in% names(options)) {
      p$control_group <- options$control_group; options$control_group <- NULL
    }
    clock <- proc.time()[[3L]]
    row <- .blank_row(id,p); row$runtime_setup <- setup
    row$estimator_version <- catalog$version[match(id,catalog$estimator)]
    adapter <- if (id %in% names(adapters)) adapters[[id]] else .adapter(id)
    .require(is.function(adapter),"Adapter must be a function")
    point <- .timed(adapter(p,options,FALSE))
    row$runtime_point <- point$elapsed
    detail <- list(point=point,methods=list())
    if (point$status=="PASS") {
      f <- point$value
      if (!is.null(f$att_st)) {
        row$cells_used <- nrow(f$att_st)
        row$min_controls <- min(f$att_st$control_count)
        row$max_controls <- max(f$att_st$control_count)
      }
      if (id=="csdid" && !is.null(f$raw$att_gt)) {
        row$cells_used <- sum(is.finite(f$raw$att_gt))
        row$cells_eligible <- length(f$raw$att_gt)
      }
      if (!is.null(f$cells_eligible)) row$cells_eligible <- f$cells_eligible
      row$estimate <- f$estimate; row$target_label <- f$target_label
      row$jackknife_support <- if (is.null(f$jackknife_support)) "changing" else f$jackknife_support
      if (!is.null(f$control_group)) row$control_group <- f$control_group
      if (!is.null(f$cohorts)) row$cohorts_used <- paste(f$cohorts,collapse=",")
      if (!is.null(f$n_used)) row$nobs <- f$n_used
      if (!is.null(f$event_window_used)) row$event_window <- paste(f$event_window_used,collapse=",")
      row$covariates_used <- paste(p$covariates,collapse=",")
      row$covariate_mode <- f$covariate_mode; row$weighting_scheme <- f$weighting_scheme
      if (is.null(row$covariate_mode)) row$covariate_mode <- "adapter-defined"
      if (is.null(row$weighting_scheme)) row$weighting_scheme <- "adapter-defined"
      if (!is.null(f$event_study) && nrow(f$event_study)>0) {
        event_rows[[id]] <- cbind(estimator=id,f$event_study)
      }
      for (method in methods) {
        result <- .timed({
          if (method=="native") .native(adapter,p,options,null)
          else if (method=="jackknife") .jack(adapter,p,options,f,null)
          else {
            .require(null==0,"RI currently supports zero sharp null only","UNSUPPORTED")
            .require(draws$status=="PASS",draws$reason,draws$status)
            .ri(method,adapter,p,options,draws$value,f)
          }
        })
        if (result$status=="PASS" && !is.null(result$value$status) && result$value$status!="PASS") {
          result$status <- result$value$status
          result$reason <- paste(unique(result$value$failures[nzchar(result$value$failures)]),collapse="; ")
        }
        detail$methods[[method]] <- result
        if (method=="native") row$runtime_native <- result$elapsed
        else if (method=="jackknife") row$runtime_jk <- result$elapsed
        else row$runtime_ri <- row$runtime_ri+result$elapsed
        if (result$status!="PASS") next
        v <- result$value
        if (method=="native") {
          row$native_se <- v$se; row$native_stat <- v$stat; row$native_df <- v$df
          row$native_p <- v$p; row$native_method <- v$method
        } else if (method=="jackknife") {
          row$jackknife_support <- v$support_policy
          row$jackknife_se <- v$se; row$jackknife_t <- v$stat
          row$jackknife_df <- v$df; row$jackknife_p <- v$p
        } else row[[paste0(method,"_p")]] <- v$p
      }
      bad <- vapply(detail$methods,function(x) x$status!="PASS",logical(1))
      row$status <- if (any(bad)) "PARTIAL" else "PASS"
      row$failure_reason <- paste(vapply(names(detail$methods)[bad],function(m)
        paste0(m,": ",detail$methods[[m]]$reason),character(1)),collapse="; ")
    } else {
      row$status <- point$status; row$failure_reason <- point$reason
    }
    if (!is.null(draws) && draws$status=="PASS") {
      row$ri_reps <- nrow(draws$value$matrix); row$assignment_space <- draws$value$space
      row$assignment_hash <- draws$value$hash
    }
    row$runtime_total <- unname(proc.time()[[3L]]-clock)+setup
    rows[[id]] <- as.data.frame(row,stringsAsFactors=FALSE)
    details[[id]] <- detail
  }
  structure(list(summary=do.call(rbind,rows),details=details,event_study=event_rows,
    assignments=if (!is.null(draws) && draws$status=="PASS") draws$value else NULL,
    assignment_status=draws,input=list(assignment_ids=p$ids,cohorts=p$cohorts,methods=methods),
    certification="LOCAL_HARNESS_TESTS_ONLY_NOT_NINE_ESTIMATOR_EXTERNAL_CERTIFICATION"),
    class="didsuite_result")
}
print.didsuite_result <- function(x,...) {
  print(x$summary[,c("estimator","estimate","native_p","jackknife_p","ri_beta_p",
                      "ri_t_native_p","status","target_label")],row.names=FALSE,...)
  invisible(x)
}
export_results <- function(x,path,format=c("csv","json")) {
  .require(inherits(x,"didsuite_result"),"Expected didsuite_result")
  format <- match.arg(format)
  if (format=="csv") utils::write.csv(x$summary,path,row.names=FALSE,na="")
  else jsonlite::write_json(list(summary=x$summary,event_study=x$event_study),path,
                            dataframe="rows",na="null",auto_unbox=TRUE,pretty=TRUE)
  invisible(path)
}

