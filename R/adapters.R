.catalog <- function() data.frame(
  estimator=c("twfe","did2s","bjs","flex","csdid","undid","didint","stacked","lpdid"),
  package=c("fasttwfeR","fastdid2sR","fastbjsR","fastflexR","fastrcsdid","suiteundidR","suitedidintR","faststackedR","fastlpdidR"),
  entrypoint=c("fasttwfe","fastdid2s","fastbjs","fastflex","fast_csdid","suite_undid","suite_didint","fit","fit"),
  native=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE),
  stringsAsFactors=FALSE)
discover_estimators <- function() {
  x <- .catalog()
  x$available <- mapply(function(pkg,fun) requireNamespace(pkg,quietly=TRUE) &&
    fun %in% getNamespaceExports(pkg),x$package,x$entrypoint)
  x$version <- vapply(x$package,function(pkg)
    if (requireNamespace(pkg,quietly=TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_,
    character(1))
  x$feature_status <- ifelse(x$available,"ENTRYPOINT_AVAILABLE_NOT_CERTIFICATION","NOT_INSTALLED_OR_ENTRYPOINT_MISSING")
  x
}
.check_options <- function(options, allowed) {
  .require(is.list(options) && (length(options)==0 ||
    (!is.null(names(options)) && !anyDuplicated(names(options)) && all(nzchar(names(options))))),
    "Estimator options must be a uniquely named list")
  .require(all(names(options) %in% allowed), paste("Unsupported options:",paste(setdiff(names(options),allowed),collapse=",")), "UNSUPPORTED")
}
.adapter <- function(id) {
  cat <- .catalog(); row <- cat[cat$estimator==id,,drop=FALSE]
  function(p,options,native=FALSE) {
    .require(requireNamespace(row$package,quietly=TRUE) &&
      row$entrypoint %in% getNamespaceExports(row$package),
      paste("Install",row$package,"with",row$entrypoint), "UNSUPPORTED")
    fun <- getExportedValue(row$package,row$entrypoint)
    .require(identical(p$cluster,p$assignment_unit),
      "This adapter requires assignment_unit == cluster", "UNSUPPORTED")
    .require((p$target_mode=="native" || id=="stacked") && p$sample_mode=="native",
      "Harmonized targets/common-support intersection not certified for this adapter", "UNSUPPORTED")
    if (id=="csdid") {
      .check_options(options,c("method","control_group","ps_fit","jackknife_support","aggregation_weights"))
      jk_support <- if (is.null(options$jackknife_support)) "changing" else options$jackknife_support
      .require(is.character(jk_support) && length(jk_support)==1L && !is.na(jk_support) &&
        jk_support %in% c("changing","fixed"),"jackknife_support must be changing or fixed","UNSUPPORTED")
      .require(p$mode=="rcs","CSDID canonical package is RCS only","UNSUPPORTED")
      .require(is.null(p$event_window),"CSDID window-restricted scalar not mapped","UNSUPPORTED")
      method <- if (is.null(options$method)) "dr" else options$method
      if (identical(method,"dripw")) method <- "dr"
      .require(length(method)==1 && method %in% c("none","reg","ipw","dr"),"Unsupported CSDID method","UNSUPPORTED")
      .require(method!="none" || !length(p$covariates),"none would ignore requested covariates","UNSUPPORTED")
      ctrl <- if (is.null(options$control_group)) p$control_group else options$control_group
      if (ctrl=="native") ctrl <- "never"
      .require(ctrl %in% c("never","not_yet"),"Unsupported CSDID control group","UNSUPPORTED")
      formula <- if (method=="none") NULL else if (!length(p$covariates)) ~1 else stats::reformulate(p$covariates)
      fit <- fun(p$data,p$outcome,p$time,".didsuite_g",p$cluster,
        weights=p$weights,est_method=method,xformla=formula,control_group=ctrl,
        inference="none",aggregation=c("simple","dynamic"),
        ps_fit=if (is.null(options$ps_fit)) "warm_start" else options$ps_fit)
      weighting <- if (is.null(options$aggregation_weights)) "cohort" else options$aggregation_weights
      .require(weighting %in% c("cohort","comparison"),"Unknown CSDID aggregation weights","UNSUPPORTED")
      allow_partial <- isTRUE(p$cluster_deletion) && jk_support=="changing"
      incomplete <- any(!is.finite(fit$att_gt))
      .require(!incomplete || allow_partial,
        "One or more required CSDID ATT(g,t) cells are unidentified","IDENTIFICATION_FAILURE")
      if (weighting=="comparison" || (incomplete && allow_partial)) {
        # Reweight the returned native ATT cells, including the sampling variability of their masses.
        m <- fit$moments; map <- fit$att_map; att <- fit$att_gt
        mass <- vapply(seq_len(nrow(map)),function(k) {
          periods <- if (weighting=="comparison") match(c(map$base[k],map$post[k]),m$times) else seq_along(m$times)
          rowSums(m$weight_sum[,periods,drop=FALSE])*(m$assignment==map$group[k])
        },numeric(length(m$regions)))
        aggregate_cells <- function(use) {
          valid <- is.finite(att) & colSums(mass)>0
          .require(allow_partial || all(valid[use]),
            "One or more required CSDID ATT(g,t) cells are unidentified","IDENTIFICATION_FAILURE")
          use <- use & valid
          if (!any(use)) return(c(estimate=NA_real_,se=NA_real_))
          counts <- mass[,use,drop=FALSE]; weights <- colSums(counts)/sum(counts)
          theta <- sum(weights*att[use])
          score <- as.vector(fit$influence$cluster_att[,use,drop=FALSE] %*% weights)/m$N +
            as.vector(counts %*% (att[use]-theta))/sum(counts)
          c(estimate=unname(theta),se=sqrt(sum(score^2)))
        }
        totals <- aggregate_cells(rep(TRUE,length(att)))
        b <- unname(totals["estimate"]); se <- unname(totals["se"])
        event_ids <- sort(unique(map$event))
        events <- data.frame(event_time=event_ids,estimate=vapply(event_ids,function(e)
          unname(aggregate_cells(map$event==e)["estimate"]),numeric(1)))
      } else {
        b <- unname(fit$aggregation$estimates["simple"])
        se <- sqrt(fit$vcov$analytical["simple","simple"])
        ev <- fit$aggregation$estimates[grepl("^event:",names(fit$aggregation$estimates))]
        events <- data.frame(event_time=as.numeric(sub("^event:","",names(ev))),estimate=unname(ev))
      }
      .require(length(b)==1 && is.finite(b),"Unsupported CSDID simple ATT","IDENTIFICATION_FAILURE")
      return(list(estimate=b,control_group=ctrl,target_label=paste("CSDID RCS",method,weighting,"weighted ATT"),
        weighting_scheme=if (weighting=="comparison") "treated mass in comparison periods" else "treated-cohort sampling mass",covariate_mode=if (method=="none") "none" else "RCS pooled base/post formula; no time freezing",
        native_se=if (native) se else NA_real_,
        native_df=Inf,native_method="CSDID cluster-summed influence/N^2 normal",
        support=names(fit$att_gt)[is.finite(fit$att_gt)],jackknife_support=jk_support,event_study=events,raw=fit))
    }
    if (id %in% c("twfe","did2s","bjs","flex")) {
      .check_options(options,character())
      .require(p$control_group=="native","Explicit alternate control rule unavailable in port","UNSUPPORTED")
      .require(is.null(p$event_window),"Scalar port does not accept event_window","UNSUPPORTED")
      fit <- fun(p$data,outcome=p$outcome,time=p$time,cluster=p$cluster,gvar=".didsuite_g",
        covariates=p$covariates,weights=p$weights,unit=p$unit,data_mode=p$mode,
        inference=if (native) "native" else "point")
      fit$native_method <- switch(id,twfe="OLS CR1 t(G-1)",did2s="two-stage corrected normal",
        bjs="BJS imputation-weight normal",flex="FLEX unconditional CR1 t(G-1)")
      fit$control_group <- if (id=="twfe") "all observations; TWFE coefficient" else "untreated observations: never and not-yet"
    } else if (id %in% c("undid","didint")) {
      .check_options(options,if (id=="didint") "ccc" else "aggregation_weights")
      .require(!native,"Native clustered VCE not available","UNSUPPORTED")
      .require(is.null(p$event_window),"Window-restricted target unsupported","UNSUPPORTED")
      args <- c(list(data=p$data,outcome=p$outcome,time=p$time,cluster=p$cluster,gvar=".didsuite_g",
        covariates=p$covariates,weights=p$weights,unit=p$unit,data_mode=p$mode,
        control_group=if (p$control_group=="native") "never" else p$control_group,inference="point"),options)
      fit <- do.call(fun,args)
    } else if (id=="stacked") {
      .check_options(options,character())
      .require(length(p$event_window)==2L,"Stacked requires its full (-pre,post) window including -1","UNSUPPORTED")
      if (p$target_mode=="harmonized") {
        .require(p$event_window[1]==0 && p$event_window[2]>=0,"Harmonized Stacked expects (0,H)","UNSUPPORTED")
        p$event_window <- c(-1,p$event_window[2])
      }
      controls <- c(never="never_treated",not_yet="not_yet_treated",both="both")
      if (p$control_group=="native") p$control_group <- "both"
      .require(p$control_group %in% names(controls),"Unsupported Stacked control group","UNSUPPORTED")
      stacked_data <- p$data
      stacked_data$.didsuite_g[stacked_data$.didsuite_g==0] <- NA_real_
      fit <- fun(stacked_data,y=p$outcome,time=p$time,cluster=p$cluster,unit=p$unit,
        gvar=".didsuite_g",covariates=p$covariates,weights=p$weights,data_mode=p$mode,
        event_window=p$event_window,control_group=unname(controls[p$control_group]),
        inference=if (native) "native" else "point")
      fit$native_df <- fit$df
      fit$native_method <- "Stacked original-cluster CR1"
      fit$support <- fit$cohorts
      fit$control_group <- p$control_group
      fit$event_window_used <- p$event_window
    } else if (id=="lpdid") {
      .check_options(options,c("weighting","covariate_lag"))
      .require(p$mode=="panel" && !is.null(p$unit),"LP-DiD requires genuine panel units; RCS unsupported","UNSUPPORTED")
      .require(is.null(p$weights),"LP-DiD observation weights unsupported","UNSUPPORTED")
      .require(p$control_group=="native","LP-DiD uses its native clean-through-horizon controls","UNSUPPORTED")
      .require(length(p$event_window)==2L,"LP-DiD requires explicit pooled event_window","UNSUPPORTED")
      fit <- do.call(fun,c(list(data=p$data,y=p$outcome,time=p$time,unit=p$unit,
        cluster=p$cluster,gvar=".didsuite_g",covariates=p$covariates,
        event_window=p$event_window,data_mode=p$mode,
        inference=if (native) "native" else "point"),options))
      if (is.na(fit$native_df)) fit$native_df <- Inf
      fit$native_method <- fit$native_vce
      fit$control_group <- "native clean through pooled horizon"
      fit$weighting_scheme <- fit$weighting
      fit$covariate_mode <- fit$covariate_timing
      fit$n_used <- fit$nobs
      fit$cohorts <- fit$cohorts_used
      fit$support <- fit$cohorts_used
      fit$event_window_used <- p$event_window
    }
    .require(length(fit$estimate)==1 && is.finite(fit$estimate),"Nonfinite point estimate","IDENTIFICATION_FAILURE")
    if (is.null(fit$weighting_scheme)) fit$weighting_scheme <- if (is.null(p$weights)) "native unweighted" else "native observation weights"
    if (is.null(fit$covariate_mode)) fit$covariate_mode <- "native supplied covariates; no harness transformation"
    fit
  }
}

