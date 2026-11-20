.abort <- function(reason, status="INVALID_INPUT") {
  stop(structure(list(message=reason, call=NULL, status=status),
                 class=c("didsuite_error","error","condition")))
}
.require <- function(ok, reason, status="INVALID_INPUT") {
  if (!isTRUE(ok)) .abort(reason, status)
}
.status <- function(e) {
  if (inherits(e, "didsuite_error")) return(e$status)
  msg <- conditionMessage(e)
  for (s in c("UNSUPPORTED","IDENTIFICATION_FAILURE","INVALID_INPUT","NUMERICAL_FAILURE"))
    if (startsWith(msg, s)) return(s)
  "NUMERICAL_FAILURE"
}
.attempt <- function(expr) tryCatch(list(status="PASS", value=expr, reason=""),
  error=function(e) list(status=.status(e), value=NULL, reason=conditionMessage(e)))
.ids <- function(x) sort(unique(as.character(x)), method="radix")
.prepare <- function(data, outcome, time, cluster, assignment_unit, unit, gvar,
                     treatment, covariates, weights, data_mode) {
  .require(is.data.frame(data) && nrow(data)>0 && !anyDuplicated(names(data)),
           "Nonempty data.frame with unique columns required")
  .require(!any(c(".didsuite_g",".didsuite_d") %in% names(data)), "Reserved column present")
  .require(!is.null(gvar) || !is.null(treatment), "Supply treatment or gvar")
  .require(!anyDuplicated(covariates), "Duplicate covariates")
  cols <- unique(c(outcome,time,cluster,assignment_unit,unit,gvar,treatment,covariates,weights))
  .require(all(cols %in% names(data)), "Requested column missing")
  .require(!anyNA(data[cols]), "Missing values: choose the sample explicitly")
  numeric <- unique(c(outcome,time,gvar,treatment,covariates,weights))
  .require(all(vapply(data[numeric], is.numeric, logical(1))) &&
             all(is.finite(as.matrix(data[numeric]))), "Model variables must be finite numeric")
  if (!is.null(weights)) .require(all(data[[weights]]>0), "Common weights must be positive")
  mode <- match.arg(data_mode, c("auto","panel","rcs"))
  if (mode=="auto") mode <- if (!is.null(unit) && anyDuplicated(data[[unit]])) "panel" else "rcs"
  if (mode=="panel") {
    .require(!is.null(unit), "Panel mode requires unit")
    .require(!anyDuplicated(data[c(unit,time)]), "Duplicate panel unit-time")
    for (col in c(cluster,assignment_unit)) .require(
      all(vapply(split(data[[col]],data[[unit]]),function(z) length(unique(z))==1L,logical(1))),
      "Panel units must nest in assignment and inference units")
  }
  ids <- .ids(data[[assignment_unit]])
  cohorts <- setNames(numeric(length(ids)),ids)
  if (!is.null(treatment)) .require(all(data[[treatment]] %in% 0:1), "Treatment must be binary")
  for (id in ids) {
    ix <- as.character(data[[assignment_unit]])==id
    if (!is.null(gvar)) {
      g <- unique(data[[gvar]][ix])
      .require(length(g)==1L && g>=0, "Cohort must be constant within assignment unit; zero means never")
    } else {
      tt <- data[[time]][ix & data[[treatment]]==1]
      g <- if (length(tt)) min(tt) else 0
    }
    expected <- as.integer(g>0 & data[[time]][ix]>=g)
    if (!is.null(treatment)) .require(all(expected==data[[treatment]][ix]),
      "Inconsistent or nonabsorbing treatment history")
    if (g>0) .require(g %in% data[[time]] && any(data[[time]][ix]<g) &&
      any(data[[time]][ix]>=g), "Treated assignment unit needs pre/post support and date in time grid")
    cohorts[id] <- g
  }
  .require(any(cohorts>0), "No treated assignment unit", "IDENTIFICATION_FAILURE")
  data$.didsuite_g <- unname(cohorts[as.character(data[[assignment_unit]])])
  data$.didsuite_d <- as.integer(data$.didsuite_g>0 & data[[time]]>=data$.didsuite_g)
  list(data=data,outcome=outcome,time=time,cluster=cluster,assignment_unit=assignment_unit,
       unit=unit,covariates=covariates,weights=weights,mode=mode,ids=ids,cohorts=cohorts)
}
.reassign <- function(p, a) {
  p$data$.didsuite_g <- a[match(as.character(p$data[[p$assignment_unit]]),p$ids)]
  p$data$.didsuite_d <- as.integer(p$data$.didsuite_g>0 & p$data[[p$time]]>=p$data$.didsuite_g)
  p$cohorts <- setNames(as.numeric(a),p$ids)
  p
}
.methods <- function(inference) {
  if (identical(inference,"point")) return(character())
  if (identical(inference,"standard")) return(c("native","jackknife"))
  if (identical(inference,"full")) return(c("native","jackknife","ri_beta","ri_t_native"))
  .require(length(inference)>0 && all(inference %in% c("native","jackknife","ri_beta","ri_t_native")),
           "Unknown inference method")
  unique(inference)
}

