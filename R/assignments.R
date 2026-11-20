# Exact set membership for RI assignment keys. R environment names are
# symbols limited to 10000 bytes; the full key must be a value, not a name.
# Short hash buckets keep lookups fast. Comparing full keys within each
# bucket prevents hash collisions from discarding distinct assignments.
.ri_key_bucket <- function(key) {
  code <- as.numeric(utf8ToInt(enc2utf8(key)))
  n <- length(code)
  index <- seq_len(n)
  forward <- sum((code * (index %% 65521 + 1)) %% 2147483647) %% 2147483647
  reverse <- sum((code * (rev(index) %% 65519 + 7)) %% 2147483629) %% 2147483629
  sprintf("h%.0f:%.0f:%d", forward, reverse, n)
}

.ri_key_has <- function(seen, key) {
  key %in% seen[[.ri_key_bucket(key)]]
}

.ri_key_add <- function(seen, key) {
  bucket <- .ri_key_bucket(key)
  keys <- seen[[bucket]]
  if (!key %in% keys) seen[[bucket]] <- c(keys, key)
  invisible(NULL)
}

# Shared canonical JSON is identical to the Python harness's v1 SHA-256 payload.
assignment_set <- function(observed, units=names(observed), reps=999L, seed=0L,
                           exact=FALSE, supplied=NULL) {
  a <- as.numeric(observed)
  if (is.null(units)) units <- as.character(seq_along(a)-1L)
  units <- as.character(units)
  .require(length(a)>0 && all(is.finite(a)) && all(a>=0) &&
             length(units)==length(a) && !anyDuplicated(units), "Invalid assignment design")
  values <- sort(unique(a)); counts <- tabulate(match(a,values), length(values))
  space <- round(exp(lgamma(length(a)+1)-sum(lgamma(counts+1))))
  .require(space>=2, "No alternative assignments", "IDENTIFICATION_FAILURE")
  key <- function(x) paste(sprintf("%.17g",x),collapse=",")
  if (inherits(supplied,"didsuite_assignments")) {
    .require(identical(supplied$units,units) && identical(supplied$observed,a),
             "Supplied assignment design differs")
    supplied <- supplied$matrix
  }
  if (!is.null(supplied)) {
    if (!is.null(colnames(supplied)))
      .require(identical(colnames(supplied),units), "Assignment column order must match units")
    draws <- as.matrix(supplied)
  } else {
    .require(length(reps)==1L && is.finite(reps) && reps>=1 && reps==floor(reps) &&
               reps<=100000, "reps must be an integer in [1,100000]")
    .require(length(seed)==1L && is.finite(seed), "Invalid seed")
    .require(!exact || space<=100000, "Exact enumeration exceeds 100000", "UNSUPPORTED")
    had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
    if (had) old <- get(".Random.seed",.GlobalEnv)
    on.exit({
      if (had) assign(".Random.seed",old,.GlobalEnv)
      else if (exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
    })
    set.seed(seed)
    n <- min(reps,space-1)
    if (space<=100000) {
      rows <- list()
      visit <- function(prefix,left) {
        if (!sum(left)) {
          if (!identical(prefix,a)) rows[[length(rows)+1L]] <<- prefix
          return(invisible(NULL))
        }
        for (j in which(left>0)) {
          nextleft <- left; nextleft[j] <- nextleft[j]-1L
          visit(c(prefix,values[j]),nextleft)
        }
      }
      visit(numeric(),counts)
      draws <- do.call(rbind,rows)
      if (!exact && n<nrow(draws)) draws <- draws[sample.int(nrow(draws),n),,drop=FALSE]
    } else {
      seen <- new.env(hash=TRUE,parent=emptyenv()); .ri_key_add(seen, key(a))
      rows <- vector("list",n); found <- 0L
      for (i in seq_len(max(1000,100*n))) {
        draw <- a[sample.int(length(a))]
        k <- key(draw)
        if (!.ri_key_has(seen, k)) {
          .ri_key_add(seen, k); found <- found+1L; rows[[found]] <- draw
          if (found==n) break
        }
      }
      .require(found==n, "Bounded unique sampler exhausted; supply assignments", "NUMERICAL_FAILURE")
      draws <- do.call(rbind,rows)
    }
  }
  .require(is.numeric(draws) && nrow(draws)>0 && ncol(draws)==length(a) &&
             all(is.finite(draws)), "Invalid assignment matrix")
  keys <- apply(draws,1,key)
  .require(!anyDuplicated(keys) && !key(a) %in% keys, "Duplicate or observed assignment")
  for (i in seq_len(nrow(draws)))
    .require(identical(sort(as.numeric(draws[i,])),sort(a)), "Changed cohort counts or dates")
  colnames(draws) <- units
  payload <- list("didsuite-assignments-v1",as.list(units),as.list(sprintf("%.17g",a)),
                  lapply(seq_len(nrow(draws)),function(i) as.list(sprintf("%.17g",draws[i,]))))
  # jsonlite's unicode escaping is made explicit for portable ASCII hashing.
  blob <- jsonlite::toJSON(payload,auto_unbox=TRUE,pretty=FALSE)
  chars <- utf8ToInt(enc2utf8(blob))
  ascii <- paste0(vapply(chars,function(ch) {
    if (ch<128) intToUtf8(ch)
    else if (ch<=65535) sprintf("\\u%04x",ch)
    else {v <- ch-65536; paste0(sprintf("\\u%04x",55296+v%/%1024),
                               sprintf("\\u%04x",56320+v%%1024))}
  },character(1)),collapse="")
  structure(list(matrix=draws,observed=a,units=units,space=space,
    hash=digest::digest(charToRaw(ascii),algo="sha256",serialize=FALSE),
    hash_algorithm="didsuite-assignments-v1/sha256",canonical_payload=ascii),
    class="didsuite_assignments")
}

