# Multi-state model
# https://github.com/cran/msm/blob/master/R/msm.R

msm <- function (formula, subject = NULL, data = list(), qmatrix, gen.inits = FALSE, 
          obstype = NULL, 
          covariates = NULL, deathexact = NULL, 
          exacttimes = FALSE, cl = 0.95, 
          center = TRUE, opt.method = "optim", hessian = NULL, 
          use.deriv = TRUE, use.expm = TRUE, analyticp = TRUE, na.action = na.omit, 
          ...) 
{
  call <- match.call()
  
  subj <- eval(substitute(subject), data, parent.frame())
  qmatrix <- crudeinits.msm(formula, subj, qmatrix, 
                            data, censor = NULL, censor.states = NULL)
  
  qmodel <- qmodel.orig <- msm.form.qmodel(qmatrix, qconstraint, 
                                           analyticp, use.expm, phase.states)
  
  #
  hmodel <- list(hidden = FALSE, models = rep(0, qmodel$nstates), 
                 nipars = 0, nicoveffs = 0, totpars = 0, ncoveffs = 0)
  emodel <- list(misc = FALSE, npars = 0, ndpars = 0, nipars = 0, 
                 nicoveffs = 0)
  hmodel$ematrix <- FALSE
  dmodel <- msm:::msm.form.dmodel(death=NULL, qmodel, hmodel)
  cmodel <- msm:::msm.form.cmodel(censor = NULL, censor.states = NULL, qmodel$qmatrix)
  #
  
  ## COVARIATES
  if (is.list(covariates)) {
    covlist <- covariates
    msm.check.covlist(covlist, qmodel)
    ter <- lapply(covlist, function(x) attr(terms(x), "term.labels"))
    covariates <- reformulate(unique(unlist(ter)))
  }
  else covlist <- NULL
  
  ## CREATE MODEL FRAME
  indx <- match(c("data", "subject", "obstrue"), names(call), 
                nomatch = 0)
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  temp[["state"]] <- as.name(all.vars(formula[[2]]))
  temp[["time"]] <- as.name(all.vars(formula[[3]]))
  
  varnames <- function(x) {
    attr(terms(x), "term.labels")
  }

  covnames <- varnames(covariates)
  
  temp[["formula"]] <- if (length(covnames) > 0) 
    reformulate(covnames)
  else ~1
  temp[["na.action"]] <- na.pass
  temp[["data"]] <- data
  mf <- eval(temp, parent.frame())
  usernames <- c(state = all.vars(formula[[2]]), 
                 time = all.vars(formula[[3]]), 
                 subject = as.character(temp$subject), 
                 obstype = as.character(substitute(obstype)), 
                 obstrue = as.character(temp$obstrue))
  
  attr(mf, "usernames") <- usernames
  
  indx <- match(c("formula", "data"), names(call), nomatch = 0)
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  mfst <- eval(temp, parent.frame())
  if (is.matrix(mfst[[1]]) && !is.matrix(mf$"(state)")) 
    mf$"(state)" <- mfst[[1]]
  if (is.factor(mf$"(state)")) {
    if (!all(grepl("^[[:digit:]]+$", as.character(mf$"(state)")))) 
      stop("state variable should be numeric or a factor with ordinal numbers as levels")
    else mf$"(state)" <- as.numeric(as.character(mf$"(state)"))
  }
  # msm.check.state(qmodel$nstates, mf$"(state)", cmodel$censor, 
  #                 hmodel)
  if (is.null(mf$"(subject)")) 
    mf$"(subject)" <- rep(1, nrow(mf))
  msm.check.times(mf$"(time)", mf$"(subject)", mf$"(state)")
  obstype <- if (missing(obstype)) 
    NULL
  else eval(substitute(obstype), data, parent.frame())
  mf$"(obstype)" <- msm:::msm.form.obstype(mf, obstype, dmodel, exacttimes)
  mf$"(obstrue)" <- msm:::msm.form.obstrue(mf, hmodel, cmodel)
  mf$"(obs)" <- seq_len(nrow(mf))
  basenames <- c("(state)", "(time)", "(subject)", "(obstype)", 
                 "(obstrue)", "(obs)")
  attr(mf, "covnames") <- setdiff(names(mf), basenames)
  attr(mf, "covnames.q") <- colnames(attr(terms(covariates), 
                                          "factors"))

  attr(mf, "ncovs") <- length(attr(mf, "covnames"))
  
  ## initcovariates
  ic <- all.vars(initcovariates)
  others <- c(covariates, misccovariates, hcovariates)
  oic <- ic[!ic %in% unlist(lapply(others, all.vars))]
  attr(mf, "icovi") <- match(oic, colnames(mf))
  if (missing(na.action) || identical(na.action, na.omit) || 
      (identical(na.action, "na.omit"))) 
    mf <- msm:::na.omit.msmdata(mf)
  else if (identical(na.action, na.fail) || (identical(na.action, 
                                                       "na.fail"))) 
    mf <- na.fail.msmdata(mf)
  else stop("na.action should be \"na.omit\" or \"na.fail\"")
  
  
  attr(mf, "npts") <- length(unique(mf$"(subject)"))
  attr(mf, "ntrans") <- nrow(mf) - attr(mf, "npts")


  covnames <- varnames(covariates)
  if (length(covnames) > 0) {
    mm.mean <- model.matrix(reformulate(covnames), mf)
    cm <- colMeans(mm.mean[duplicated(mf$"(subject)", fromLast = TRUE), 
                           , drop = FALSE], na.rm = T)
    cm["(Intercept)"] <- 0
  }
  else cm <- NULL
  attr(mf, "covmeans") <- cm
  mf.agg <- msm:::msm.form.mf.agg(list(data = list(mf = mf), qmodel = qmodel, 
                                 hmodel = hmodel, cmodel = cmodel))
  mf$"(pcomb)" <- msm:::msm.form.hmm.agg(mf)
  
  if (inherits(misccovariates, "formula")) {
    if (!emodel$misc) 
      stop("misccovariates supplied but no ematrix")
    hcovariates <- lapply(ifelse(rowSums(emodel$imatrix) > 
                                   0, deparse(misccovariates), deparse(~1)), as.formula)
  }
  mm.cov <- msm:::msm.form.mm.cov(list(data = list(mf = mf), covariates = covariates, 
                                 center = center))
  mm.cov.agg <- msm:::msm.form.mm.cov.agg(list(data = list(mf.agg = mf.agg), 
                                         covariates = covariates, hmodel = hmodel, cmodel = cmodel, 
                                         center = center))
  mm.mcov <- msm:::msm.form.mm.mcov(list(data = list(mf = mf), misccovariates = misccovariates, 
                                   emodel = emodel, center = center))
  mm.hcov <- msm:::msm.form.mm.hcov(list(data = list(mf = mf), hcovariates = hcovariates, 
                                   qmodel = qmodel, hmodel = hmodel, center = center))
  mm.icov <- msm:::msm.form.mm.icov(list(data = list(mf = mf), initcovariates = initcovariates, 
                                   hmodel = hmodel, center = center))
  if (!is.null(covlist)) {
    cri <- msm:::msm.form.cri(covlist, qmodel, mf, mm.cov)
  }
  else cri <- NULL
  qcmodel <- if (ncol(mm.cov) > 1) 
    msm:::msm.form.covmodel(mf, mm.cov, constraint, covinits, cm, 
                      qmodel$npars, cri)
  else {
    if (!is.null(constraint)) 
      warning("constraint specified but no covariates")
    list(npars = 0, ncovs = 0, ndpars = 0)
  }
  if (!emodel$misc || is.null(misccovariates)) 
    ecmodel <- list(npars = 0, ncovs = 0)
  if (!is.null(misccovariates)) {
    if (!emodel$misc) {
      warning("misccovariates have been specified, but misc is FALSE. Ignoring misccovariates.")
    }
    else {
      ecmodel <- msm.form.covmodel(mf, mm.mcov, miscconstraint, 
                                   misccovinits, cm, emodel$npars, cri = NULL)
      hcovariates <- msm.misccov2hcov(misccovariates, emodel)
      hcovinits <- msm.misccovinits2hcovinits(misccovinits, 
                                              hcovariates, emodel, ecmodel)
    }
  }
  if (!is.null(hcovariates)) {
    if (hmodel$mv) 
      stop("hcovariates not supported for multivariate hidden Markov models")
    hmodel <- msm.form.hcmodel(hmodel, mm.hcov, hcovinits, 
                               hconstraint)
    if (emodel$misc) 
      hmodel$covconstr <- msm.form.hcovconstraint(miscconstraint, 
                                                  hmodel)
  }
  else if (hmodel$hidden) {
    npars <- if (hmodel$mv) 
      colSums(hmodel$npars)
    else hmodel$npars
    hmodel <- c(hmodel, list(ncovs = rep(rep(0, hmodel$nstates), 
                                         npars), ncoveffs = 0))
    class(hmodel) <- "hmodel"
  }
  if (!is.null(initcovariates)) {
    if (hmodel$hidden) 
      hmodel <- msm.form.icmodel(hmodel, mm.icov, initcovinits)
    else warning("initprobs and initcovariates ignored for non-hidden Markov models")
  }
  else if (hmodel$hidden) {
    hmodel <- c(hmodel, list(nicovs = rep(0, hmodel$nstates - 
                                            1), nicoveffs = 0, cri = ecmodel$cri))
    class(hmodel) <- "hmodel"
  }
  if (hmodel$hidden && !emodel$misc) {
    hmodel$constr <- msm.form.hconstraint(hconstraint, hmodel)
    hmodel$covconstr <- msm.form.hcovconstraint(hconstraint, 
                                                hmodel)
  }
  if (hmodel$hidden) 
    hmodel$ranges <- msm.form.hranges(hranges, hmodel)
  if (hmodel$hidden) 
    hmodel <- msm.form.initprobs(hmodel, initprobs, mf)
  
  
  p <- msm:::msm.form.params(qmodel, qcmodel, emodel, hmodel, fixedpars)
  msmdata <- list(mf = mf, mf.agg = mf.agg, mm.cov = mm.cov, 
                  mm.cov.agg = mm.cov.agg, mm.mcov = mm.mcov, mm.hcov = mm.hcov, 
                  mm.icov = mm.icov)
  if (p$fixed) 
    opt.method <- "fixed"
  if (is.null(hessian)) 
    hessian <- !p$fixed
  p <- msm:::msm.optim(opt.method, p, hessian, use.deriv, msmdata, 
                 qmodel, qcmodel, cmodel, hmodel, ...)
  if (p$fixed) {
    p$foundse <- FALSE
    p$covmat <- NULL
  }
  else {
    p$params <- msm.rep.constraints(p$params, p, hmodel)
    hess <- if (hessian) 
      p$opt$hessian
    else p$information
    if (!is.null(hess) && all(!is.na(hess)) && all(!is.nan(hess)) && 
        all(is.finite(hess)) && all(eigen(hess)$values > 
                                    0)) {
      p$foundse <- TRUE
      p$covmat <- matrix(0, nrow = p$npars, ncol = p$npars)
      p$covmat[p$optpars, p$optpars] <- solve(0.5 * hess)
      p$covmat <- p$covmat[!duplicated(abs(p$constr)), 
                           !duplicated(abs(p$constr)), drop = FALSE][abs(p$constr), 
                                                                     abs(p$constr), drop = FALSE]
      p$ci <- cbind(p$params - qnorm(1 - 0.5 * (1 - cl)) * 
                      sqrt(diag(p$covmat)), p$params + qnorm(1 - 0.5 * 
                                                               (1 - cl)) * sqrt(diag(p$covmat)))
      p$ci[p$fixedpars, ] <- NA
      for (i in 1:2) p$ci[, i] <- gexpit(p$ci[, i], p$ranges[, 
                                                             "lower", drop = FALSE], p$ranges[, "upper", drop = FALSE])
    }
    else {
      p$foundse <- FALSE
      p$covmat <- p$ci <- NULL
      if (!is.null(hess)) 
        warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.")
    }
  }
  p$estimates.t <- p$params
  p$estimates.t <- msm.inv.transform(p$params, hmodel, p$ranges)
  if (any(p$plabs == "p") && p$foundse) {
    p.se <- p.se.msm(x = list(data = msmdata, qmodel = qmodel, 
                              emodel = emodel, hmodel = hmodel, qcmodel = qcmodel, 
                              ecmodel = ecmodel, paramdata = p, center = center), 
                     covariates = if (center) 
                       "mean"
                     else 0)
    p$ci[p$plabs %in% c("p", "pbase"), ] <- as.numeric(unlist(p.se[, 
                                                                   c("LCL", "UCL")]))
  }
  if (p$foundse && any(p$plabs == "initp")) 
    p <- initp.ci.msm(p, cl)
  msmobject <- list(call = match.call(), minus2loglik = p$lik, 
                    deriv = p$deriv, estimates = p$params, estimates.t = p$estimates.t, 
                    fixedpars = p$fixedpars, center = center, covmat = p$covmat, 
                    ci = p$ci, opt = p$opt, foundse = p$foundse, data = msmdata, 
                    qmodel = qmodel, emodel = emodel, qcmodel = qcmodel, 
                    ecmodel = ecmodel, hmodel = hmodel, cmodel = cmodel, 
                    pci = pci, paramdata = p, cl = cl, covariates = covariates, 
                    misccovariates = misccovariates, hcovariates = hcovariates, 
                    initcovariates = initcovariates)
  attr(msmobject, "fixed") <- p$fixed
  class(msmobject) <- "msm"
  msmobject <- msm.form.output(msmobject, "intens")
  q <- qmatrix.msm(msmobject, covariates = (if (center) 
    "mean"
    else 0))
  msmobject$Qmatrices$baseline <- q$estimates
  msmobject$QmatricesSE$baseline <- q$SE
  msmobject$QmatricesL$baseline <- q$L
  msmobject$QmatricesU$baseline <- q$U
  if (hmodel$hidden) {
    msmobject$hmodel <- msm.form.houtput(hmodel, p, msmdata)
  }
  if (emodel$misc) {
    msmobject <- msm.form.output(msmobject, "misc")
    e <- ematrix.msm(msmobject, covariates = (if (center) 
      "mean"
      else 0))
    msmobject$Ematrices$baseline <- e$estimates
    msmobject$EmatricesSE$baseline <- e$SE
    msmobject$EmatricesL$baseline <- e$L
    msmobject$EmatricesU$baseline <- e$U
  }
  msmobject$msmdata[!names(msmobject$msmdata) == "mf"] <- NULL
  msmobject$sojourn <- sojourn.msm(msmobject, covariates = (if (center) 
    "mean"
    else 0))
  msmobject
}





### FORM QMODEL ####
msm.form.qmodel <- function(qmatrix, qconstraint=NULL, analyticp=TRUE, use.expm=FALSE, phase.states=NULL)
{
  #msm.check.qmatrix(qmatrix)
  nstates <- dim(qmatrix)[1]
  qmatrix <- msm.fixdiag.qmatrix(qmatrix)
  if (is.null(rownames(qmatrix)))
    rownames(qmatrix) <- colnames(qmatrix) <- paste("State", seq(nstates))
  else if (is.null(colnames(qmatrix))) colnames(qmatrix) <- rownames(qmatrix)
  imatrix <- ifelse(qmatrix > 0, 1, 0)
  inits <- t(qmatrix)[t(imatrix)==1]
  npars <- sum(imatrix)
  ## for phase-type models, leave processing qconstraint until after phased Q matrix has been formed in msm.phase2qmodel
  if (!is.null(qconstraint) && is.null(phase.states)) { 
    if (!is.numeric(qconstraint)) stop("qconstraint should be numeric")
    if (length(qconstraint) != npars)
      stop("baseline intensity constraint of length " ,length(qconstraint), ", should be ", npars)
    constr <- match(qconstraint, unique(qconstraint))
  }
  else
    constr <- 1:npars
  ndpars <- max(constr)
  ipars <- t(imatrix)[t(lower.tri(imatrix) | upper.tri(imatrix))]
  graphid <- paste(which(ipars==1), collapse="-")
  
  if (analyticp && graphid %in% names(msm:::.msm.graphs[[paste(nstates)]])) {
    iso <- msm:::.msm.graphs[[paste(nstates)]][[graphid]]$iso
    perm <- msm:::.msm.graphs[[paste(nstates)]][[graphid]]$perm
    qperm <- order(perm)
  }
  else {
    iso <- 0
    perm <- qperm <- NA
  }
  qmodel <- list(nstates=nstates, iso=iso, perm=perm, qperm=qperm,
                 npars=npars, imatrix=imatrix, qmatrix=qmatrix, inits=inits,
                 constr=constr, ndpars=ndpars, expm=as.numeric(use.expm))
  class(qmodel) <- "msmqmodel"
  qmodel
}

### FIX DIAg QMATRIX ####
msm.fixdiag.qmatrix <- function(qmatrix)
{
  diag(qmatrix) <- 0
  diag(qmatrix) <- - rowSums(qmatrix)
  qmatrix
}

### DEATH MODEL ####
msm.form.dmodel <- function (death, qmodel, hmodel) 
{
  nstates <- qmodel$nstates
  statelist <- if (nstates == 2) 
    "1, 2"
  else if (nstates == 3) 
    "1, 2, 3"
  else paste("1, 2, ... ,", nstates)
  if (is.null(death)) 
    death <- FALSE
  if (is.logical(death) && death == TRUE) 
    states <- nstates
  else if (is.logical(death) && death == FALSE) 
    states <- numeric(0)
  else if (!is.numeric(death)) 
    stop("Exact death states indicator must be numeric")
  else if (length(setdiff(death, 1:nstates)) > 0) 
    stop("Exact death states indicator contains states not in ", 
         statelist)
  else states <- death
  ndeath <- length(states)
  if (hmodel$hidden) {
    mods <- if (is.matrix(hmodel$models)) 
      hmodel$models[1, ]
    else hmodel$models
    if (!all(mods[states] == match("identity", .msm.HMODELS))) 
      stop("Exact death states should have the identity hidden distribution hmmIdent()")
    obs <- ifelse(hmodel$npars[states] > 0, hmodel$pars[hmodel$parstate %in% 
                                                          states], states)
  }
  else obs <- states
  if (any(states %in% transient.msm(qmatrix = qmodel$qmatrix))) 
    stop("Not all the \"death\" states are absorbing states")
  list(ndeath = ndeath, states = states, obs = obs)
}

### CHECK COVLIST ####
msm.check.covlist <- function (covlist, qemodel) 
{
  check.numnum <- function(str) length(grep("^[0-9]+-[0-9]+$", 
                                            str)) == length(str)
  num <- sapply(names(covlist), check.numnum)
  if (!all(num)) {
    badnums <- which(!num)
    plural1 <- if (length(badnums) > 1) 
      "s"
    else ""
    plural2 <- if (length(badnums) > 1) 
      "e"
    else ""
    badnames <- paste(paste("\"", names(covlist)[badnums], 
                            "\"", sep = ""), collapse = ",")
    badnums <- paste(badnums, collapse = ",")
    stop("Name", plural1, " ", badnames, " of \"covariates\" formula", 
         plural2, " ", badnums, " not in format \"number-number\"")
  }
  for (i in seq_along(covlist)) if (!inherits(covlist[[i]], 
                                              "formula")) 
    stop("\"covariates\" should be a formula or list of formulae")
  trans <- sapply(strsplit(names(covlist), "-"), as.numeric)
  tm <- if (inherits(qemodel, "msmqmodel")) 
    "transition"
  else "misclassification"
  qe <- if (inherits(qemodel, "msmqmodel")) 
    "qmatrix"
  else "ematrix"
  imat <- qemodel$imatrix
  for (i in seq(length = ncol(trans))) {
    if (imat[trans[1, i], trans[2, i]] != 1) 
      stop("covariates on ", names(covlist)[i], " ", tm, 
           " requested, but this is not permitted by the ", 
           qe, ".")
  }
}

### CHECK TIME ####

msm.check.times <- function (time, subject, state = NULL) 
{
  final.rows <- !is.na(subject) & !is.na(time)
  if (!is.null(state)) {
    nas <- if (is.matrix(state)) 
      apply(state, 1, function(x) all(is.na(x)))
    else is.na(state)
    final.rows <- final.rows & !nas
    state <- if (is.matrix(state)) 
      state[final.rows, , drop = FALSE]
    else state[final.rows]
  }
  final.rows <- which(final.rows)
  time <- time[final.rows]
  subject <- subject[final.rows]
  subj.num <- match(subject, unique(subject))
  nobspt <- table(subj.num)
  if (any(nobspt == 1)) {
    badsubjs <- unique(subject)[nobspt == 1]
    andothers <- if (length(badsubjs) > 3) 
      " and others"
    else ""
    if (length(badsubjs) > 3) 
      badsubjs <- badsubjs[1:3]
    badlist <- paste(badsubjs, collapse = ", ")
    plural <- if (length(badsubjs) == 1) 
      ""
    else "s"
    has <- if (length(badsubjs) == 1) 
      "has"
    else "have"
    warning("Subject", plural, " ", badlist, andothers, " only ", 
            has, " one complete observation")
  }
  ind <- tapply(seq_along(subj.num), subj.num, length)
  imin <- tapply(seq_along(subj.num), subj.num, min)
  imax <- tapply(seq_along(subj.num), subj.num, max)
  adjacent <- (ind == imax - imin + 1)
  if (any(!adjacent)) {
    badsubjs <- unique(subject)[!adjacent]
    andothers <- if (length(badsubjs) > 3) 
      " and others"
    else ""
    if (length(badsubjs) > 3) 
      badsubjs <- badsubjs[1:3]
    badlist <- paste(badsubjs, collapse = ", ")
    plural <- if (length(badsubjs) == 1) 
      ""
    else "s"
    stop("Observations within subject", plural, " ", badlist, 
         andothers, " are not adjacent in the data")
  }
  orderedpt <- !tapply(time, subj.num, is.unsorted)
  if (any(!orderedpt)) {
    badsubjs <- unique(subject)[!orderedpt]
    andothers <- if (length(badsubjs) > 3) 
      " and others"
    else ""
    if (length(badsubjs) > 3) 
      badsubjs <- badsubjs[1:3]
    badlist <- paste(badsubjs, collapse = ", ")
    plural <- if (length(badsubjs) == 1) 
      ""
    else "s"
    stop("Observations within subject", plural, " ", badlist, 
         andothers, " are not ordered by time")
  }
  if (!is.null(state)) {
    if (is.matrix(state)) 
      state <- apply(state, 1, paste, collapse = ",")
    prevsubj <- c(-Inf, subj.num[seq_along(subj.num) - 1])
    prevtime <- c(-Inf, time[1:length(time) - 1])
    prevstate <- c(-Inf, state[1:length(state) - 1])
    sametime <- final.rows[subj.num == prevsubj & prevtime == 
                             time & prevstate != state]
    badlist <- paste(paste(sametime - 1, sametime, sep = " and "), 
                     collapse = ", ")
    if (any(sametime)) 
      warning("Different states observed at the same time on the same subject at observations ", 
              badlist)
  }
  invisible()
}

### OBSERVATION TYPE ####
msm.form.obstype <- function (mf, obstype, dmodel, exacttimes) 
{
  if (!is.null(obstype)) {
    if (!is.numeric(obstype)) 
      stop("obstype should be numeric")
    if (length(obstype) == 1) 
      obstype <- rep(obstype, nrow(mf))
    else if (length(obstype) != nrow(mf)) 
      stop("obstype of length ", length(obstype), ", should be length 1 or ", 
           nrow(mf))
    if (any(!obstype[duplicated(mf$"(subject)")] %in% 1:3)) 
      stop("elements of obstype should be 1, 2, or 3")
  }
  else if (!is.null(exacttimes) && exacttimes) 
    obstype <- rep(2, nrow(mf))
  else {
    obstype <- rep(1, nrow(mf))
    if (dmodel$ndeath > 0) {
      dobs <- if (is.matrix(mf$"(state)")) 
        apply(mf$"(state)", 1, function(x) any(x %in% 
                                                 dmodel$obs))
      else (mf$"(state)" %in% dmodel$obs)
      obstype[dobs] <- 3
    }
  }
  obstype
}

msm.form.obstrue <- function (mf, hmodel, cmodel) 
{
  obstrue <- mf$"(obstrue)"
  if (!is.null(obstrue)) {
    if (!hmodel$hidden) {
      warning("Specified obstrue for a non-hidden model, ignoring.")
      obstrue <- rep(1, nrow(mf))
    }
    else if (!is.numeric(obstrue) && !is.logical(obstrue)) 
      stop("obstrue should be logical or numeric")
    else {
      if (is.logical(obstrue) || (all(na.omit(obstrue) %in% 
                                      0:1) && !any(is.na(obstrue)))) {
        if (!is.null(ncol(mf$"(state)")) && ncol(mf$"(state)") > 
            1) 
          stop("obstrue must contain NA or the true state for a multiple outcome HMM, not an 0/1 indicator")
        obstrue <- ifelse(obstrue, mf$"(state)", 0)
      }
      else {
        if (!all(na.omit(obstrue) %in% 0:hmodel$nstates)) {
          stop("Interpreting \"obstrue\" as containing true states, but it contains values not in 0,1,...,", 
               hmodel$nstates)
        }
        obstrue[is.na(obstrue)] <- 0
      }
    }
  }
  else if (hmodel$hidden) 
    obstrue <- rep(0, nrow(mf))
  else obstrue <- mf$"(state)"
  if (cmodel$ncens > 0) {
    for (i in seq_along(cmodel$censor)) obstrue[obstrue == 
                                                  cmodel$censor[i] & obstrue > 0] <- cmodel$states[cmodel$index[i]]
  }
  obstrue
}
