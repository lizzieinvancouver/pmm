###
#
# Definition of the corIntra correlation structure
# for PGLS.
#

Initialize.corIntra <- function(object, data, ...)
{
  ## The same as in Initialize corStruct:
  form <- formula(object)
  ## Obtaining the group information, if any
  if (!is.null(getGroupsFormula(form))) {
    attr(object, "groups") <- getGroups(object, form, data = data)
    attr(object, "Dim") <- Dim(object, attr(object, "groups"))
  } else { # no groups
    attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
  }
  ## Obtaining the covariate(s)
  attr(object, "covariate") <- getCovariate(object, data = data)
  
  ## Specific to corPhyl:
  intra <- attr(object, "intra")
  #if (is.null(data)) data <- parent.frame()
  if (is.null(rownames(data))) {
    warning("No rownames supplied in data frame, data taken to be in the same order than in matrix")
    attr(object, "index") <- 1:dim(data)[1]
  } else {
    index <- match(rownames(data), rownames(intra))
    if (any(is.na(index))) {
      warning("Rownames in data frame do not match matrix names; data taken to be in the same order as in matrix")
      attr(object, "index") <- 1:dim(data)[1]
    } else {
      attr(object, "index") <- index
    }
  }
  object
}


corIntra <- function(value, inter, intra, form = ~1, fixed = FALSE)
{
  if (value < 0 || value > 1)
    stop("the value of delta must be between 0 and 1.")
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  attr(value, "inter") <- inter
  attr(value, "intra") <- intra
  class(value) <- c("corIntra", "corPhyl", "corStruct")
  value
}

corMatrix.corIntra <- function(object, covariate = getCovariate(object), ...)
{
    if (!any(attr(object, "index")))
      stop("object has not been initialized")
    matinter <- attr(object, "inter")
    matintra <- attr(object, "intra")
    index <- attr(object, "index")
    matinter <- matinter[index, index]
    matintra <- matintra[index, index]

    mat <- object[1]*matinter+(1-object[1])*matintra
    diag(mat) <- 1
    mat
  }

coef.corIntra <- function(object, unconstrained = TRUE, ...)
{
  if (unconstrained) {
    if (attr(object, "fixed")) return(numeric(0))
    else return(object[1])
  }
  aux <- object[1]
  names(aux) <- "delta"
  aux
}