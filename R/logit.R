#' @rdname logit
#' @export
#'
#' @title Perform a logit or inverse logit transformation
#'
#' @description
#' This function accepts a numeric vector and outputs the logit or inverse logit
#' of it.
#'
#' @usage
#' logit(x)
#' invlogit(x)
#'
#' @param x A numeric vector
#'
#' @return A numeric vector of the same length as that input
#'
#' @aliases logit invlogit
#'
#' @examples
#' x <- c(seq(from=0.1, to=0.9, by=0.1))
#' y <- logit(x)
logit <- function(x) {
    return(qlogis(x))
}

#' @rdname logit
#' @export
invlogit <- function(x) {
    return(plogis(x))
}
