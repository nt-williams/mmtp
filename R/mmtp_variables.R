#' R6 class for an lcmmtp_variables
#'
#' @export
mmtp_variables <- R6::R6Class(
    "mmtp_variables",
    public = list(
        W = NULL,
        L = NULL,
        A = NULL,
        M = NULL,
        Z = NULL,
        Y = NULL,
        cens = NULL,
        tau = 1,
        initialize = function(W, A, Z, M, Y, cens) {
            checkmate::assertCharacter(A, len = 1)
            checkmate::assertCharacter(Y, len = 1)

            self$A <- A
            self$Y <- Y

            checkmate::assertList(W, types = "character", len = 2)
            self$W <- W

            if (!missing(cens)) {
                checkmate::assertCharacter(cens, len = 1)
                self$cens <- cens
            }

            checkmate::assertCharacter(M)
            checkmate::assertCharacter(Z)

            self$Z <- Z
            self$M <- M

            invisible(self)
        },
        #' Get all parent nodes for a variable
        history = function(var = c("A", "Z", "M", "Y")) {
            switch(
                match.arg(var),
                A = private$parents_A(),
                Z = private$parents_Z(),
                M = private$parents_M(),
                Y = private$parents_Y()
            )
        },
        all_vars = function() {
            unique(c(self$W$trt, self$W$outcome, self$A, self$Z, self$M, self$Y, self$cens))
        }
    ),
    private = list(
        parents_A = function() {
            self$W$trt
        },
        parents_Z = function() {
            unique(c(self$W$trt, self$W$outcome, self$A))
        },
        parents_M = function() {
            unique(c(self$W$trt, self$W$outcome, self$A, self$Z))
        },
        parents_Y = function() {
            unique(c(self$W$trt, self$W$outcome, self$A, self$Z, self$M))
        }
    )
)
