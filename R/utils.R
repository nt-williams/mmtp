g <- glue::glue

G <- function(var, g, level, cens = TRUE) {
    out <- (as.numeric(var == level) * cens) / g #b(g)
    out[is.na(out)] <- 0
    out
}

density_ratios <- function(pred, cens = TRUE) {
    (pred * cens) / (1 - pmin(pred, 0.999))
}

K_p <- function(G, l, u) {
    if (l > u) {
        return(rep(1, nrow(G)))
    }

    out <- apply(as.matrix(G[, g("lcmmtp_Gp_A{l:u}")]), 1, prod)
    out
}

K_s <- function(G, l, u) {
    if (l > u) {
        return(rep(1, nrow(G)))
    }

    out <- apply(as.matrix(G[, g("lcmmtp_Gs_A{l:u}")]), 1, prod)
    out
}

H <- function(G, l, u) {
    if (l > u) {
        return(rep(1, nrow(G)))
    }

    out <- apply(as.matrix(G[, g("lcmmtp_G_M{l:u}")]), 1, prod)
    out
}

Sum <- function(x) Reduce(`+`, x)

b <- function(x) {
    pmax(x, .01)
}
