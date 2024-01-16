CrossFit_D_Lt <- function(x, t, Folds, control) {

    if (t == x$vars$tau) {
        x$augmented[[g("lcmmtp_D_Z{t+1}")]] <- x$augmented[[x$vars$Y]]
        x$augmented[[g("lcmmtp_Q_Z{t+1}")]] <- x$augmented[[x$vars$Y]]
    }

    cfd <- list()
    for (v in 1:Folds$V) {
        Tr   <- Folds$Tr(x$augmented, v)
        Tr_a <- x$augment(Tr, t)
        P_a  <- x$augment(Folds$P(x$augmented, v), t)

        Tr_a <- Tr_a[x$observed(Tr_a, t), ]

        o <- x$observed(P_a, t, T)

        g_Mt <- g_Ast <- g_Apt <- g_Ct <- vector("numeric", nrow(P_a))

        P_a[[g("lcmmtp_Q_L{t}")]][o] <- CrossFit(
            Tr_a[Tr_a[[g("lcmmtp_med_{t}")]] == Tr_a[[x$vars$M]], ],
            P_a[o, ],
            g("lcmmtp_D_Z{t+1}"),
            c(g("lcmmtp_med_{t:x$vars$tau}"), x$vars$history("M")),
            ifelse(t == x$vars$tau, x$type, "continuous"),
            control$learners_QZ,
            control$folds_QZ
        )
        P_a[[g("lcmmtp_Q_L{t}")]] <- 0

        Tr_Ap <- x$stack_data(Folds$Tr(x$data, v), Folds$Tr(x$shifted_aprime, v), t)

        g_Apt[o] <- CrossFit(Tr_Ap[x$observed(Tr_Ap, t, T), ],
                             P_a[o, ],
                             "tmp_lcmmtp_stack_indicator",
                             c(x$vars$history("A"), x$vars$A),
                             "binomial",
                             control$learners_trt,
                             control$folds_trt)

        Tr_As <- x$stack_data(Folds$Tr(x$data, v), Folds$Tr(x$shifted_astar, v), t)

        g_Ast[o] <- CrossFit(Tr_As[x$observed(Tr_As, t, T), ],
                             P_a[o, ],
                             "tmp_lcmmtp_stack_indicator",
                             c(x$vars$history("A"), x$vars$A),
                             "binomial",
                             control$learners_trt,
                             control$folds_trt)

        # Create pooled data for g_Mt fit
        m_fit_data <- x$augment(Folds$Tr(x$data, v), t)
        m_fit_data <- m_fit_data[x$observed(m_fit_data, t, T), ]
        m_fit_data[["lcmmtp_pseudo_m_fit"]] <-
            as.numeric(m_fit_data[[g("lcmmtp_med_{t}")]] == m_fit_data[[x$vars$M[t]]])

        g_Mt[o] <- CrossFit(m_fit_data,
                            P_a[o, ],
                            "lcmmtp_pseudo_m_fit",
                            c(x$vars$history("M"), g("lcmmtp_med_{t}")),
                            "binomial",
                            control$learners_mediator,
                            control$folds_mediator)

        P_a[[g("lcmmtp_Gp_A{t}")]] <- density_ratios(g_Apt, x$observed(P_a, t))
        P_a[[g("lcmmtp_Gs_A{t}")]] <- density_ratios(g_Ast, x$observed(P_a, t))
        P_a[[g("lcmmtp_G_M{t}")]] <- G(P_a[[x$vars$M[t]]],
                                       g_Mt,
                                       P_a[[g("lcmmtp_med_{t}")]],
                                       x$observed(P_a, t))

        P_a[[g("lcmmtp_D_L{t}")]] <- D_Lt(P_a, t, x$vars$tau)
        cfd[[v]] <- P_a
    }

    cfd <- Reduce(rbind, cfd)
    data.table::setorder(cfd, "lcmmtp_ID")

    x$augmented <- cfd
}
