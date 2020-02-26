#Functions used

#Like ifelse but can take multiple conditions and outputs
ifelse_ <- function(...) {
    if (...length() %% 2 == 0) stop("ifelse_ must have an odd number of arguments: pairs of test/yes, and one no.")
    out <- ...elt(...length())
    if (...length() > 1) {
        if (!is.atomic(out) && !is.factor(out)) stop("The last entry to ifelse_ must be atomic or factor.")
        if (length(out) == 1) out <- rep(out, length(..1))
        n <- length(out)
        for (i in seq_len((...length() - 1)/2)) {
            test <- ...elt(2*i - 1)
            yes <- ...elt(2*i)
            if (length(yes) == 1) yes <- rep(yes, n)
            if (length(yes) != n || length(test) != n) stop("All entries must have the same length.")
            if (!is.logical(test)) stop(paste("The", (2*i - 1), "entry to ifelse_ must be logical."))
            if (!is.atomic(yes) && !is.factor(yes)) stop(paste("The", (2*i), "entry to ifelse_ must be atomic or factor."))
            pos <- which(test)
            out[pos] <- yes[pos]
        }
    }
    else {
        if (!is.atomic(out) && !is.factor(out)) stop("The first entry to ifelse_ must be atomic or factor.")
    }
    return(out)
}

#Get the 512 effects 
est_512_models <- function(A, Xmat, Y, models_512) {
    sapply(1:512, function(i) {
        Xmat0 <- Xmat[,models_512[i,]]
        .lm.fit(cbind(A, Xmat0, 1), y = Y)$coef[1]
    })
}

#Perform a replication of the simulation (see code for details)
est_rep <- function(rep, method = "kn", n = 200) {
    if (method == "ng") {
        #All this to replicate covariate + treatment distribution of K&N but
        #generating predictors first and then treatment to ensure DGP is
        #faithful to implied DAG (i.e., w/ confounders rather than mediators)
        G1 = sample(LETTERS[1:5], n, T, prob = c(5/50, 4/50, 32/50, 4/50, 5/50))
        X1 = runif(n, 
                   ifelse_(G1=="A", 1, G1=="B", 0, G1=="C", 1, G1=="D", 5, 0), 
                   ifelse_(G1=="A", 6, G1=="B", 1, G1=="C", 5, G1=="D", 6, 5))
        X2 = runif(n, 
                   ifelse_(G1=="A", 5, G1=="B", 1, G1=="C", 1, G1=="D", 1, 0), 
                   ifelse_(G1=="A", 6, G1=="B", 5, G1=="C", 5, G1=="D", 5, 1))
        
        #P is propensity score.
        P = 1*(X1 > 1)*(X2 > 1) - .5*(X1 > 1)*(X1 < 5)*(X2 > 1)*(X2 < 5)
        
        #Have to use sample() rather than rbinom to ensure exactly 100 T and 
        #100 C are assigned. Could use rbinom if this isn't imperative,
        #but doing so makes it harder to standardize the number of dropped
        #pairs accross replications. If using sample, P(A=1) depends on who
        #was assigned to treatment before, so P =/= PS for "last" units to be
        #assigned.
        
        #A <- rbinom(n, 1, P)
        A <- rep(0, n)
        A[P == 1] <- 1
        A[P == .5] <- sample(c(rep(1, n/2 - sum(P==1)), rep(0, n/2 - sum(P==0))), sum(P == .5), FALSE)
    }
    else {
        A <- sample(rep(0:1, each = n/2), n, FALSE)
        X1 <- ifelse(A==1, runif(n, 1,6), runif(n, 0,5))
        X2 <- ifelse(A==1, runif(n, 1,6), runif(n, 0,5))
    }
    tind <- which(A==1)
    cind <- which(A==0)
    d <- data.frame(A,X1,X2)
    
    d$Y <- 2*A + X1 + X2 + rnorm(n)
    Xmat <- cbind(X1, X2, X1^2, X2^2, X1*X1, X1^3, X2^3, X1^2*X2, X1*X2^2)
    models_512 <- as.matrix(do.call(expand.grid, lapply(1:9, function(x) c(F,T))))
    dmat <- as.matrix(d)
    
    #Container for matches
    ps_matches <- mahal_matches <- data.frame(T = tind, C = 0L, dist = 0)
    
    #Estimate effects with PSM
    ps <- glm(A ~ X1 + X2, data = d, family = binomial)$fitted #LR PS
    
    #Create distance matrix of PS between T and C
    psdistmat <- matrix(nrow = length(tind), ncol = length(cind))
    for (t in 1:nrow(psdistmat)) {
        psdistmat[t, ] <- abs(ps[tind[t]] - ps[cind])
    }
    
    #Get nearest neighbor match w/o replacement
    for (t in 1:nrow(ps_matches)) {
        #For each T (indexed by t), find C that hasn't already been matched that
        #has smallest distance from that T
        ps_matches[t, "C"] <- cind[which.min(psdistmat[t, !cind %in% ps_matches[["C"]]])]
        
        #Record what the distance was
        ps_matches[t, "dist"] <- psdistmat[t, which(cind == ps_matches[t, "C"])]
        
        #Note, starting from the top of ps_matches, but order is random because
        #ps_matches is not sorted on the PS
    }
    
    #Order ps_matches from largest to smallest distance
    ps_matches <- ps_matches[order(ps_matches$dist, decreasing = TRUE),]
    
    #Estimate effects in progressively pruned datasets
    psm_out <- lapply(0:(nrow(ps_matches) - 10), function(drop) {
        #drop is how many pairs are dropped (in K&N, indexed by how many observations
        #were dropped)
       if (drop == 0) keep_units <- rep(seq_len(nrow(d)))
       else {
         keep_units <- c(tind[tind %in% ps_matches[-c(1:drop), "T"]],
                        cind[cind %in% ps_matches[-c(1:drop), "C"]])
       }
        d0 <- d[keep_units,,drop = FALSE]
        Xmat0 <- Xmat[keep_units,,drop = FALSE]
        
        #Estimate the 512 models
        effects <- est_512_models(d0$A, Xmat0, d0$Y, models_512)
        
        #effects[1] is the simple difference in means in the matched sample
        return(c(variance = var(effects), max = max(effects), diff = effects[1],
                 X1_imbal = mean(Xmat0[d0$A==1,1]) - mean(Xmat0[d0$A==0,1]),
                 X2_imbal = mean(Xmat0[d0$A==1,2]) - mean(Xmat0[d0$A==0,2])))
    })
    
    invcovX <- solve(cov(dmat[,c("X1", "X2")]))
    mahalmat <- matrix(nrow = length(tind), ncol = length(cind))
    for (t in 1:nrow(mahalmat)) {
        mahalmat[t, ] <- sqrt(mahalanobis(dmat[cind,c("X1", "X2")], dmat[tind[t],c("X1", "X2")],
                                          cov = invcovX, inverted = TRUE))
    }
    for (t in 1:nrow(mahal_matches)) {
        mahal_matches[t, "C"] <- cind[which.min(mahalmat[t, !cind %in% mahal_matches[["C"]]])]
        mahal_matches[t, "dist"] <- mahalmat[t, which(cind == mahal_matches[t, "C"])]
    }
    mahal_matches <- mahal_matches[order(mahal_matches$dist, decreasing = TRUE),]
    mahal_out <- lapply(0:(nrow(mahal_matches) - 10), function(drop) {
        if (drop == 0) keep_units <- rep(seq_len(nrow(d)))
        else {
            keep_units <- c(tind[tind %in% mahal_matches[-c(1:drop), "T"]],
                        cind[cind %in% mahal_matches[-c(1:drop), "C"]])
        }
        d0 <- d[keep_units,,drop = FALSE]
        Xmat0 <- Xmat[keep_units,,drop = FALSE]
        
        effects <- est_512_models(d0$A, Xmat0, d0$Y, models_512)
        return(c(variance = var(effects), max = max(effects), diff = effects[1],
                 X1_imbal = mean(Xmat0[d0$A==1,1]) - mean(Xmat0[d0$A==0,1]),
                 X2_imbal = mean(Xmat0[d0$A==1,2]) - mean(Xmat0[d0$A==0,2])))
    })
    return(list(psm_out = psm_out, mahal_out = mahal_out))
}

get_out <- function(...) {
    dplyr::bind_rows(lapply(list(...), function(res) {
        dplyr::bind_rows(setNames(lapply(c("psm_out", "mahal_out"), function(m) {
            data.frame(ind = seq_along(res[[1]][[m]]),
                       variance = sapply(seq_along(res[[1]][[m]]), function(dr) {
                           mean(sapply(res, function(r) r[[m]][[dr]]["variance"]))
                       }),
                       max_est = sapply(seq_along(res[[1]][[m]]), function(dr) {
                           mean(sapply(res, function(r) r[[m]][[dr]]["max"]))
                       }),
                       est = sapply(seq_along(res[[1]][[m]]), function(dr) {
                           mean(sapply(res, function(r) r[[m]][[dr]]["diff"]))
                       }),
                       X1_imbal = sapply(seq_along(res[[1]][[m]]), function(dr) {
                           mean(sapply(res, function(r) r[[m]][[dr]]["X1_imbal"]))
                       }),
                       X2_imbal = sapply(seq_along(res[[1]][[m]]), function(dr) {
                           mean(sapply(res, function(r) r[[m]][[dr]]["X2_imbal"]))
                       }))
        }), c("psm_out", "mahal_out")), .id = "method")
    }), .id = "data_method")
}