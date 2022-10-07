#2GP_functions.r
#--------------dnorm_log----------------------------------------------
# Probability density function-
dnorm_log <- function(x, mu, sigma){ #sigma is standard deviation
    p <- log(1.0/sqrt(2*pi*sigma^2)) -0.5*((x-mu)/sigma)^2
    p

}
#--------------C_ii-------------------------------------------------
# Autocovariance
C_ii <- function(par, d){
    A <- exp(par[3])
    B <- exp(par[4])
    v <- par[1]
    w <- par[2]
    C_ii_u <- pi^0.5*v^2*1.0/sqrt(A)* exp(-A*d^2*0.25)
    C_ii_v <- pi^0.5*w^2*1.0/sqrt(B)* exp(-B*d^2*0.25)

    C_ii <- C_ii_u + C_ii_v
    C_ii
}
#--------------C_12-------------------------------------------------
# Cross covariance
C_12 <- function(par, d){
    A_1 <- exp(par[3])
    A_2 <- exp(par[4])
    v_1 <- par[1]
    v_2 <- par[2]
    miu <- par[5]

    Sigma <- A_1*(A_1 + A_2)^(-1)*A_2

    C_12_u <- (2*pi)^0.5*v_1*v_2*1.0/sqrt(A_1+A_2) * exp(-Sigma*(d-miu)^2*0.5)
   
    C_12_u
}
#--------------C_21-------------------------------------------------
C_21 <- function(par, d){
    A_1 <- exp(par[3])
    A_2 <- exp(par[4])
    v_1 <- par[1]
    v_2 <- par[2]
    miu <- par[5]

    Sigma <- A_1*(A_1 + A_2)^(-1)*A_2

    C_21_u <- (2*pi)^0.5*v_1*v_2*1.0/sqrt(A_1+A_2) * exp(-Sigma*(d+miu)^2*0.5)
    C_21_u
}
#--------------Covariance-------------------------------------------------
Covariance <- function(par, t1, t2){ #parameters:; [v1 w1 f1 g1 Beta1 v2 w2 f2 g2 Beta2 mu];
    par1 <- par[1:5] 
    par2 <- par[6:10] 
    par3 <- c(par[1], par[6], par[3], par[8], par[11])  
    N1 <- length(t1)
    N2 <- length(t2)
    dim <- N1 + N2 
    I <- mat.or.vec(dim, dim) 
    Cv <- mat.or.vec(dim, dim) 
    C11 <- mat.or.vec(N1, N1)
    C22 <- mat.or.vec(N2, N2)
    C12 <- mat.or.vec(N1, N2) 
    C21 <- mat.or.vec(N2, N1)
    for (i in 1:N1){
        for (j in 1:N1){
            C11[i,j] <- C_ii(par1, t1[i]-t1[j] )
            if (i==j){
                I[i,j] <- exp(par1[5])^2;  
           }
        }
    }
    for (i in 1:N2){
        for (j in 1:N2){
            C22[i,j] <- C_ii(par2, t2[i]-t2[j])
            if (i==j){
                I[(N1+i),(N1+j)] <- exp(par2[5])^2
            }
        }
    }
    for (i in 1:N1){
        for (j in 1:N2){
            C12[i,j] <- C_12(par3, t1[i]-t2[j] )
        }
    }
    for (i in 1:N2){
        for (j in 1:N1){
            C21[i,j] <- C_21(par3, t2[i]-t1[j])
        }
    }
    Cv[1:N1, 1:N1] <- C11           
    Cv[(N1+1):dim, (N1+1):dim] <- C22 

    Cv[(N1+1):dim, 1:N1] <- C21       
    Cv[1:N1, (N1+1):dim] <- C12      

    Cv <- Cv + I
    Cv
}

#--------------log_lik-------------------------------------------------
log_lik <- function(par, Y, t1, t2, prior){

    N <- length(Y)
    C <- Covariance(par, t1, t2)
    C_I <- solve(C)
    prior_v1 <- dnorm_log(par[1], prior[1], prior[2])
    prior_v2 <- dnorm_log(par[6], prior[1], prior[2])

    prior_w1 <- dnorm_log(par[2], prior[3], prior[4])
    prior_w2 <- dnorm_log(par[7], prior[3], prior[4])

    prior_f1 <- dnorm_log(par[3], prior[5], prior[6])
    prior_f2 <- dnorm_log(par[8], prior[5], prior[6])

    prior_g1 <- dnorm_log(par[4], prior[7], prior[8])
    prior_g2 <- dnorm_log(par[9], prior[7], prior[8])

    prior_b1 <- dnorm_log(par[5], prior[9], prior[10])
    prior_b2 <- dnorm_log(par[10], prior[9], prior[10])

    prior_mu <- dnorm_log(par[11], prior[11], prior[12])
    
    M <- t(Y) %*% C_I %*% Y

    v <- prior_v1 + prior_v2
    w <- prior_w1 + prior_w2
    f <- prior_f1 + prior_f2
    g <- prior_g1 + prior_g2
    beta <- prior_b1 + prior_b2

    SVD <- svd(C) 
    A <- sum(log(SVD$d))
    log_lkh <- c(-(-0.5*A - 0.5*M - 0.5*N*log(2*pi) + v + w + f + g + beta + prior_mu))
    log_lkh
}

Optim <- function(par, fn, gr = NULL, ..., method = c("Nelder-Mead", 
"BFGS", "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf, 
control = list(), hessian = FALSE) 
{
    fn1 <- function(par) fn(par, ...)
    gr1 <- if (!is.null(gr)) 
    function(par) gr(par, ...)
    method <- match.arg(method)
    if ((length(lower) > 1L || length(upper) > 1L || lower[1L] != 
    -Inf || upper[1L] != Inf) && method != "L-BFGS-B") {
        warning("bounds can only be used with method L-BFGS-B")
        method <- "L-BFGS-B"
    }
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
    length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100L, 
    abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
    beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
    factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    nmsC <- names(con)
    if (method == "Nelder-Mead") 
    con$maxit <- 5000
    if (method == "SANN") {
        con$maxit <- 10000
        con$REPORT <- 100
    }
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (con$trace < 0) 
    warning("read the documentation for 'trace' more carefully")
    else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 
    0) 
    stop("'trace != 0' needs 'REPORT >= 1'")
    if (method == "L-BFGS-B" && any(!is.na(match(c("reltol", 
    "abstol"), namc)))) 
    warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
    npar <- length(par)
    if (npar == 1 && method == "Nelder-Mead") 
    warning("one-diml optimization by Nelder-Mead is unreliable: use optimize")
    lower <- as.double(rep(lower, , npar))
    upper <- as.double(rep(upper, , npar))
    res <- .Internal(optim(par, fn1, gr1, method, con, lower, 
    upper))
    names(res) <- c("par", "value", "counts", "convergence", 
    "message")
    nm <- names(par)
    if (!is.null(nm)) 
    names(res$par) <- nm
    names(res$counts) <- c("function", "gradient")
    if (hessian) {
        hess <- .Internal(optimhess(res$par, fn1, gr1, con))
        hess <- 0.5 * (hess + t(hess))
        if (!is.null(nm)) 
        dimnames(hess) <- list(nm, nm)
        res$hessian <- hess
    }
    res
}
