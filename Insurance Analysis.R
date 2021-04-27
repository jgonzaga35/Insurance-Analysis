# Preparation
data <- read.csv("Loss.csv", header=TRUE)
claims <- as.vector(data$X6.09610025985779)

# Install Packages
install.packages("stats4")
install.packages("MASS")
install.packages("fitdistrplus")
install.packages("ggplot2")
install.packages("fgof")
install.packages("dgof")
install.packages("actuar")
require(stats4)
require(MASS)
require(fitdistrplus)
require(ggplot2)
require(fgof)
require(dgof)
require(actuar)

# Use MLE to fit an appropriate accident severity
# distribution for individual claims. 
# Estimate the model parameters for a given model and
# present the fitted model
fitW <- fitdist(claims, "weibull", method = "mle")
fitln <- fitdist(claims, "lnorm", method = "mle")
fitg <- fitdist(claims, "gamma", method = "mle")
fitP <- fitdist(claims, "pareto", method = "mle", 
                start=list(scale=0.01,shape=500))
# Present estimates
fitW$estimate
fitln$estimate
fitg$estimate
fitP$estimate

# Graphs
cdfcomp(list(fitW, fitg, fitln, fitP), 
        legendtext=c("Weibull", "Gamma", "Lognormal", "Pareto"),
        plotstyle = "ggplot")+geom_line(linetype="dashed",size=2)+
    ggtitle("Empirical and Theoretical CDFS")+labs(x="Claims")+
    theme(legend.position = c(0.8,0.4),legend.text = element_text(size=15))+
    theme(plot.title = element_text(size=16,hjust=0.5,face="bold"))+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
denscomp(breaks=100, list(fitW, fitg, fitln, fitP), 
         legendtext=c("Weibull", "Gamma", "Lognormal", "Pareto"),
         plotstyle = "ggplot")+geom_line(linetype="dashed",size=2)+
    theme_bw()+ggtitle("Histogram of Claims")+
    theme(legend.position = c(0.8,0.4),legend.text = element_text(size=15))+
    labs(x="Claims",size=15)+
    theme(plot.title = element_text(size=16,hjust=0.5,face="bold"))+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
qqcomp(list(fitW, fitg, fitln, fitP), 
       legendtext=c("Weibull", "Gamma", "Lognormal", "Pareto"),
       plotstyle="ggplot")+theme_bw()+theme(legend.position = c(0.8,0.4),
                                            legend.text = element_text(size=15))+
    theme(plot.title = element_text(size=16,hjust=0.5,face="bold"))+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
ppcomp(list(fitW, fitg, fitln, fitP), 
       legendtext=c("Weibull", "Gamma", "Lognormal", "Pareto"),
       plotstyle="ggplot")+theme_bw()+theme(legend.position = c(0.8,0.4),
                                            legend.text = element_text(size=15))+
    theme(plot.title = element_text(size=16,hjust=0.5,face="bold"))+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# Test statistics
gofstat(list(fitW, fitg, fitln, fitP), 
        fitnames=c("Weibull", "Gamma", "Lognormal", "Pareto"))

##################
# TASK 2         #
##################
# NO REINSURANCE #
##################
# Approximate the ruin probability within 5 years
install.packages("ruin")
library(ruin)
install.packages("actuar")
library(actuar)
install.packages("heavy")
library(heavy)

## Helper function to determine if a surplus process has
## dropped below 0 at any point in time
is_ruin <- function(X){
    if (length(X[X<0]) >= 1) {
        return(1)
    } else {
        return(0)
    }
}

## Simulates surplus process for no reinsurance 
simulateRuin <-
    function(nsim=100000,Tup=60,U0=40,theta=0.375,lambda=1,alpha=2,beta=0.5,plot.ruin=TRUE){
        set.seed(100)
        # Calculate Premium according to expected value principle
        P <- (1+theta)*lambda*alpha/beta
        intervals <- Tup
        
        # Generate claim numbers for each time period and each simulation
        X <- matrix(rpois(nsim*intervals,lambda),nrow=intervals,ncol=nsim)
        
        # For each period with at least 1 claim, find the sum of amount of each claim
        npos <- sum(X>0)
        X[X>0] <- rgamma(npos,alpha*X[X>0],beta)
        
        # Compute surplus process
        X <- rbind(U0,P-X)
        X <- apply(X,2,cumsum)
        
        # Calculate ruin probabilities and find average
        ruin_count <- apply(X,2,is_ruin)
        print(mean(ruin_count))
        
        if (plot.ruin == TRUE) {
            t <- seq(0,Tup,length=intervals+1)
            plot(t,X[,1],xlim=c(0,Tup),ylim=c(min(X,0),max(X)),type="n",ylab="Insurer's captial (\u20AC)",
                 xlab="Time (months)",main="Surplus process",xaxs="i") 
            
            for(j in 1:nsim)
                if (all(X[,j] >= 0)) {
                    lines(t,X[,j], lwd=0.01)
                } else {
                    lines(t,X[,j], lwd=0.01, col="red")
                }
            abline(h=0,col="red")
        }
    }
# 5 years
simulateRuin()
# 1 year
simulateRuin(Tup=12)

## Calculating adjustment coefficient
# Create objective function
eqR <- function(r) {
    loading = 0.375
    lambda = 1
    alpha = 2
    beta = 0.5
    expec_Y = alpha/beta
    premium = (1 + loading)*lambda*expec_Y
    return(premium*r - lambda*((beta/(beta-r))^2 - 1))
}

# Create a function for finding the root of eqR for theta
R <- uniroot(eqR,lower=0.0001,upper=0.5)$root

# Calculate exponential upper bound
upper_bound <- exp(-R*40)

# Print values of R and upper bound
R
upper_bound

# Task 2 Part B: With Reinsurance
# Reinsurance Product A
# Proportional Reinsurance with 
# premium loading factor of 50% and 
# the direct insurer retains alpha = 0.7 of each claim
proportionalRuin <-
    function(retention=0.7,nsim=100000,Tup=60,U0=40,theta=0.375,lambda=1,alpha=2,beta=0.5,plot.ruin=TRUE){
        set.seed(100)
        # Calculate Premium according to expected value principle
        P <- ((1+theta)*lambda*alpha/beta) - ((1+0.5)*lambda*(1-retention)*alpha/beta)
        intervals <- Tup
        
        # Generate claim numbers for each time period and each simulation
        X <- matrix(rpois(nsim*intervals,lambda*Tup/intervals),nrow=intervals,ncol=nsim)
        
        # For each period with at least 1 claim, find the sum of amount of each claim
        npos <- sum(X>0)
        X[X>0] <- rgamma(npos,alpha*X[X>0],beta)
        X[X>0] <- X[X>0]*retention
        # Compute surplus process
        X <- rbind(U0,P-X)
        X <- apply(X,2,cumsum)
        
        # Calculate ruin probabilities and find average
        ruin_count <- apply(X,2,is_ruin)
        print(mean(ruin_count))
        
        if (plot.ruin == TRUE) {
            t <- seq(0,Tup,length=intervals+1)
            plot(t,X[,1],xlim=c(0,Tup),ylim=c(min(X,0),max(X)),type="n",ylab="Insurer's captial (\u20AC)",
                 xlab="Time (months)",main="Surplus process",xaxs="i") 
            
            for(j in 1:nsim)
                if (all(X[,j] >= 0)) {
                    lines(t,X[,j], lwd=0.01)
                } else {
                    lines(t,X[,j], lwd=0.01, col="red")
                }
            abline(h=0,col="red")
        }
    }
# 5 years
proportionalRuin()

# 1 year
proportionalRuin(Tup=12)

######### FIND RANGE OF ALPHA FOR PROPORTIONAL REINSURANCE ##############
# Find range of alpha
find_range_of_alpha <- function() {
    surplus <- function(alpha) {
        loading <- 0.375
        claim_freq <- 1
        expec_Y <- 4
        premium <- (1 + loading)*claim_freq*expec_Y
        
        reinsurer_loading <- 0.5
        reinsurer_P <- (1 + reinsurer_loading)*claim_freq*(1 - alpha)*expec_Y
        
        expec_claim_amount <- alpha* expec_Y * claim_freq
        return(premium - reinsurer_P - expec_claim_amount)
    }
    retention <- c()
    net <- c()
    x<-0
    range <- c()
    for (i in seq(from=0,to=1,by=0.0001)){
        retention[x] <- i
        net[x] <- surplus(alpha=i)
        x = x + 1
    }
    data <- cbind(retention, net)
    positive_values <- data[which(net > 0)]
    range <- c(min(positive_values), max(positive_values))
    print(range)
    plot(data, main="Range of Retention Rate",
         ylab="Net Earnings",xlab="Retention Rate")
    abline(h=0,col="red")  
}
find_range_of_alpha()

###########################################################################
########### FINDING OPTIMAL ALPHA THAT MAXIMISES R #######################
require(actuar)
find_optimal_alpha <- function() {
    # Objective Function
    eqR <- function(r,alpha){
        return (mgfgamma(r*alpha,shape=2,rate=0.5)-1-r*(6*alpha-0.5))
    }    
    
    fR <- function(x){
        uniroot(eqR,lower=0.001,upper=0.49,alpha=x)$root
    }
    # Value of R given retention of 0.7 and upper bound
    print(c(fR(0.7),exp(-40*fR(0.7))))
    
    alpha <- c()
    Adjustment_Coefficient <- c()
    a = 1
    for (i in seq(0.25038:1, by=0.001)) {
        alpha[a] = i
        Adjustment_Coefficient[a] = fR(i)
        a = a + 1
    }
    ans = cbind(alpha,Adjustment_Coefficient)
    plot(ans, main="Maximum Adjustment Coefficient",
         ylab="Lundberg Adjustment Coefficient (R)",
         xlab="Retention Rate (alpha)")
    
    ## Iterating over from retention of 0.7 to 0.999
    R <- optimise(fR,interval=c(0.25,1),maximum=TRUE)
    
    maximised_R <- R$objective
    optimal_alpha <- R$maximum
    upper_bound <- exp(-R$objective * 40)
    ans = c(optimal_alpha, maximised_R, upper_bound)
    
    return(ans)
}
find_optimal_alpha()

# Value of R given d = 6:
fR(0.7)
# Upper bound
exp(-40*fR(0.7))

# Task 2 Part B
# an Excess of Loss (EoL) reinsurance with a limit d = 6
# and the reinsurance company charges a premium loading factor
# of 50% for this EoL reinsurance.

# Find the approximated ruin probabilities within 5 years
nonProportionalRuin <-
    function(limit=6,nsim=1,loading=0.5,Tup=60,U0=40,theta=0.375,lambda=1,alpha=2,beta=0.5,plot.ruin=TRUE){
        set.seed(100)
        # Calculate Premium according to expected value principle
        # Formula for Stop Loss Premium
        integrand <- function(x) (x-limit)*(dgamma(x,shape=alpha,rate=beta))
        
        value <- integrate(integrand,lower=limit,upper=Inf)$value
        pH <- (1+loading)*lambda*value
        
        P = ((1+theta)*lambda*alpha/beta) - pH
        intervals <- Tup
        
        # Generate claim numbers for each time period and each simulation
        X <- matrix(rpois(nsim*intervals,lambda),nrow=intervals,ncol=nsim)
        
        # For each period with at least 1 claim, find the sum of amount of each claim
        for (j in 1:nsim) {
            for (i in 1:Tup) {
                number_of_claims <- X[i,j]
                if (number_of_claims > 0) {
                    total = 0
                    for (k in 1:number_of_claims) {
                        claim_amount = rgamma(1,alpha,beta)
                        if (claim_amount > limit) {
                            claim_amount = limit
                            
                        }
                        total = total + claim_amount
                    }
                    X[i,j] = total
                }
            }
        }
        
        # Compute surplus process
        X <- rbind(U0,P-X)
        
        X <- apply(X,2,cumsum)
        
        # Calculate ruin probabilities and find average
        ruin_count <- apply(X,2,is_ruin)
        print(mean(ruin_count))
        
        if (plot.ruin == TRUE) {
            t <- seq(0,Tup,length=intervals+1)
            
            plot(t,X[,1],xlim=c(0,Tup),ylim=c(min(X,0),max(X)),type="n",
                 ylab="Insurer's captial (\u20AC)",
                 xlab="Time (months)",main="Surplus process",xaxs="i") 
            
            for(j in 1:nsim)
                if (all(X[,j] >= 0)) {
                    lines(t,X[,j], lwd=0.01)
                } else {
                    lines(t,X[,j], lwd=0.01, col="red")
                }
            abline(h=0,col="red")
        }
    }
# 5 years
nonProportionalRuin()
# 1 year
nonProportionalRuin(Tup=12)

# Find range of d
find_range_of_d <- function() {
    surplus_np <- function(d) {
        theta <- 0.375
        lambda <- 1
        alpha <- 2
        beta <- 0.5
        premium <- ((1+theta)*lambda*alpha/beta)
        loading <- 0.5
        
        # Find pi_h
        integrand <- function(x) {
            (x-d)*(dgamma(x,alpha,beta))
        }
        value <- integrate(integrand,lower=d,upper=Inf)
        premium_H <- (1+loading)*lambda*value$value
        
        # Total expected claim amount
        expec_Y <- lambda*(alpha/beta)
        
        # Expected amount of claims covered by reinsurer
        expec_hY <- value$value
        return(premium - premium_H - (expec_Y - expec_hY))
    }
    limit <- c()
    net <- c()
    x<-0
    range <- c()
    for (i in seq(from=1,to=1000,by=0.001)){
        limit[x] <- i
        net[x] <- surplus_np(d=i)
        x = x + 1
    }
    total <- cbind(limit, net)
    positive_values <- total[which(net >= 0)]
    range <- c(min(positive_values), max(positive_values))
    print(range)
    plot(total, main="Range of Limit",ylab="Net Earnings",
         xlab="Limit")
    abline(h=0,col="red")  
}
find_range_of_d()

# Optimizing adjustment coefficient
# Non-Proportional
premium_H = function(d) {
    integrand <- function(x) {
        (x-d)*(dgamma(x,2,0.5))
    }
    value <- integrate(integrand,lower=d,upper=Inf)
    ans <- 1.5*value$value
    return(ans)
}
expected = function(r,d) {
    integrand2 <- function(x) {
        0.25*x*exp(-0.5*x)*exp(r*x)
    }
    value2 = integrate(integrand2,lower=0,upper=d)
    
    integrand3 <- function(x) {
        exp(r*d)*0.25*x*exp(-0.5*x)
    }
    value3 <- integrate(integrand3,lower=d,upper=Inf)
    
    return(value2$value+value3$value)
}

eqR = function(r,d) {
    return(1 + r*(5.5 - premium_H(d)) - expected(r,d))
    
}

fR <- function(x){
    uniroot(eqR,lower=0.001,upper=0.5,d=x)$root
}
R <- optimise(fR,interval=c(0,10),maximum=TRUE)

i = 1.1
b = 1
d = c()
Adjustment_Coefficient <- c()
while(i < 10) {
    d[b] = i
    Adjustment_Coefficient[b] = fR(i)
    b = b + 1
    i = i+0.01
}
total = cbind(d,Adjustment_Coefficient)
plot(total, main="Maximum Adjustment Coefficient",
     ylab="Lundberg Adjustment Coefficient (R)",
     xlab="Limit (d)")
#upper bound
upper_bound_np = exp(-40*R$objective)
upper_bound_np

# Value of R given d = 6:
fR(6)
# Upper bound
exp(-40*fR(6))



