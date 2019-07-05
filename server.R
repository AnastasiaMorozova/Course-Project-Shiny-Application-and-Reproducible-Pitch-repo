#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
#install.packages(mvtnorm)
library(mvtnorm) # for random numbers

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    # if we vary the theoretical R squared and then define reduced form parameters
    R_sq<-input$rsq;  # theoretical first stage regression R squared 
    
    set.seed(2019)
    Rep<-1000;   
    N<-100;           
    K<-1;             # number of endogenous variables      K=1
    q<-input$number;
    ## parameters of the model
    
    sig2_z=1        # variance of the instrument
    VC_z<-diag(q)*sig2_z;  # VC matrix of the instruments
    sig2_u=1        # variance of red.form error term
    sig2_eps=1      # variance of the struct. eq. error term
    cor_eps_u=input$endogeneity  # correlation between epsilon and u ro={-0.9,-0.8,-0.5,-0.25,-0.1,-0.05,0,0.1,0.5,0.8,0.9}
    cov_eps_u=sqrt(sig2_eps)*sqrt(sig2_u)*cor_eps_u  # covariance between epsilon and u
    # according to the formulas in Hahn and Hausman (2002), Hahn et al. (2004)
    w12<-cov_eps_u*(-1)  
    b0<-2*w12
    sig2_v<-(b0^2)*sig2_u+sig2_eps+2*cov_eps_u*b0
    Sigma_reduced_matrix<-array(c(sig2_v=1,w12,w12,sig2_u=1),dim=c(2,2)); 
    
    # predefine the reduced form parameters
    pi0=NULL
    for(i in 1:q) {
      # i-th element of `pi` is according to the Rsquared
      pi0[i] <- sqrt(R_sq/(q*(1-R_sq)));
     
    }
  
    ################################################################################
    ##                              SPACE HOLDERS                                 ##
    ################################################################################
    # For conventional estimators
    bols= matrix(NA,nrow=Rep,ncol=1);
    sols = matrix(NA,nrow=Rep,ncol=3);
    b2sls = matrix(NA,nrow=Rep,ncol=1);
    s2sls = matrix(NA,nrow=Rep,ncol=3);
   
    # F-test and Wald test
    FandWald=matrix(NA,nrow=Rep,ncol=5);
    colnames(FandWald)=c("F-stat","RuleF","Wald-stat","RuleWald","Check F")
   
    
    ################################################################################
    ##                          SIMULATION LOOP                                   ##
    ################################################################################
    for (r in 1:Rep){;
      #                             simulation of the data             
      
      #######  error generation options
      { # if generated errors from the reduced forms using multivariate normal Hahn et al. (2004)
        # 
        err_r<-rmvnorm(N, c(0,0), Sigma_reduced_matrix)
        v1_r<-err_r[,1];
        u_r<-err_r[,2];
        eps_r<-v1_r-b0*u_r;
      }
      
      
      ZOrig_r<-rmvnorm(N, rep(0, q), VC_z);
      # IV regression model
      xOrig_r<-ZOrig_r%*%pi0+u_r;
      yOrig_r<-xOrig_r%*% b0+eps_r;
      
      # STATISTICS INFEASIBLE
      # Wald statistics
      Wald_stat<-(t(pi0)%*%solve(as.numeric(sig2_u)*solve(t(ZOrig_r)%*%ZOrig_r))%*%pi0)
      # F test statistics
      F_stat<-Wald_stat/q
      
      #                           estimation part and testing
      
      # estimation preparation, centering         REWRITE
      Z_r <- ZOrig_r - array(rep(colMeans(ZOrig_r),each = N), dim=c(N,q));  # check how to create vector of ones
      # hist(Z_r[,1])
      x_r <- xOrig_r - mean(xOrig_r);
      y_r <- yOrig_r - mean(yOrig_r);
      # First stage F test 
      zz <- t(Z_r)%*%Z_r
      zx <- t(Z_r)%*%x_r
      pi_est <- solve(zz,zx)
      
      first_stage_resid_r <- x_r-Z_r%*%pi_est  # first stage residual 
      s2_resid_r<-(t(first_stage_resid_r) %*% first_stage_resid_r)/(N-q) # estimated variance of the error term
      # first stage Wald and F statistics 
      Wald1_r<-(t(pi_est)%*%solve(as.numeric(s2_resid_r)*solve(t(Z_r)%*%Z_r))%*%pi_est)
      F1_r<-Wald1_r/q
      FandWald[r,1]<-F1_r
      FandWald[r,2]<-F_stat
      FandWald[r,3]<-Wald1_r
      FandWald[r,4]<-Wald_stat
      FandWald[r,5]<-as.numeric(summary(lm(x_r~Z_r))$fstatistic[1])
     
     
      #                              estimation OLS  
      
      xx <- t(x_r)%*%x_r
      xy <- t(x_r)%*%y_r
      beta_ols_r <- solve(xx,xy)
      bols[r,1]<-beta_ols_r;           # OLS estimate
      resid_ols_r <- y_r-x_r%*%beta_ols_r  # first stage residual 
      s2_res_ols_r<-(t(resid_ols_r) %*% resid_ols_r)/(N-K) # estimated variance of the error term
      sols[r,1]<-beta_ols_r-b0;       # for mean and median bias
      sols[r,2]<-(beta_ols_r-b0)^2;   # for MSE and RMSE
      sols[r,3]<-abs(beta_ols_r-b0);  # for median absolute deviation
  
      #                         estimation 2SLS with all instruments
      Pz_r<-Z_r%*%solve(t(Z_r)%*%Z_r)%*%t(Z_r)
      x_2sls_fit<-Pz_r%*%x_r
      beta_2sls_r<-solve(t(x_r)%*%Pz_r%*%x_r)%*%t(x_r)%*%Pz_r%*%y_r
      b2sls[r,1]<-beta_2sls_r;           # OLS estimate
      y_2sls_fit<-x_r%*%beta_2sls_r;
      resid_2sls_r <- y_r - x_r%*%beta_2sls_r;
      VC1 <- ((t(resid_2sls_r)%*%resid_2sls_r)/(N-K))%*%solve(t(x_r)%*%Pz_r%*%x_r);
      
      s2sls[r,1]<-beta_2sls_r-b0;       # for mean and median bias
      s2sls[r,2]<-(beta_2sls_r-b0)^2;   # for MSE and RMSE
      s2sls[r,3]<-abs(beta_2sls_r-b0);  # for median absolute deviation
      
      
     }  # end of the simulation loop
    
    # Bias, mean Bias, median Bias
    mean<-c(mean(bols),mean(b2sls))
    meanBias<-c(mean(sols[,1]),mean(s2sls[,1]))
    medianBias<-c(median(sols[,1]),median(s2sls[,1]))
    
    
    # Median absolute deviation
  MAD1<-c(median(sols[,3]),median(s2sls[,3]))
  MAD2<-c(median(abs(bols-mean(bols))), median(abs(b2sls-mean(b2sls))))
    # MSE, RMSE
    MSE<-c(mean(sols[,2]),mean(s2sls[,2]))
    RMSE<-sqrt(MSE);
    
    # generate bins based on input$bins from ui.R
  
    
    # draw the histogram with the specified number of bins
    x=bols
    if (input$choise==1){
      x=bols
      col=1
    } else {
      x=b2sls
      col=2
    }
    
    hist(x, col = 'lightgreen', border = 'white',main='Histogram of the chosen estimator', xlim=c(min(1.25*min(x),b0),max(1.25*max(x),-b0)))
    abline(v = mean(x), lwd = 2)
    abline(v = b0, col = "magenta", lwd = 4)
    
    output$betatrue<-renderPrint(b0)
    output$Bias<-renderPrint(meanBias[col])
    output$MSEest<-renderPrint(MSE[col])
    output$MADest<-renderPrint(MAD1[col])
  })
  
})
