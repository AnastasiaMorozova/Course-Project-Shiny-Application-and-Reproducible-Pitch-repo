#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Endogeneity: Simulation study"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3('Parameters to choose for the Monte Carlo study'),
      
      numericInput('endogeneity', 'Numeric input: degree of endogeneity, defined by correlation between error terms', 0.5, min = 0.1, max = 0.8, step = 0.1),
      sliderInput("number",
                   "Number of instrumental variables, K:",
                   min = 1,
                   max = 20,
                   value = 1),
      sliderInput("rsq",
                  "Theoretical R squared for the first stage regression:",
                  min = 0.1,
                  max = 0.9,
                  value = 0.8),
      checkboxGroupInput("choise", "Checkbox",
                         c("OLS estimator" = "1",
                           "2SLS estimator"= "2"),selected="1")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Theory",
                 h5('This application reproduces the simulation study illustrating the properties of the ordinary least squares and the two-stage least squares estimators in the presence of the endogeneity problem.
'),
                 h2('Endogeneity'),
                 h5('In econometrics, endogeneity broadly refers to situations in which an explanatory variable is correlated with the error term. The reasons for appearance of endogeneity problem are manifold, it can appear because of simultaneity of the dependent and explanatory variables, when an unobserved or omitted variable is confounding both independent and dependent variables, or when independent variables are measured with error.'), 

                 h5('If the independent variable is correlated with the error term in a regression model then the estimate of the regression coefficient in an ordinary least squares (OLS) regression is biased. For correcting the bias applied researchers often use instrumental variable regression and two-stage least squares (2SLS) estimation in particular.
                    If a valid instrument (i.e. correlation between the instrument and the structural error term is equal to zero) is available, consistent estimates of the structural effect may still be obtained. An instrument is a variable that does not itself belong in the explanatory equation but is correlated with the endogenous explanatory variables, conditional on the value of other covariates.'),
                 h2('Model'),
                 withMathJax(helpText("Structual equation with one endogenous explanatory variable X: $$Y_i=\\beta X_i + \\varepsilon_i$$")),
                 div("$$cov(X_i,\\varepsilon_i) \\neq 0$$"),
                 div("Structural effect of interest: $$\\beta$$ "),
                 div("First-stage regression with K instrumental variables Z: $$	X_i=Z_i'\\pi + u_i$$"),
                 div("$$Z_i=(Z_{i1},Z_{i2},\\dots, Z_{iK})', cov(Z_i,\\varepsilon_i) = 0, cov(X_i,Z_i) \\neq 0$$"),
                 div("Martix notation (N - number of observations): $$Y=(Y_1,\\dots,Y_N)',X=(X_1,\\dots,X_N)', \\varepsilon=(\\varepsilon_1,\\dots,\\varepsilon_N)', Z=(Z_1',\\dots,Z_N')'$$")
        ),
        tabPanel("Properties of the OLS",
                 h4('Formula:'),
                 div("$$\\widehat{\\beta}_{OLS}=(X'X)^{-1}X'Y=\\beta+(X'X)^{-1}X'\\varepsilon$$"),
                 h4('Properties of the estimator under endogeneity:'),
                 h5('- can be inconsistent estimator of the true structural effect'),
                 div("- biased estimator of the true structural effect, i.e. $$ E[\\widehat{\\beta}_{OLS}] \\neq \\beta$$"),
                 h4('*Note for the output:'),                 
                 h6('One can see from the simulation study that OLS estimator is biased under the edogeneity by comparing the mean estimate with the true value of the structural parameter (red and black line), as well as checking the mean bias and MSE of the estimator.
                    The median absolute deviation (MAD) describes the variability of the OLS estimates.')
        ),
        tabPanel("Properties of the 2SLS",
                 h3('Formula:'),
                 div("$$\\widehat{\\beta}_{2SLS}=(X'Z(Z'Z)^{-1}Z'X)^{-1}X'Z(Z'Z)^{-1}Z'Y=\\beta+(X'Z(Z'Z)^{-1}Z'X)^{-1}X'Z(Z'Z)^{-1}Z'\\varepsilon$$"),
                 h4('Properties of the estimator:'),
                 h5('- consistent estimator of the true structural effect'),
                 div("- asymptotically unbiased estimator of the true structural effect"),
                 h4('Properties of the 2SLS estimator in the special cases:'),
                 h5('- If the instruments are only weakly correlated with the endogenous variable X (weak instruments problem), the 2SLS estimator becomes biased and inconsistent.'),
                 h5('- Increasing the number of instrumental variables, K, leads to a higher asymptotic efficiency of the 2SLS estimator (i.e. smaller asymptotic variance). However, it may lead to an increased bias of the estimator in the finite samples.'),
                 h4('*Note for the output:'),                  
                 h6('One can see from the simulation study that 2SLS estimator is biased or unbiased for a particular design by comparing the mean estimate with the true value of the structural parameter (red and black line), as well as checking the mean bias and MSE of the estimator.
                    The median absolute deviation (MAD) describes the variability of the 2SLS estimates.')
        ),
        tabPanel("Simulation design",
                 h2('Model'),
                 div("Structual equation with one endogenous explanatory variable X: $$Y_i=\\beta X_i + \\varepsilon_i= \\beta Z_i'\\pi + v_i, v_i= \\beta u_i+\\varepsilon_i$$"),
                 div("First-stage regression with K instrumental variables Z: $$	X_i=Z_i'\\pi + u_i$$"),
                 h2('Initial parameters:'),
                 h5('- Number of replications: 1000. Number of observations, N: 100. Number of instruments, K: chosen by a user from 1 to 20.'),
                 div("- Intruments: $$Z_i \\sim N(0,I_K)$$"),
                 div("- Error term of reduced form equations: $$ \\begin{pmatrix} v_i \\ u_i\\end{pmatrix} ' \\sim N \\bigg( 0, \\begin{bmatrix}  1 & w_{12} ;\\  w_{12} &  1 \\end{bmatrix} \\bigg), w_{12}=0.5*\\beta$$"),
                 div("- Degree of endogeneity, $$corr(\\varepsilon_i,u_i)$$ defined by the user in the range between 0.1 and 0.8."),
                 div("- Value of the structural parameter beta is set such that $$Var(\\varepsilon_i)=1$$"),
                 div("- Strength of the insumental variables set is defined by the first-stage theoretical R-squared: $$R_f^2$$ defined by the user in the range between 0.1 (weak instruments) and 0.9 (strong instruments)."),
                 div("- First stage parameter for every instrument 1,..., K: $$\\pi_k=\\sqrt{\\frac{R_f^2}{K(1-R_f^2)}}$$")
        ),
        tabPanel("Output",
                 h3('Distribution plot'),
                 plotOutput("distPlot"),
                 h5('Black line - mean value of the estimates. The true value of the structural parameter in the model is given by the red line and equal to'),
                 verbatimTextOutput("betatrue"),
                 h4('Mean bias of the estimator:'),
                 verbatimTextOutput("Bias"),
                 h4('MSE of the estimator:'),
                 verbatimTextOutput("MSEest"),
                 h4('Median absolute deviation of the estimator:'),
                 verbatimTextOutput("MADest")
                 )
    )
    )
  )
))
