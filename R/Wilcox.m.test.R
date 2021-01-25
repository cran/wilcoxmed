#' 1-Sample  Wilcoxon Sign Rank Hypothesis Test for Medians
#'
#' @name Wilcox.m.test
#' @description This function allows the user to conduct the 1-Sample  Wilcoxon Sign Rank Hypothesis
#' Test for Medians using the probability values from the exact
#' distribution of W.
#'
#' @usage Wilcox.m.test(dat, m_h0, alpha = 0.05,
#' alternative=c('greater', 'lesser', 'noteq'), normal_approx=FALSE)
#' @param dat data vector relating to the sample the user is
#' performing the hypothesis test for
#' @param m_h0 The value of the median as specified by the null hypothesis H_0
#' @param alpha The significance level of the hypothesis test (default = 0.05)
#' @param alternative The sign of the alternative hypothesis.
#' e.g 'greater' - H_1:m>m_h0 , 'lesser' - H_1:m<m_h0, 'noteq' - H_1:m!=m_h0
#' @param normal_approx Should the normal approximation test be applied? (default = FALSE)
#' @return Prints out the results of the tests, and returns 3 values- test statistic,
#' p-value, and the significance level of the test, alpha
#' @references Peter J. Bickel and Kjell A. Doksum (1973). \emph{Mathematical Statistics:
#' Basic Ideas and Selected Topics}. Prentice Hall.
#' @examples
#' ##Given some data: 3, 4, 7, 10, 4, 12, 1, 9, 2, 15
#' ##If we want to test the hypotheses H_0: m=5 against H_1: m>5
#' ##without using normal approximation:
#' vec = c(3, 4, 7, 10, 4, 12, 1, 9, 2, 15)
#' res = Wilcox.m.test(dat = vec, m_h0 = 5,
#' alternative = 'greater', normal_approx = FALSE)
#'
#' ##If we want to apply the normal approximation(Z-test), with the same hypotheses:
#' res = Wilcox.m.test(dat = vec, m_h0 = 5,
#' alternative = 'greater', normal_approx = TRUE)
#' @details This hypothesis test allows breaking of ties, and the number of
#' ties broken is also reflected in the printed results.
#' @seealso \code{\link{wilcox.test}} for the same tests applied to 2 sample problems
#' but is not able to break ties
#' @export
#'
Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)

Wilcox.m.test <- function(dat, m_h0 ,alpha = 0.05,
                          alternative=c('greater', 'lesser', 'noteq'),
                          normal_approx=FALSE){

  console_length = getOption('width')
  n <- length(dat)
  Ri <- sign(dat-m_h0)*rank(abs(dat-m_h0))
  ties <-length(abs(dat-m_h0)) - length(unique(rank(abs(dat-m_h0))))
  W_obs <- sum(Ri[which(Ri>0)])
  if(!normal_approx){
    if(alternative == 'greater'){
      p_value <- W_stat(n, W_obs, 'geq')
      cat('\n', rep('', floor((console_length-42)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m >', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic W =', W_obs,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }
    else if(alternative == 'lesser'){
      p_value <- W_stat(n, W_obs, 'leq')
      cat('\n', rep('', floor((console_length-42)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m <', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic W =', W_obs,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }
    else if(alternative == 'noteq'){
      p_value <- min(2*W_stat(n, W_obs, 'geq'),1)
      cat('\n', rep('', floor((console_length-42)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m !=', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic W =', W_obs,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }

  }
  else if(normal_approx){
    EW <- n*(n+1)/4
    SDW <- sqrt(n*(n+1)*(2*n+1)/24)
    Z_stat  <- (W_obs - EW)/SDW
    if(alternative == 'greater'){
      p_value <- stats::pnorm(Z_stat, lower.tail = FALSE)
      cat('\n', rep('', floor((console_length-60)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test - Normal Approximation', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m >', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic Z =', Z_stat ,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }
    else if(alternative == 'lesser'){
      p_value <- stats::pnorm(Z_stat, lower.tail = TRUE)
      cat('\n', rep('', floor((console_length-60)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test - Normal Approximation', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m <', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic Z =', Z_stat ,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }
    else if(alternative == 'noteq'){
      p_value <- min(2*stats::pnorm(Z_stat, lower.tail = TRUE),1)
      cat('\n', rep('', floor((console_length-60)/2)),
          ' 1-Sample Median Wilcoxon Sign-Rank Test - Normal Approximation', '\n','\n',
          rep('', floor((console_length-10)/2)), 'Ties Broken =', ties, '\n','\n',
          rep('', floor((console_length-60)/2)),
          'Null hypothesis H_0: m =', m_h0, rep('', 5),
          'Alternative Hypothesis H_1: m !=', m_h0,'\n','\n',
          rep('', floor((console_length-60)/2)),
          'Test Statistic Z =', Z_stat ,rep('', 5),
          'p-value = ', p_value , rep('', 5), 'alpha =', alpha,'\n','\n')
      if(p_value<=alpha){
        cat(rep('', floor((console_length-30)/2)),
            'Test Result: Reject H_0' , '\n','\n')
      }
      else{cat(rep('', floor((console_length-30)/2)),
               'Test Result: Do not Reject H_0' , '\n','\n')}
    }
  }
  if(!normal_approx){
    invisible((list('p-value' = p_value, 'Test Statistic W' =  W_obs, 'alpha' = alpha)))
  }
  else{
    invisible((list('p-value' = p_value, 'Test Statistic Z' =  Z_stat, 'alpha' = alpha)))
  }
}


