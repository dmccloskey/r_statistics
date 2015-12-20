from .r_dependencies import *
from .r_base import r_base
from math import sqrt
class r_svd(r_base):
    def calculate_pcaMethods_svd(self):
        '''calculate svd using pcaMethods'''
        pass;
    def calculate_pcaMethods_robustSvd(self):
        '''calculate robustSvd using pcaMethods
        The robust SVD of the matrix is x = u d v'.
            d A vector containing the singular values of x.
            u A matrix whose columns are the left singular vectors of x.
            v A matrix whose columns are the right singular vectors of x
        
        '''
        pass;
    def detect_outliers_svd(self):
        '''detect outliers using pcaMethods

        ## Load a complete sample metabolite data set and mean center the data
        data(metaboliteDataComplete)
        mdc <- prep(metaboliteDataComplete, center=TRUE, scale="none")
        ## Now create 5% of outliers.
        cond <- runif(length(mdc)) < 0.05;
        mdcOut <- mdc
        mdcOut[cond] <- 10
        ## Now we do a conventional SVD and a robustSvd on both, the original and the
        ## data with outliers.
        resSvd <- svd(mdc)
        resSvdOut <- svd(mdcOut)
        resRobSvd <- robustSvd(mdc)
        resRobSvdOut <- robustSvd(mdcOut)
        ## Now we plot the results for the original data against those with outliers
        ## We can see that robustSvd is hardly affected by the outliers.
        plot(resSvd$v[,1], resSvdOut$v[,1])
        plot(resRobSvd$v[,1], resRobSvdOut$v[,1])
        '''
        pass;