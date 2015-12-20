from .r_dependencies import *
from .r_base import r_base
class r_svm(r_base):
    def calculate_svm(self):
        '''Support Vector Machine
        requires package e1071
        
        Documentation:
        tune.svm(train.x, train.y = NULL, data = list(), validation.x = NULL,
            validation.y = NULL, ranges = NULL, predict.func = predict,
            tunecontrol = tune.control(), ...)
        tune.control(random = FALSE, nrepeat = 1, repeat.aggregate = min,
            sampling = c("cross", "fix", "bootstrap"), sampling.aggregate = mean,
            sampling.dispersion = sd,
            cross = 10, fix = 2/3, nboot = 10, boot.size = 9/10, best.model = TRUE,
            performances = TRUE, error.fun = NULL)
        best.tune(...)
        svm(x, y = NULL, scale = TRUE, type = NULL, kernel = "radial",
            degree = 3, gamma = if (is.vector(x)) 1 else 1 / ncol(x),
            coef0 = 0, cost = 1, nu = 0.5,
            class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
            shrinking = TRUE, cross = 0, probability = FALSE, fitted = TRUE, seed = 1L)
        
        
        '''
        # TODO
        # Call to R
        try:
            # format the data into R objects
            # call tune
            r_statement = ('concentrations_m = matrix(concentrations, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted))); 
            ans = robjects.r(r_statement);
            return
        except Exception as e:
            print(e);
            exit(-1);