from .r_dependencies import *
from .r_base import r_base
class r_decisionTree(r_base):
    def calculate_decisionTree(self):
        '''calculate a decision tree
        requires the rpart package
        e.g. cart_model <- rpart(y  x1 + x2, data=as.data.frame(cbind(y,x1,x2)), method="class")
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