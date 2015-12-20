from .r_dependencies import *
from .r_base import r_base
class r_geneticAlgorithm(r_base):
    def calculate_geneticAlgorithm(self):
        '''genetic algorithm
        requires the GA or genalg package
        e.g. http://www.r-bloggers.com/genetic-algorithms-a-simple-r-example/
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