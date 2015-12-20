from .r_dependencies import *
from .r_base import r_base
class r_correlation(r_base):
    def calculate_r_correlation(self,
            data_var_I,
            data_var_O,
            y="NULL",
            use="everything",
            method = "spearman"
            ):
        '''Calculate the correlation matrix
        INPUT:
        data_var_I = name of the R matrix
        y = vector or matrix with compatible dimensions to x (default = NULL, i.e., y=x)
        use =
        method = "pearson","kendall","spearman"
        OUTPUT:
        data_O = correlation matrix
        data_var_O = name of the R object containing the correlation matrix
        '''
        try:
            r_statement = ('%s <- cor(%s, y=%s, method = "%s")' %(
                data_var_O,ata_var_I,mvr_model_I,method))
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans); #dim 1 = features, dim 2 = comps
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);