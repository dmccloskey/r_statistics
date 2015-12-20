from .r_dependencies import *
from .r_base import r_base
from math import sqrt
class r_outlier(r_base):

    def detect_outliers_pca(self):
        '''detect outliers using pcaMethods
        alternatively:
        1. pca
        2. robustPca
        3. plot loadings of pca vs. loadings of robustPca'''