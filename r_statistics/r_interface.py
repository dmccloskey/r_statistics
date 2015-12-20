from .r_base import r_base
# r modules
from .r_anova import r_anova
from .r_dataNormalization import r_dataNormalization
from .r_missingValues import r_missingValues
from .r_pca import r_pca
from .r_pls import r_pls
from .r_statistics import r_statistics

class r_interface(
        r_anova,
        r_dataNormalization,
        r_missingValues,
        r_pca,
        r_pls,
        r_statistics):
    '''conveniency class that wraps all r class'''
    pass;