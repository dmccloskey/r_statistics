from .r_dependencies import *

class r_base():
    def __init__(self,R_packages=None,bioconductor_packages=None,upgrade_bioconductor=False):
        self.stats = importr('stats');
        self.tools = importr('tools');
        try:
            if R_packages is None:
                R_packages = ["Amelia",
                    "mixOmics",
                    "pls",
                    "spls",
                    "caret",
                    "coin",
                    "rpart",
                    "e1071",
                    "class",
                    "cluster",
                    "randomForest",
                    "DAAG",
                    "boot",
                    "devtools"];
            if bioconductor_packages is None:
                bioconductor_packages = ["Biobase",
                    "LMGene",
                    "pcaMethods",
                    "ropls",]
            self.import_RPackages(R_packages,bioconductor_packages,upgrade_bioconductor);
        except Exception as e:
            print(e);

    def import_RPackages(self,
            R_packages_I=[],
            bioconductor_packages_I=[],
            upgrade_bioconductor_I=False
            ):
        '''load required R packages
        NOTE: must be run as administrator if packages need to be installed!
        INPUT:
        R_packages_I = list of R packages
        bioconductor_packages_I = list of bioconductor packages
        upgrade_bioconductor_I = boolean
        '''

        #lib_loc = '"C:/Users/Douglas/Documents/Douglas/R/win-library/3.0"';
        #r_statement = ('library("Amelia",lib.loc = %s)' % lib_loc);

        for package in R_packages_I:
            self.library_R(package);
            
        #Upgrade bioconductor
        if upgrade_bioconductor_I:
            self.upgrade_bioconductor();

        for package in bioconductor_packages_I:
            self.library_bioconductor(package);

    def upgrade_bioconductor(self,source_I='https://bioconductor.org/biocLite.R'):
        '''upgrade/install bioconductor'''
        
        r_statement = ('source("%s")'%(source_I))
        ans = robjects.r(r_statement);
        try:
            r_statement = ('biocLite("BiocUpgrade",ask=FALSE)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('biocLite(ask=FALSE)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);

    def library_R(self,library_I):
        '''Load R library
        INPUT:
        library_I
        '''
        try:
            r_statement = ('library("%s")' %(library_I));
            ans = robjects.r(r_statement);
            r_statement = ('require(%s)'%(library_I));
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("%s",dependencies=TRUE)'%(library_I));
                ans = robjects.r(r_statement);
                r_statement = ('library("%s")'%(library_I));
                ans = robjects.r(r_statement);
                r_statement = ('require(%s)'%(library_I));
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);

    def library_bioconductor(self,library_I,
            source_I='https://bioconductor.org/biocLite.R'):
        '''Load bioconductor library
        INPUT:
        library_I
        source_I = bioconductor source url
        '''
        try:
            r_statement = ('library("%s")' %(library_I));
            ans = robjects.r(r_statement);
            r_statement = ('require(%s)'%(library_I));
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('source("%s")'%(source_I));
                ans = robjects.r(r_statement);
                r_statement = ('biocLite("%s",ask=FALSE)'%(library_I));
                ans = robjects.r(r_statement);
                r_statement = ('library("%s")'%(library_I));
                ans = robjects.r(r_statement);
                r_statement = ('require(%s)'%(library_I));
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);

    def make_matrixFromList(self,list_I,nrows_I,ncolumns_I,matrix_O):
        '''R commands to make a matrix from a list of data
        INPUT:
        list_I = list of data
        nrows_I = int, # of rows
        ncolumns_I = int, # of columns
        OUTPUT:
        matrix_O = name of the R matrix in the workspace'''
        
        try:
            # convert lists to R matrix
            matrix_r = '';
            for c in list_I:
                matrix_r = (matrix_r + ',' + str(c));
            matrix_r = matrix_r[1:];
            r_statement = ('matrix_tmp = c(%s)' % matrix_r);
            ans = robjects.r(r_statement);
            r_statement = ('%s = matrix(matrix_tmp, nrow = %s, ncol = %s, byrow = TRUE)' %(matrix_O,nrows_I,ncolumns_I)); 
            ans = robjects.r(r_statement);
            # cleanup
            variables_delete = ['matrix_tmp'];
            self.remove_workspaceVariables(variables_I=variables_delete);
        except Exception as e:
            print(e);
            exit(-1);
    def make_factorsFromList(self,list_I,factors_O):
        '''R commands to make a factors from list
        INPUT:
        list_I = list of data
        OUTPUT:
        factors_O = name of the R factor in the workspace
        '''
        try:
            factor_r = '';
            for c in list_I:
                factor_r = (factor_r + ',' + '"' + str(c) + '"');
            factor_r = factor_r[1:];
            r_statement = ('factors = c(%s)' % factor_r);
            ans = robjects.r(r_statement);
            r_statement = ('%s = factor(factors)' % factors_O); 
            ans = robjects.r(r_statement);
            # cleanup
            variables_delete = ['factors'];
            self.remove_workspaceVariables(variables_I=variables_delete);
        except Exception as e:
            print(e);
            exit(-1);
    def make_dummyMatrixFromFactors(self,factors_I,dummy_O):
        '''R commands to make a dummy factor matrix for use with plsda
        INPUT:
        factors_I = name of the R factor variable in the workspace
        OUTPUT:
        dummy_O = name of the R dummy factor matrix in the workspace
        '''
        try:
            # make the dummy response for PLS-DA
            r_statement = ('%s <- I(model.matrix(~y-1, data.frame(y = %s)))' %(dummy_O,factors_I)); 
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def transpose_matrix(self,matrix_I,matrix_O):
        '''R commands to transpose a matrix in the workspace
        INPUT:
        matrix_I = name of the R matrix in the workspace
        OUTPUT:
        matrix_O = name of the transposed R matrix in the workspace
        '''
        
        try:
            r_statement = ('%s = t(%s)' %(matrix_I,matrix_O)); 
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def clear_workspace(self,):
        '''Clear the R workspace'''
        try:
            # clear the workspace
            r_statement = ('rm(list = ls())');
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def remove_workspaceVariables(self,variables_I):
        '''R commands to remove specific workspace variables
        INPUT:
        variables_I = [] of strings, list of workspace variables to remove
        '''
        try:
            #concatenate all variables
            variables_str = ','.join(variables_I);
            #remove the variables
            r_statement = ('rm(%s)' %(variables_str));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_pValueCorrected(self,pvalue_I,pvalue_O,method_I = "bonferroni"):
        '''calculate the corrected p-value
        INPUT:
        pvalue_I = float, uncorrected p-value
                   OR string, name of the r-workspace variable
        method_I = string, method name
        pvalue_O = string, name of the r-workspace variable
        OUTPUT:
        pvalue_corrected_O = float, corrected p-value
        DESCRIPTION:
        https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

        p.adjust(p, method = p.adjust.methods, n = length(p))

        p.adjust.methods
        # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
        #   "fdr", "none")
        '''
        pvalue_corrected_O = None;
        try:
            if type(pvalue_I)==type(''):
                r_statement = ('%s = p.adjust(%s, method = "%s"' %(pvalue_O,pvalue_I,method_I)); 
            else:
                r_statement = ('%s = p.adjust(as.numeric(%s), method = "%s")' %(pvalue_O,pvalue_I,method_I)); 
            ans = robjects.r(r_statement);
            pvalue_corrected_O = ans[0]; 
        except Exception as e:
            print(e);
            exit(-1);
        return pvalue_corrected_O;
    def make_vectorFromList(self,list_I,vector_O):
        '''
        make an R vector from a python list
        supports conversion of the following types:
            floats
            integers
            strings
            booleans
        INPUT:
        list_I = list of floats or strings or booleans
        OUTPUT:
        vector_O = string, name of the R workspace variable

        TODO:
        add in additional/better tests for numpy types
        '''
        try:
            if type(list_I[0])==type(1.0) or type(list_I[0])==type(1)\
                or str(type(list_I[0]))=="<class 'numpy.int64'>"\
                or str(type(list_I[0]))=="<class 'numpy.int8'>"\
                or str(type(list_I[0]))=="<class 'numpy.int16'>"\
                or str(type(list_I[0]))=="<class 'numpy.int32'>"\
                or str(type(list_I[0]))=="<class 'numpy.float16'>"\
                or str(type(list_I[0]))=="<class 'numpy.float32'>"\
                or str(type(list_I[0]))=="<class 'numpy.float64'>":
                list_str = '';
                for c in list_I:
                    list_str = (list_str + ',' + str(c));
                list_str = list_str[1:];
                r_statement = ('%s = c(%s)' % (vector_O,list_str));
                ans = robjects.r(r_statement);
            elif type(list_I[0])==type('') or str(type(list_I[0]))=="<class 'numpy.str_'>":
                list_str = '';
                for c in list_I:
                    list_str = (list_str + ',' + '"' + c + '"');
                list_str = list_str[1:];
                r_statement = ('%s = c(%s)' % (vector_O,list_str));
                ans = robjects.r(r_statement);
            elif type(list_I[0])==type(True):
                list_str = '';
                for c in list_I:
                    if c:
                        list_str = (list_str + ',' + "TRUE");
                    else:
                        list_str = (list_str + ',' + "FAlSE");
                list_str = list_str[1:];
                r_statement = ('%s = c(%s)' % (vector_O,list_str));
                ans = robjects.r(r_statement);
            else:
                print('list type not recognized.');
        except Exception as e:
            print(e);
            exit(-1);

