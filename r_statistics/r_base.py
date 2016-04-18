from .r_dependencies import *

class r_base():
    def __init__(self):
        self.stats = importr('stats');
        self.tools = importr('tools');
        try:
            self._import_RPackages();
        except Exception as e:
            print(e);

    def _import_RPackages(self,
            R_packages_I=[],
            bioconductor_packages_I=[]
            ):
        '''load required R packages
        NOTE: must be run as administrator if packages need to be installed!
        INPUT:
        R_packages_I = list of R packages
        bioconductor_packages_I = list of bioconductor packages
        '''

        #lib_loc = '"C:/Users/Douglas/Documents/Douglas/R/win-library/3.0"';
        #r_statement = ('library("Amelia",lib.loc = %s)' % lib_loc);

        #Amelia (missing value imputation)
        try:
            r_statement = ('library("Amelia")');
            ans = robjects.r(r_statement);
            r_statement = ('require(Amelia)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("Amelia",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("Amelia")');
                ans = robjects.r(r_statement);
                r_statement = ('require(Amelia)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #Upgrade bioconductor
        r_statement = ('source("https://bioconductor.org/biocLite.R")');
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
        #Biobase
        try:
            r_statement = ('library("Biobase")');
            ans = robjects.r(r_statement);
            r_statement = ('require(Biobase)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('source("https://bioconductor.org/biocLite.R")');
                ans = robjects.r(r_statement);
                r_statement = ('biocLite("Biobase",ask=FALSE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("Biobase")');
                ans = robjects.r(r_statement);
                r_statement = ('require(Biobase)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #mixOmics (splsda, plsda, pca, clustering)
        try:
            r_statement = ('library("mixOmics")');
            ans = robjects.r(r_statement);
            r_statement = ('require(mixOmics)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("mixOmics",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("mixOmics")');
                ans = robjects.r(r_statement);
                r_statement = ('require(mixOmics)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #pls (pls, plsda, oplsda)
        try:
            r_statement = ('library("pls")');
            ans = robjects.r(r_statement);
            r_statement = ('require(pls)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("pls",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("pls")');
                ans = robjects.r(r_statement);
                r_statement = ('require(pls)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #spls (spls, splsda)
        try:
            r_statement = ('library("spls")');
            ans = robjects.r(r_statement);
            r_statement = ('require(spls)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("spls",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("spls")');
                ans = robjects.r(r_statement);
                r_statement = ('require(spls)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #caret (utilities for regression packages including pls and spls, and boosting methods)
        try:
            r_statement = ('library("caret")');
            ans = robjects.r(r_statement);
            r_statement = ('require(caret)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("caret",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("caret")');
                ans = robjects.r(r_statement);
                r_statement = ('require(caret)');
                ans = robjects.r(r_statement);   
            except Exception as e:
                print(e);
        ##RVAideMemoire (utilities for regression packages including pls and spls)
        #try:
        #    r_statement = ('library("RVAideMemoire")');
        #    ans = robjects.r(r_statement);
        #    r_statement = ('require(RVAideMemoire)');
        #    ans = robjects.r(r_statement);
        #except:
        #    r_statement = ('install.packages("RVAideMemoire",dependencies=TRUE)');
        #    ans = robjects.r(r_statement);
        #    r_statement = ('library("RVAideMemoire")');
        #    ans = robjects.r(r_statement);
        #    r_statement = ('require(RVAideMemoire)');
        #    ans = robjects.r(r_statement);
        #ropls (pls,opls,oplsda,plsda, pca, clustering)
        try:
            r_statement = ('library("ropls")');
            ans = robjects.r(r_statement);
            r_statement = ('require(ropls)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('source("https://bioconductor.org/biocLite.R")');
                ans = robjects.r(r_statement);
                r_statement = ('biocLite("ropls",ask=FALSE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("ropls")');
                ans = robjects.r(r_statement);
                r_statement = ('require(ropls)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #pcaMethods (missing value and pca analysis)
        try:
            r_statement = ('library("pcaMethods")');
            ans = robjects.r(r_statement);
            r_statement = ('require(pcaMethods)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('source("https://bioconductor.org/biocLite.R")');
                ans = robjects.r(r_statement);
                r_statement = ('biocLite("pcaMethods",ask=FALSE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("pcaMethods")');
                ans = robjects.r(r_statement);
                r_statement = ('require(pcaMethods)');
                ans = robjects.r(r_statement);    
            except Exception as e:
                print(e);
        #exactRankTests (package for computing exact rank tests with tied observations)
        #TODO: switch to 'coin' package
        try:
            r_statement = ('library("exactRankTests")');
            ans = robjects.r(r_statement);
            r_statement = ('require(exactRankTests)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("exactRankTests",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("exactRankTests")');
                ans = robjects.r(r_statement);
                r_statement = ('require(exactRankTests)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #rpart (package for decision trees)
        try:
            r_statement = ('library("rpart")');
            ans = robjects.r(r_statement);
            r_statement = ('require(rpart)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("rpart",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("rpart")');
                ans = robjects.r(r_statement);
                r_statement = ('require(rpart)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #e1071 (package for svm, naive bayes)
        try:
            r_statement = ('library("e1071")');
            ans = robjects.r(r_statement);
            r_statement = ('require(e1071)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("e1071",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("e1071")');
                ans = robjects.r(r_statement);
                r_statement = ('require(e1071)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #class (package for K-Nearest Neighbors)
        try:
            r_statement = ('library("class")');
            ans = robjects.r(r_statement);
            r_statement = ('require(class)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("class",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("class")');
                ans = robjects.r(r_statement);
                r_statement = ('require(class)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #cluster (package for kmeans)
        try:
            r_statement = ('library("cluster")');
            ans = robjects.r(r_statement);
            r_statement = ('require(cluster)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("cluster",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("cluster")');
                ans = robjects.r(r_statement);
                r_statement = ('require(cluster)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #randomForest (package for randomForest)
        try:
            r_statement = ('library("randomForest")');
            ans = robjects.r(r_statement);
            r_statement = ('require(randomForest)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("randomForest",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("randomForest")');
                ans = robjects.r(r_statement);
                r_statement = ('require(randomForest)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #DAAG (package for cross validation of linear models)
        try:
            r_statement = ('library("DAAG")');
            ans = robjects.r(r_statement);
            r_statement = ('require(DAAG)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("DAAG",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("DAAG")');
                ans = robjects.r(r_statement);
                r_statement = ('require(DAAG)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #boot (package for cross validation of glm)
        try:
            r_statement = ('library("boot")');
            ans = robjects.r(r_statement);
            r_statement = ('require(boot)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('install.packages("boot",dependencies=TRUE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("boot")');
                ans = robjects.r(r_statement);
                r_statement = ('require(boot)');
                ans = robjects.r(r_statement);
            except Exception as e:
                print(e);
        #LMGene (glog)
        try:
            r_statement = ('library("LMGene")');
            ans = robjects.r(r_statement);
            r_statement = ('require(LMGene)');
            ans = robjects.r(r_statement);
        except:
            try:
                r_statement = ('source("https://bioconductor.org/biocLite.R")');
                ans = robjects.r(r_statement);
                r_statement = ('biocLite("LMGene",ask=FALSE)');
                ans = robjects.r(r_statement);
                r_statement = ('library("LMGene")');
                ans = robjects.r(r_statement);
                r_statement = ('require(LMGene)');
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

