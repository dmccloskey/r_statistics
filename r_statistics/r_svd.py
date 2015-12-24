from .r_dependencies import *
from .r_base import r_base
from math import sqrt
class r_svd(r_base):
    def calculate_pcaMethods_svd(self,
            data_var_I,
            data_var_O,):
        '''calculate svd using pcaMethods'''
        try:
            r_statement = ('%s <- svd(%s)' %(
                data_var_O,data_var_I))
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_pcaMethods_robustSvd(self,
            data_var_I,
            data_var_O,):
        '''calculate robustSvd using pcaMethods
        The robust SVD of the matrix is x = u d v'.
            d A vector containing the singular values of x.
            u A matrix whose columns are the left singular vectors of x.
            v A matrix whose columns are the right singular vectors of x
        
        '''
        try:
            r_statement = ('%s <- robustSvd(%s)' %(
                data_var_O,data_var_I))
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def extract_svd(self,
            data_var_I,):
        '''extract svd matrices

        d	
        a vector containing the singular values of x, of length min(n, p).
        u	
        a matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu)
        v	
        a matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension c(p, nv)
        
        '''
        try:
            r_statement = ('%s$u' %(data_var_I))
            ans = robjects.r(r_statement);
            u = numpy.array(ans);
            #r_statement = ('diag(%s$d)' %(data_var_I))
            r_statement = ('%s$d' %(data_var_I))
            ans = robjects.r(r_statement);
            d = numpy.array(ans);
            r_statement = ('%s$v' %(data_var_I))
            ans = robjects.r(r_statement);
            v = numpy.array(ans);
            return u,d,v;
        except Exception as e:
            print(e);
            exit(-1);
    
    def detect_outliers_svd(self,data_I):
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
        
        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=['sample_name_abbreviation'],
            data_IO=[],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        sna = column_variables['sample_name_abbreviation'];
        nsna_unique,sna_unique = listdict.get_uniqueValues('sample_name_abbreviation');
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues(concentrations,na_str_I="NA");
        if mv==0:
            # Call to R
            try:
                # clear the workspace
                self.clear_workspace();
                # convert lists to R matrix
                self.make_matrixFromList(concentrations,len(cn_sorted),len(sns_sorted),'concentrations_m');
                # convert to the transpose
                self.transpose_matrix('concentrations_mt','concentrations_m');
                # calculate svd
                self.calculate_pcaMethods_svd('concentrations_mt','resSVD');
                # extract svd matrices
                u,d,v = self.extract_svd('resSVD');
                # calculate robustSVD
                self.calculate_pcaMethods_robustSvd('concentrations_mt','resRobustSVD');
                # extract out robustSVD matrices
                ur,dr,vr = self.extract_svd('resRobustSVD');

                #plot(resSvd$v[,1], resSvdOut$v[,1])
                #plot(resRobSvd$v[,1], resRobSvdOut$v[,1])
                #s$u %*% D %*% t(s$v) #  X = U D V'
                #robjects.r('resSVD$u %*% diag(resSVD$d) %*% t(resSVD$v)')
                #import numpy as np
                #import matplotlib.pyplot as plt

                #plt.scatter(v[:,1], vr[:,1]);
                #plt.show()
                #x = numpy.array(robjects.r('concentrations_mt'))
                #xsvd = numpy.multiply(numpy.multiply(u,numpy.diag(d)),numpy.transpose(v));

            except Exception as e:
                print(e);
                exit(-1);
            return u,d,v,ur,dr,vr;
        else:
            print('missing values found!');
    def calculate_svd(self,data_I,svd_method_I,svd_options_I={}):
        '''calculate svd
        '''
        
        u_O,d_O,v_O = [],[],[];
        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=['sample_name_abbreviation'],
            data_IO=[],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        sna = column_variables['sample_name_abbreviation'];
        nsna_unique,sna_unique = listdict.get_uniqueValues('sample_name_abbreviation');
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues(concentrations,na_str_I="NA");
        if mv==0:
            # Call to R
            try:
                # clear the workspace
                self.clear_workspace();
                # convert lists to R matrix
                self.make_matrixFromList(concentrations,len(cn_sorted),len(sns_sorted),'concentrations_m');
                # convert to the transpose
                self.transpose_matrix('concentrations_mt','concentrations_m');
                # calculate svd
                if svd_method_I == 'svd':
                    self.calculate_pcaMethods_svd('concentrations_mt','resSVD');
                elif svd_method_I == 'robustSvd':
                    self.calculate_pcaMethods_robustSvd('concentrations_mt','resSVD');
                else:
                    print('svd method not recognized');
                    return u_O,d_O,v_O;
                # extract svd matrices
                u,d,v = self.extract_svd('resSVD');
                # reformat svd output
                u_O,d_O,v_O = self.reformat_udv(u,d,v,cn_sorted,sns_sorted,cgn,sna,svd_method_I,svd_options_I);
            except Exception as e:
                print(e);
                exit(-1);
            return u_O,d_O,v_O;
        else:
            print('missing values found!');
            return u_O,d_O,v_O ;

    def reformat_udv(self,u,d,v,cn,sns,cgn,sna,svd_method_I,svd_options_I):
        '''reformat SVD u, d, and v matrices to listDicts
        INPUT:
        OUTPUT
        u_O
        d_O
        v_O
        '''
        # extractou out the U matrix
        u_O = [];
        cnt=0;
        for r in range(u.shape[0]):
            for c in range(u.shape[1]):
                data_tmp = {};
                data_tmp['sample_name_short'] = sns[r];
                data_tmp['sample_name_abbreviation'] = sna[r];
                data_tmp['u_matrix'] = u[r,c];
                data_tmp['singular_value_index'] = c+1;
                data_tmp['svd_method'] = svd_method_I;
                data_tmp['svd_options'] = {
                    };
                u_O.append(data_tmp);
                cnt+=1;
        # extract out the V matrix
        v_O = [];
        cnt=0;
        for r in range(v.shape[0]):
            for c in range(v.shape[1]): #comp
                data_tmp = {};
                data_tmp['component_name'] = cn[r]; #need to double check
                data_tmp['component_group_name'] = cgn[r];
                data_tmp['v_matrix'] = v[r,c]; #need to double check
                data_tmp['singular_value_index'] = c+1;
                data_tmp['svd_method'] = svd_method_I;
                data_tmp['svd_options'] = {
                    };
                v_O.append(data_tmp);
                cnt+=1;
        # extract out d vector
        d_O = [];
        cnt=0;
        for r in range(d.shape[0]):
            data_tmp = {};
            data_tmp['d_vector'] = d[r];
            data_tmp['d_fraction'] = d[r]/numpy.sum(d);
            data_tmp['singular_value_index'] = r+1;
            data_tmp['svd_method'] = svd_method_I;
            data_tmp['svd_options'] = {
                };
            d_O.append(data_tmp);
            cnt+=1;
        return u_O,d_O,v_O;