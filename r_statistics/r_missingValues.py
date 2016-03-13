from .r_dependencies import *
from .r_base import r_base
class r_missingValues(r_base):
    def calculate_missingValues(self,data_I,n_imputations_I = 1000):
        '''calculate missing values using a bootstrapping approach
        as implemented in the AmeliaII R package
        https://cran.r-project.org/web/packages/Amelia/Amelia.pdf
        http://r.iq.harvard.edu/docs/amelia/amelia.pdf
        1000 imputations (default) are computed and averaged to generate
        the resulting data without missing values
        INPUT:
        data_I = listDict
        n_imputations_I = integer, number of imputations
        OUTPUT:
        sns_O = [] of string, sample_name_short
        cn_O = [] of string, component_name
        cc_O = [] of float, calculated_concentration
        '''

        # convert data dict to matrix filling in missing values
        # with 'NA'
        sns = []
        cgn = []
        for d in data_I:
                sns.append(d['sample_name_short']);
                cgn.append(d['component_name']);
        sns_sorted = sorted(set(sns))
        cgn_sorted = sorted(set(cgn))
        concentrations = ['NA' for r in range(len(sns_sorted)*len(cgn_sorted))];
        cnt = 0;
        for c in cgn_sorted:
                for s in sns_sorted:
                    for d in data_I:
                        if d['sample_name_short'] == s and d['component_name'] == c:
                            if d['calculated_concentration']:
                                concentrations[cnt] = d['calculated_concentration'];
                                break;
                    cnt = cnt+1

        sns_O = [];
        cn_O = [];
        cc_O = [];
        # check if there were any missing values in the data set in the first place
        mv = 0;
        for c in concentrations:
            if c=='NA':
                mv += 1;
        if mv>0:
            # Call to R
            try:
                # convert lists to R objects
                # concentrations_R = robjects.FloatVector(concentrations);
                concentrations_r = '';
                for c in concentrations:
                    concentrations_r = (concentrations_r + ',' + str(c));
                concentrations_r = concentrations_r[1:];
                r_statement = ('concentrations = c(%s)' % concentrations_r);
                ans = robjects.r(r_statement);
                r_statement = 'concentrations_log = log(concentrations)'
                ans = robjects.r(r_statement);
                r_statement = ('concentrations_m = matrix(concentrations_log, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cgn_sorted),len(sns_sorted)));
                ans = robjects.r(r_statement);
                r_statement = ('a.out = amelia(concentrations_m, m=%s)' % n_imputations_I);
                ans = robjects.r(r_statement);
                # extract out data matrices
                concentrations_2d = numpy.zeros((len(cgn_sorted)*len(sns_sorted),n_imputations_I));
                for n in range(n_imputations_I):
                    cnt_data = 0;
                    for d in ans.rx2('imputations')[0]:
                        #NOTE: the matrix is linearized by column (opposite of how the matrix was made)
                        concentrations_2d[cnt_data,n] = exp(d);
                        cnt_data = cnt_data + 1;
                # calculate the average
                concentrations_1d_ave = numpy.zeros((len(cgn_sorted)*len(sns_sorted)));
                cnt_imputations = 0;
                for n in range(len(cgn_sorted)*len(sns_sorted)):
                    concentrations_1d_ave[n] = numpy.average(concentrations_2d[n][:]);
                # convert array back to dict
                #data_O = [];
                for c in range(len(sns_sorted)):
                    for r in range(len(cgn_sorted)):
                        #data_tmp = {};
                        #data_tmp['sample_name_short'] = sns_sorted[c]
                        #data_tmp['component_name'] = cgn_sorted[r]
                        #data_tmp['calculated_concentration'] = concentrations_1d_ave[c*len(sns_sorted)+r];
                        #data_O.append(data_tmp);
                        if isinstance(concentrations_1d_ave[c*len(cgn_sorted)+r], (int, float)) and not numpy.isnan(concentrations_1d_ave[c*len(cgn_sorted)+r]):
                            sns_O.append(sns_sorted[c]);
                            cn_O.append(cgn_sorted[r]);
                            cc_O.append(concentrations_1d_ave[c*len(cgn_sorted)+r]);

            # expand the array
            #concentrations_2d_ave = numpy.zeros((len(cgn_sorted),len(sns_sorted))); # transpose of input
            #for c in range(len(sns_sorted)):
            #    for r in range(len(cgn_sorted)):
            #        concentrations_2d_ave[r][c] = concentrations_1d_ave[c*len(sns_sorted)+r];

            except Exception as e:
                print(e);
                exit(-1);
        else:
            for r in range(len(cgn_sorted)):
                for c in range(len(sns_sorted)):
                    #data_O.append(data_tmp);
                    sns_O.append(sns_sorted[c]);
                    cn_O.append(cgn_sorted[r]);
                    cc_O.append(concentrations[r*len(sns_sorted)+c]);

        # reformat the matrix of element averages to dict
        return sns_O,cn_O,cc_O;
    def calculate_missingValues_pcaMethods(self,
            data_I,
            pca_model_I,
            pca_method_I,):
        '''Calculate missing values using pca methods
        INPUT:
        pca_method_I = string, pcaMethods method to use
            llsImput = missing value estimation using local least squares
            nipals
            bpca
            ppca
            svdImpute

        
        '''
        try:
            r_statement = ('imputed <- completeObs(pc)' % n_imputations_I);
            ans = robjects.r(r_statement);
            # convert array back to dict
            data_O = [];
            for c in range(len(sns_sorted)):
                for r in range(len(cgn_sorted)):
                    data_tmp = {};
                    data_tmp['sample_name_short'] = sns_sorted[c]
                    data_tmp['component_name'] = cgn_sorted[r]
                    data_tmp['calculated_concentration'] = None;

        except Exception as e:
            print(e);
            exit(-1);
        return data_O;

    def pcaMethods_llsImpute(self,
            input_var_I,
            output_var_I,
            k,
            center,
            imputeMissingValues,
            correlation,
            allVariables="TRUE",
            maxSteps="100"):
        '''Missing value estimationg using local least squares (LLS)
        INPUT:
        input_var_I = name of the R workspace matrix variable
            variables (e.g., genes) in rows; observations (e.g., samples) in rows
            missing values denoted as "NA"
        k = string int, cluster size
        center = string boolean, default, TRUE
        imputeMissingValues = string boolean, default, TRUE
        correlation = string, "pearson","kendall","spearman"
        allVariables = string boolean, default=FALSE

        OUTPUT:
        output_var_I = name fo the R workspace variable
        
        '''
        try:
            r_statement = ('%s <- llsImpute(%s, k = %s, completeObs=%s, correlation="%s", allVariables=%s)' %(
                output_var_I,input_var_I,k,center,imputeMissingValues,correlation,allVariables,maxSteps));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_pcaMethods_missingValues(self,
                pcamethods_var_I,
                imputed_O
                ):
        '''Imput missing data from pcamethods object
        INPUT:
        pcamethods_var_I = name of the R pcamethods object variable
        OUTPUT:
        imputed_O = R workspace variable
        data_O = numpy array of the full dataset
        '''
        try:
            r_statement = ('%s <- completeObs(%s)' %(imputed_O,pcamethods_var_I));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;

        except Exception as e:
            print(e);
            exit(-1);

    
    def calculate_missingValues(self,data_I,n_imputations_I = 1000,geometric_imputation_I=True):
        '''calculate missing values using a bootstrapping approach
        as implemented in the AmeliaII R package
        https://cran.r-project.org/web/packages/Amelia/Amelia.pdf
        http://r.iq.harvard.edu/docs/amelia/amelia.pdf
        1000 imputations (default) are computed and averaged to generate
        the resulting data without missing values
        INPUT:
        data_I = listDict
        n_imputations_I = integer, number of imputations
        geometric_imputation_I = boolean, apply a log normalization prior to imputing missing values
                    default: True, (better performance found for data that has not been previously logx transformed)
        OUTPUT:
        data_O = listDict
        '''

        ## convert data dict to matrix filling in missing values
        ## with 'NA'
        ##SPLIT1:
        #listdict = listDict(data_I);
        #concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList(
        #    row_label_I='component_name',
        #    column_label_I='sample_name_short',
        #    value_label_I='calculated_concentration',
        #    row_variables_I=['component_group_name'],
        #    column_variables_I=[
        #        'time_point','experiment_id'],
        #    data_IO=[],
        #    na_str_I="NA");
        #cgn = row_variables['component_group_name'];
        #tps = column_variables['time_point']
        #eis = column_variables['experiment_id']
        ## check if there were any missing values in the data set in the first place
        #mv = 0;
        #mv = listdict.count_missingValues(concentrations,na_str_I="NA");

        #SPLIT2:
        #define constants that should probably be inputs
        row_variables_I = ['component_name','component_group_name'];
        column_variables_I = ['sample_name_short','time_point','experiment_id'];
        value_label_I ='calculated_concentration';
        na_str_I="NA";
        #make the pandas dataframe
        listdict_pd = pd.DataFrame(data_I);
        listdict_pivot = listdict_pd.pivot_table(values=value_label_I,index = row_variables_I,columns = column_variables_I)
        #check for missing values
        mv = listdict_pivot.size - listdict_pivot.count().get_values().sum();
        #fill values with 'NA', convert to 1d numpy array, convert to list
        concentrations = list(listdict_pivot.fillna(na_str_I).get_values().ravel());
        #extract out rows and column variables
        #group = list(listdict_pd.groupby(row_variables).groups.keys());
        #group.sort();
        #cn_sorted = [g[0] for g in group];
        #cgn = [g[1] for g in group];
        row_variables = {}
        for i,rv in enumerate(row_variables_I):
            row_variables[rv] = [];
            for g in listdict_pivot.index.unique():
                row_variables[rv].append(g[i]);
        cn_sorted = row_variables['component_name']
        cgn = row_variables['component_group_name']
        # columns are in the same order as they were initialized during the pivot
        #sns_sorted = [v[0] for v in listdict_pivot.columns.get_values()];
        #tps = [v[1] for v in listdict_pivot.columns.get_values()];
        #eis = [v[2] for v in listdict_pivot.columns.get_values()];
        column_variables = {}
        for i,cv in enumerate(column_variables_I):
            column_variables[cv] = [];
            for g in listdict_pivot.columns.unique():
                column_variables[cv].append(g[i]);
        sns_sorted = column_variables['sample_name_short'];
        tps = column_variables['time_point'];
        eis = column_variables['experiment_id'];

        if mv>0:
            # Call to R
            try:
                # clear the workspace
                self.clear_workspace();
                #TODO: refactor into individual function
                # convert lists to R objects
                # concentrations_R = robjects.FloatVector(concentrations);
                concentrations_r = '';
                for c in concentrations:
                    concentrations_r = (concentrations_r + ',' + str(c));
                concentrations_r = concentrations_r[1:];
                r_statement = ('concentrations = c(%s)' % concentrations_r);
                ans = robjects.r(r_statement);
                if geometric_imputation_I:
                    r_statement = 'concentrations_log = log(concentrations)'
                    ans = robjects.r(r_statement);
                r_statement = ('concentrations_m = matrix(concentrations_log, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted)));
                ans = robjects.r(r_statement);
                r_statement = ('a.out = amelia(concentrations_m, m=%s)' % n_imputations_I);
                ans = robjects.r(r_statement);

                ##SPLIT1
                ## extract out data matrices
                #concentrations_2d = numpy.zeros((len(cn_sorted)*len(sns_sorted),n_imputations_I));
                #for n in range(n_imputations_I):
                #    cnt_data = 0;
                #    for d in ans.rx2('imputations')[0]:
                #        #NOTE: the matrix is linearized by column (opposite of how the matrix was made)
                #        if geometric_imputation_I:
                #            concentrations_2d[cnt_data,n] = exp(d);
                #        else:
                #            concentrations_2d[cnt_data,n] = d;
                #        cnt_data = cnt_data + 1;
                ## calculate the average
                #concentrations_1d_ave = numpy.zeros((len(cn_sorted)*len(sns_sorted)));
                #cnt_imputations = 0;
                #for n in range(len(cn_sorted)*len(sns_sorted)):
                #    concentrations_1d_ave[n] = numpy.average(concentrations_2d[n][:]);

                #SPLIT 2:
                # extract out data matrices
                imputations = numpy.array([numpy.array(i) for i in ans.rx2('imputations')]); #dim imputations, rows, columns
                if geometric_imputation_I: imputations = numpy.exp(imputations);
                concentrations_2d = numpy.mean(imputations,axis=0); #dim rows, columns
                # convert array back to dict
                data_O = [];
                for c in range(len(sns_sorted)):
                    for r in range(len(cn_sorted)):
                        #concentration_imputed = concentrations_1d_ave[c*len(cn_sorted)+r];
                        concentration_imputed = concentrations_2d[r,c]
                        if isinstance(concentration_imputed, (int, float)) and not numpy.isnan(concentration_imputed):
                            #sns_O.append(sns_sorted[c]);
                            #cn_O.append(cn_sorted[r]);
                            #cc_O.append(concentrations_1d_ave[c*len(cn_sorted)+r]);
                            data_tmp = {};
                            data_tmp['analysis_id'] = data_I[0]['analysis_id']
                            data_tmp['experiment_id'] = eis[c]
                            data_tmp['time_point'] = tps[c]
                            data_tmp['sample_name_short'] = sns_sorted[c]
                            data_tmp['component_name'] = cn_sorted[r]
                            data_tmp['component_group_name'] = cgn[r]
                            data_tmp['calculated_concentration'] = concentration_imputed
                            data_tmp['calculated_concentration_units'] = data_I[0]['calculated_concentration_units']
                            data_tmp['imputation_method']='ameliaII';
                            data_tmp['imputation_options']={
                                'n_imputations':n_imputations_I,
                                'geometric_imputation':geometric_imputation_I}
                            data_O.append(data_tmp);
                        else:
                            print('NaN found in AmeliaII output.');
                # expand the array
                #concentrations_2d_ave = numpy.zeros((len(cn_sorted),len(sns_sorted))); # transpose of input
                #for c in range(len(sns_sorted)):
                #    for r in range(len(cn_sorted)):
                #        concentrations_2d_ave[r][c] = concentrations_1d_ave[c*len(sns_sorted)+r];
                return data_O;
            except Exception as e:
                print(e);
                exit(-1);
        else:
            return None;
