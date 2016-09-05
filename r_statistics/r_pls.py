from .r_dependencies import *
from .r_base import r_base
from math import sqrt
class r_pls(r_base):
    def calculate_plsda_mixomics(self,data_I,factor_I= "sample_name_abbreviation",
                        ncomp=5,max_iter=500,
                        tol = 1e-3, near_zero_var = "TRUE",
                        method_predict="all",validation="Mfold",
                        folds = 10, progressBar = "FALSE"):
        '''Perform PLS-DA

        mixOmics pslda methods:
        plsda(X, Y, ncomp = 3, max.iter = 500, tol = 1e-06, near.zero.var = TRUE)
        perf(object,
            method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
            validation = c("Mfold", "loo"),
            folds = 10, progressBar = FALSE, near.zero.var = FALSE, ...)
        predict(object, newdata, method = c("all", "max.dist",
            "centroids.dist", "mahalanobis.dist"), ...)

        INPUT:
        data_I = [{}], list of data dicts from .csv or a database
        factor_I = column of the data dict to use as the factor
                (default = "sample_name_abbreviation")

        STATUS: Working, need to finalize extracting the scores, loadings, and model performance for output

        '''
        
        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList_pd(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=[factor_I],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        factor = column_variables[factor_I];
        factors_unique = listdict.get_uniqueValues_list(factor_I);
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues_pivotTable();
        if mv==0:
            # Call to R
            try:
                # convert lists to R matrix
                concentrations_r = '';
                for c in concentrations:
                    concentrations_r = (concentrations_r + ',' + str(c));
                concentrations_r = concentrations_r[1:];
                r_statement = ('concentrations = c(%s)' % concentrations_r);
                ans = robjects.r(r_statement);
                r_statement = ('concentrations_m = matrix(concentrations, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted))); 
                ans = robjects.r(r_statement);
                # convert lists to R factors
                factor_r = '';
                for c in factor:
                    factor_r = (factor_r + ',' + '"' + str(c) + '"');
                factor_r = factor_r[1:];
                r_statement = ('factors = c(%s)' % factor_r);
                ans = robjects.r(r_statement);
                r_statement = ('factors_m = matrix(factors, nrow = %s, ncol = %s, byrow = TRUE)' %(len(sns_sorted),1));
                ans = robjects.r(r_statement);
                r_statement = ('factors_f = factor(factors_m)'); 
                ans = robjects.r(r_statement);
                # call plsda
                r_statement = ('result = plsda(t(concentrations_m), factors_f,\
                    ncomp=%s, max.iter=%s, tol=%s, near.zero.var=%s)' %(ncomp,max_iter,tol,near_zero_var));
                    #need to send in the transpose
                ans = robjects.r(r_statement);
                ##plot the plsda results
                #r_statement = ('plotIndiv(result, ind.names = factors_f, plot.ellipse = TRUE, add.legend =TRUE)');
                #ans = robjects.r(r_statement);
                #get the plsda results
                #get the variates
                variates = ans.rx2("variates"); 
                #get the X and Y variates (i.e, scores in pca)
                #columns of X are the principle variants
                scores_x = np.array(variates.rx2('X'));# variates_x.shape
                scores_y = np.array(variates.rx2('Y'));# variates_y.shape
                #get the loadings
                loadings = ans.rx2("loadings")
                #get the X and Y loadings
                #columns of X are the principle loadings
                loadings_x = np.array(loadings.rx2('X'))
                loadings_y = np.array(loadings.rx2('Y'));# variates_y.shape
                # VIP
                r_statement = ('result_vip = vip(result)');
                ans = robjects.r(r_statement);
                # get the VIPs
                #columns of vip are the vip's of the loadings
                vip = np.array(ans);
                # cross validate the plsda model
                r_statement = ('result_perf = perf(result, method.predict="%s",\
                    validation="%s", folds=%s, progressBar = %s, near.zero.var=%s)'\
                        %(method_predict,validation,folds,progressBar,"FALSE"));
                ans = robjects.r(r_statement);
                # get the results of the cross validation
                # columns error_rate are the principle variants
                # rows are the methods used
                #      "all" = "max.dist","centroids.dist", "mahalanobis.dist"
                error_rate = np.array(ans.rx2('error.rate'));
                error_rate_var = np.zeros_like(error_rate);
                for row_cnt in range(error_rate.shape[0]):
                    for col_cnt in range(error_rate.shape[1]):
                        if row_cnt == 0:
                            error_rate_var[row_cnt,col_cnt]=error_rate[row_cnt,col_cnt];
                        else:
                            error_rate_var[row_cnt,col_cnt]=error_rate[row_cnt,col_cnt]-error_rate[row_cnt-1,col_cnt];
                # extract out scores
                data_scores = [];
                cnt=0;
                for r in range(scores_x.shape[0]):
                    for c in range(scores_x.shape[1]):
                        data_tmp = {};
                        data_tmp['sample_name_short'] = sns_sorted[r];
                        data_tmp['score'] = scores_x[r,c];
                        data_tmp['axis'] = c+1;
                        #data_tmp['var_proportion'] = var_proportion[c];
                        #data_tmp['var_cumulative'] = var_cumulative[c];
                        data_tmp['error_rate'] = error_rate[c,0];
                        data_tmp['factor'] = factor[r];
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_scores.append(data_tmp);
                        cnt+=1;
                # extract out loadings
                data_loadings = [];
                cnt=0;
                for r in range(loadings_x.shape[0]):
                    for c in range(loadings_x.shape[1]):
                        data_tmp = {};
                        data_tmp['component_name'] = cn_sorted[r]; #need to double check
                        data_tmp['component_group_name'] = cgn[r];
                        data_tmp['loadings'] = loadings_x[r,c]; #need to double check
                        data_tmp['axis'] = c+1;
                        data_tmp['vip'] = vip[r,0];
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_loadings.append(data_tmp);
                        cnt+=1;
                # extract out the model performance statistics
                data_perf = [];
            except Exception as e:
                print(e);
                exit(-1);
            return data_scores,data_loadings
        else:
            print('missing values found!')

    def calculate_perf_mixOmics(self,
            perf_mixOmics_O = 'perf.o',
            near_zero_var = "FALSE",
            method_predict="all",validation="Mfold",
            folds = 10, progressBar = "FALSE"
            ):
        '''CV a mixOmics model'''
        try:
            r_statement = ('%s = perf(result, method.predict="%s",\
                validation="%s", folds=%s, progressBar = %s, near.zero.var=%s)'\
                    %(cv_mixOmics_O,method_predict,validation,folds,progressBar,"FALSE"));
            ans = robjects.r(r_statement);
            
        except Exception as e:
            print(e);
            exit(-1);

    def extract_per_mixOmics(self,
            perf_mixOmics_I = 'perf.o',
            ):
        '''Extract the perf of a mixOmics model'''
        try:
            r_statement = ('%s'
                % (perf_mixOmics_I));
            ans = robjects.r(r_statement);
            # get the results of the cross validation
            # columns error_rate are the principle variants
            # rows are the methods used
            #      "all" = "max.dist","centroids.dist", "mahalanobis.dist"
            error_rate = np.array(ans.rx2('error.rate'));
            error_rate_var = np.zeros_like(error_rate);
            for row_cnt in range(error_rate.shape[0]):
                for col_cnt in range(error_rate.shape[1]):
                    if row_cnt == 0:
                        error_rate_var[row_cnt,col_cnt]=error_rate[row_cnt,col_cnt];
                    else:
                        error_rate_var[row_cnt,col_cnt]=error_rate[row_cnt,col_cnt]-error_rate[row_cnt-1,col_cnt];
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_oplsda_ropl(self,data_I,factor_I= "sample_name_abbreviation",
            predI = "NA",
            orthoI = 0,
            algoC = "default",
            crossvalI = 10,
            log10L = "FALSE",
            permI = 10,
            scaleC = "center",
            subset = "NULL",
            printL = "FALSE",
            plotL = "FALSE"):
        '''Perform (O)PLS(-DA)
        rpols methods:
        http://bioconductor.org/packages/release/bioc/manuals/ropls/man/ropls.pdf
        opls(x,
            y = NULL,
            predI = NA,
            orthoI = 0,
            algoC = c("default", "nipals", "svd")[1],
            crossvalI = 7,
            log10L = FALSE,
            permI = 10,
            scaleC = c("center", "pareto", "standard")[3],
            subset = NULL,
            printL = TRUE,
            plotL = TRUE,
            .sinkC = NULL,
            ...)

        INPUT:
        data_I = [{}], list of data dicts from .csv or a database
        factor_I = column of the data dict to use as the factor for (O)-PLS-DA
                (default = "sample_name_abbreviation")
                if "NULL" pca is performed
        orthoI = 0, PLS
                 "NA", 1, 2, 3, ..., OPLS

        NOTE: supports only two factors
        USE CASE: biomarker discover (control vs. 1 condition);

        STATUS: Working, need to finalize extracting the scores, loadings, and model performance for output

        '''
        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList_pd(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=[factor_I],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        factor = column_variables[factor_I];
        factors_unique = listdict.get_uniqueValues_list(factor_I);
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues_pivotTable();
        if mv==0:
            # Call to R
            try:
                # convert lists to R matrix
                concentrations_r = '';
                for c in concentrations:
                    concentrations_r = (concentrations_r + ',' + str(c));
                concentrations_r = concentrations_r[1:];
                r_statement = ('concentrations = c(%s)' % concentrations_r);
                ans = robjects.r(r_statement);
                r_statement = ('concentrations_m = matrix(concentrations, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted))); 
                ans = robjects.r(r_statement);
                # convert lists to R factors
                factor_r = '';
                for c in factor:
                    factor_r = (factor_r + ',' + '"' + str(c) + '"');
                factor_r = factor_r[1:];
                r_statement = ('factors = c(%s)' % factor_r);
                ans = robjects.r(r_statement);
                #r_statement = ('factors_m = matrix(factors, nrow = %s, ncol = %s, byrow = TRUE)' %(len(sns_sorted),1));
                #ans = robjects.r(r_statement);
                #r_statement = ('factors_f = factor(factors_m)'); 
                #ans = robjects.r(r_statement);
                r_statement = ('factors_f = factor(factors)'); 
                ans = robjects.r(r_statement);
                # call (o)plsda
                #TODO
                r_statement = ('result = opls(t(concentrations_m), y=factors_f,\
                        predI = %s,\
                        orthoI = %s,\
                        algoC = %s,\
                        crossvalI = %s,\
                        log10L = %s,\
                        permI = %s,\
                        scaleC = %s,\
                        subset = %s,\
                        printL = %s,\
                        plotL = %s)' %(predI,orthoI,algoC,crossvalI,log10L,permI,scaleC,subset,printL,plotL));
                    #need to send in the transpose
                ans = robjects.r(r_statement);
                #get the plsda results
                #get the scores
                scores_x = np.array(ans.rx2('scoreMN'));
                #get the loadings
                loadings_x = np.array(ans.rx2('loadingMN'));
                if orthoI !=0:
                    #get the scores
                    scoresortho_x = np.array(ans.rx2('orthoScoreMN'));
                    #get the loadings
                    loadingsortho_x = np.array(ans.rx2('orthoLoadingMN'));
                # get the VIPs
                vip = np.array(ans.rx2("vipVn")); 
                ans = robjects.r(r_statement);
                # get the variance of each component
                var_std_x = np.array(ans.rx2("xSdVn")); 
                var_std_y = np.array(ans.rx2("ySdVn")); 
                var_proportion = np.power(var_std_x,2);
                var_proportion_y = np.power(var_std_x,2);
                # get model summary
                #TODO:
                #---------------------------------
                modelDF = ans.rx2("modelDF");
                summaryDF = ans.rx2("summaryDF");
                # extract out scores
                data_scores = [];
                cnt=0;
                for r in range(scores_x.shape[0]):
                    for c in range(scores_x.shape[1]):
                        data_tmp = {};
                        data_tmp['sample_name_short'] = sns_sorted[r];
                        data_tmp['score'] = scores_x[r,c];
                        data_tmp['axis'] = c+1;
                        #data_tmp['var_proportion'] = var_proportion[c];
                        #data_tmp['var_cumulative'] = var_cumulative[c];
                        data_tmp['error_rate'] = error_rate[c,0];
                        data_tmp['factor'] = factor[r];
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_scores.append(data_tmp);
                        cnt+=1;
                # extract out loadings
                data_loadings = [];
                cnt=0;
                for r in range(loadings_x.shape[0]):
                    for c in range(loadings_x.shape[1]):
                        data_tmp = {};
                        data_tmp['component_name'] = cn_sorted[r]; #need to double check
                        data_tmp['component_group_name'] = cgn[r];
                        data_tmp['loadings'] = loadings_x[r,c]; #need to double check
                        data_tmp['axis'] = c+1;
                        data_tmp['vip'] = vip[r,0];
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_loadings.append(data_tmp);
                        cnt+=1;
                # extract out the model performance statistics
                data_perf = [];
                #---------------------------------
                #END TODO:
            except Exception as e:
                print(e);
                exit(-1);
            return data_scores,data_loadings
        else:
            print('missing values found!')
    def calculate_mvr(self,data_I,
            pls_model_I = 'PLS-DA',
            response_I = None,
            factor_I= "sample_name_abbreviation",
            ncomp = 5,
            Y_add = "NULL",
            scale = "TRUE",
            #validation = "LOO",
            validation = "CV",
            segments = 10,
            method = "cppls",
            stripped = "FALSE",
            lower = 0.5,
            upper = 0.5, 
            trunc_pow = "FALSE", 
            weights = "NULL",
            p_method = "fdr",
            nperm = 999):
        '''Perform PLS(-DA)

        pls is used for principle component regression, partial least squares regression, and Canonical powered partial least squares regression
            NOTE: oscorespls implements the NIPALS algorithm, which is also used in many other packages including scikit-learn
        pls is used for PLS-DA via converting the factor vector into a dummy response variable (the plsda function in caret automates this task)
        spls is used for sparse versions of pls
        caret provides utility functions for pls and spls (e.g., calculating the vip)


        mvr(formula, ncomp, Y.add, data, subset, na.action,
        method = pls.options()$mvralg,
        scale = FALSE, validation = c("none", "CV", "LOO"),
        model = TRUE, x = FALSE, y = FALSE, ...)
        plsr(..., method = pls.options()$plsralg)
        cppls(..., Y.add, weights, method = pls.options()$cpplsalg)
        pcr(..., method = pls.options()$pcralg)

        where formula = y ~ X and data = data.frame with y and X

        where segments is the number of segments for cross validation of type "CV"

        where methods = the algorithmn used 
            "kernelpls", "widekernelpls", "simpls": partial least squares regression
            "oscorespls": orthonogonal scores pls (uses the NIPALS algorithm)
                two-factor only
            "cppls": Canonical Powered Partial Least Squares 
                PLS
                PPLS
                PLS-DA (dummy discrete response y)
                PPLS-DA
                CPLS
                CPPLS
            "svdpc": principle component regression

        INPUT:
        data_I = [{}], list of data dicts from .csv or a database
        pls_model_I = name of the pls method (e.g., plsda)
        response_I = column of the data dict to use as the response factor for PLS
                    (general PLS has not yet been implemented, so response_I is not currently used)
        factor_I = column of the data dict to use as the factor for (O)-PLS-DA
                (default = "sample_name_abbreviation")

        TODO: add in permutation test
        MVA.test(X, Y, cmv = FALSE, ncomp = 5, kout = 7, kinn = 8, model = c("PLSR",
            "CPPLS", "PLS-DA", "PPLS-DA", "LDA", "QDA", "PLS-DA/LDA", "PLS-DA/QDA",
            "PPLS-DA/LDA","PPLS-DA/QDA"), Q2diff = 0.05, lower = 0.5, upper = 0.5,
            Y.add = NULL, weights = rep(1, nrow(X)), set.prior = FALSE,
            crit.DA = c("plug-in", "predictive", "debiased"), p.method = "fdr",
            nperm = 999,...)


        '''
        

        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList_pd(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=[factor_I],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        factor = column_variables[factor_I];
        factors_unique = listdict.get_uniqueValues_list(factor);
        # check the ncomp
        if len(factors_unique)<ncomp:
            ncomp = len(factors_unique);
        if len(factor)<segments:
            segments = int(0.75*len(factor));
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues_pivotTable();
        if mv==0:
            # Call to R
            try:
                # clear the workspace
                self.clear_workspace();
                # convert lists to R matrix
                self.make_matrixFromList(concentrations,len(cn_sorted),len(sns_sorted),'concentrations_m');
                # convert to the transpose
                self.transpose_matrix('concentrations_mt','concentrations_m');
                # convert lists to R factors
                if factor_I:
                    self.make_factorsFromList(factor,'factors_f');
                    # make the dummy response for PLS-DA
                    self.make_dummyMatrixFromFactors('factors_f','dummy');
                    # make the R dataframe
                    r_statement = ('dataframe <- data.frame(concentrations_mt = concentrations_mt, factors_f = factors_f)'); 
                    ans = robjects.r(r_statement);
                    r_statement = ('dataframe$dummy <- dummy'); 
                    ans = robjects.r(r_statement);
                elif response_I:
                    #TODO: response and responses_sorted
                    response_r = '';
                    for c in response:
                        response_r = (response_r + ',' + '"' + str(c) + '"');
                    response_r = response_r[1:];
                    r_statement = ('responses = c(%s)' % response_r);
                    ans = robjects.r(r_statement);
                    r_statement = ('responses_m = matrix(responses, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted)));  
                    ans = robjects.r(r_statement);
                    r_statement = ('responses_mt = t(concentrations_m)');  
                    ans = robjects.r(r_statement);
                    # make the R dataframe
                    r_statement = ('dataframe <- data.frame(concentrations_mt = t(concentrations_m), responses_mt = t(responses_m))');  
                    ans = robjects.r(r_statement);
                # make the formula
                if factor_I:
                    ##works as well (requires dummy and concentrations_m to be in the global environment)
                    #r_statement = ('fit <- lm(dummy ~ concentrations_m)');
                    #ans = robjects.r(r_statement);
                    fit = 'dummy ~ concentrations_mt';
                elif response_I:
                    fit = 'responses_mt ~ concentrations_mt';
                #call mvr
                self.call_mvr('result',fit,
                              ncomp=ncomp,scale=scale,validation=validation,
                              segments=segments,method=method,lower=lower,
                              upper=upper,weights=weights
                              );
                # plot the pls results
                #r_statement = ('plot(result,plottype = "scores")');
                #ans = robjects.r(r_statement);
                #get the plsda results
                data_scores,data_loadings,data_loadings_factors = self.calculate_mvr_scoresAndLoadings('result',
                            pls_model_I,
                            sns_sorted,
                            factor,
                            factors_unique,
                            cn_sorted,
                            cgn,
                            scale,
                            lower,
                            upper,
                            weights,
                            method,);
                #get the coefficients
                data_coefficients = self.calculate_mvr_coefficients('result',
                            pls_model_I,
                            factors_unique,
                            cn_sorted,
                            cgn,
                            scale,
                            lower,
                            upper,
                            weights,
                            method,
                            );
                # get the VIPs
                data_vip = self.calculate_mvr_vip('result',
                                pls_model_I,
                                factors_unique,
                                cn_sorted,
                                cgn,
                                scale,
                                lower,
                                upper,
                                weights,
                                method,);
                # get the model cross validation statistics
                data_perf = self.calculate_mvr_performance(
                            'result',
                            pls_model_I,
                            ncomp,
                            Y_add,
                            scale,
                            validation,
                            segments,
                            method,
                            stripped,
                            lower,
                            upper, 
                            trunc_pow, 
                            weights,
                            p_method,
                            nperm);
            except Exception as e:
                print(e);
                exit(-1);
            return data_scores,data_loadings,data_perf,data_vip,data_coefficients,data_loadings_factors
        else:
            print('missing values found!');

    def call_mvr(self,
            mvr_O = 'result',
            fit = 'dummy ~ concentrations_mt',      
            data = 'dataframe',           
            ncomp = 5,
            scale = "TRUE",
            validation = "CV",
            segments = 10,
            method = "cppls",
            lower = 0.5,
            upper = 0.5, 
            weights = "NULL",
            trunc_pow = "FALSE", 
            stripped = "FALSE",
            Y_add = "NULL",
                 ):
        '''call R mvr'''
        try:
            #call mvr
            if validation == "CV":
                r_statement = ('%s = mvr(formula = %s, ncomp = %s, data = %s, scale = %s, validation = "%s", segments = %s, method = "%s", lower = %s, upper = %s,  weights = %s)'\
                        %(mvr_O,fit,ncomp,data,scale,validation,segments,method,lower,upper,weights));
            elif validation == "LOO":
                ##works as well (requires dummy and concentrations_m to be in the global environment)
                #r_statement = ('result = mvr(fit, %s, data = dataframe, scale = %s, validation = "%s", method = "%s", lower = %s, upper = %s,  weights = %s)'\
                #        %(ncomp,scale,validation,method,lower,upper,weights));
                # (requires dummy and concentration to be named variables in the dataframe)
                r_statement = ('%s = mvr(formula = %s, ncomp = %s, data = %s, scale = %s, validation = "%s", method = "%s", lower = %s, upper = %s,  weights = %s)'\
                        %(mvr_O,fit,ncomp,data,scale,validation,method,lower,upper,weights));
            # adding Y.add and trunc.pow as input to mvr throws an error:
            #r_statement = ('%s = mvr(%s, %s,\
            #        data = dataframe,\
            #        Y.add = %s,\
            #        scale = %s,\
            #        validation = "%s",\
            #        segments = %s,\
            #        method = "%s",\
            #        stripped = %s,\
            #        lower = %s,\
            #        upper = %s,\
            #        trunc.pow = %s,\
            #        weights = %s)' %(mvr_O,fit,ncomp,Y_add,scale,validation,segments,method,stripped,lower,upper,trunc_pow,weights));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
            
    def calculate_mvr_vip(self,
            mvr_model_I,
            pls_model_I,
            factors_unique,
            cn_sorted,
            cgn,
            scale,
            lower,
            upper,
            weights,
            method,
            ):
        '''Calculate the VIPs of the mvr model

        requires caret

        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_vip

        '''
        try:
            vip,vip_reduced=self.extract_mvr_vip(mvr_model_I)
            # extract out VIP
            data_vip = [];
            for r in range(vip.shape[0]):
                for c in range(vip.shape[1]):
                    data_tmp = {};
                    data_tmp['response_name'] = factors_unique[r];
                    data_tmp['component_name'] = cn_sorted[c]; #need to double check
                    data_tmp['component_group_name'] = cgn[c];
                    data_tmp['pls_vip'] = vip[r,c];
                    data_tmp['pls_model'] = pls_model_I;
                    data_tmp['pls_method'] = method;
                    data_tmp['pls_options'] = {'pls_scale':scale,
                                                'lower':lower,
                                                'upper':upper,
                                                'weights':weights,
                        };
                    data_vip.append(data_tmp);
            for c in range(vip.shape[1]):
                data_tmp = {};
                data_tmp['response_name'] = 'all';
                data_tmp['component_name'] = cn_sorted[c]; #need to double check
                data_tmp['component_group_name'] = cgn[c];
                data_tmp['pls_vip'] = vip_reduced[c];
                data_tmp['pls_model'] = pls_model_I;
                data_tmp['pls_method'] = method;
                data_tmp['pls_options'] = {'pls_scale':scale,
                                            'lower':lower,
                                            'upper':upper,
                                            'weights':weights,
                    };
                data_vip.append(data_tmp);
        except Exception as e:
            print(e);
            exit(-1);
        return data_vip;
    def extract_mvr_vip(self,
            mvr_model_I,
            ):
        '''Calculate the VIPs of the mvr model

        requires caret

        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_vip

        '''
        try:
            # get the VIPs
            r_statement = ('varImp(%s)'%(mvr_model_I)); #requires caret
            ans = robjects.r(r_statement);
            vip = np.array(ans); #dim 1 = factors/responses, dim 2 = features
            vip_reduced = np.zeros_like(vip[0,:]);
            for j in range(vip.shape[1]):
                vip_reduced[j]=vip[:,j].sum();
            # extract out VIP
            return vip,vip_reduced;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_mvr_scoresAndLoadings(self,
            mvr_model_I,
            pls_model_I,
            sns_sorted,
            factor,factors_unique,
            cn_sorted,
            cgn,
            scale,
            lower,
            upper,
            weights,
            method,
            ):
        '''extract out mvr scores and loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            scores_x,scores_y,var_proportion,var_cumulative,loadings_x,loadings_y,cor_x,cor_y=self.extract_mvr_scoresAndLoadings(mvr_model_I)
            # extract out scores
            data_scores = [];
            cnt=0;
            for r in range(scores_x.shape[0]):
                for c in range(scores_x.shape[1]):
                    data_tmp = {};
                    data_tmp['sample_name_short'] = sns_sorted[r];
                    data_tmp['response_name'] = factor[r];
                    data_tmp['score'] = scores_x[r,c];
                    data_tmp['score_response'] = scores_y[r,c]
                    data_tmp['axis'] = c+1;
                    data_tmp['var_proportion'] = var_proportion[c];
                    data_tmp['var_cumulative'] = var_cumulative[c];
                    data_tmp['pls_model'] = pls_model_I;
                    data_tmp['pls_method'] = method;
                    data_tmp['pls_options'] = {'pls_scale':scale,
                                                'lower':lower,
                                                'upper':upper,
                                                'weights':weights,
                        };
                    #data_tmp['experiment_id'] = experiment_ids[cnt];
                    #data_tmp['time_point'] = time_points[cnt];
                    data_scores.append(data_tmp);
                    cnt+=1;
            # extract out loadings
            data_loadings = [];
            cnt=0;
            for r in range(loadings_x.shape[0]):
                for c in range(loadings_x.shape[1]): #comp
                    data_tmp = {};
                    data_tmp['component_name'] = cn_sorted[r]; #need to double check
                    data_tmp['component_group_name'] = cgn[r];
                    data_tmp['loadings'] = loadings_x[r,c]; #need to double check
                    data_tmp['axis'] = c+1;
                    data_tmp['correlations'] = cor_x[r,c];
                    data_tmp['pls_model'] = pls_model_I;
                    data_tmp['pls_method'] = method;
                    data_tmp['pls_options'] = {'pls_scale':scale,
                                                'lower':lower,
                                                'upper':upper,
                                                'weights':weights,
                        };
                    #data_tmp['experiment_id'] = experiment_ids[cnt];
                    #data_tmp['time_point'] = time_points[cnt];
                    data_loadings.append(data_tmp);
                    cnt+=1;
            # extract out loadings
            data_loadings_factors = [];
            cnt=0;
            for r in range(loadings_y.shape[0]):
                for c in range(loadings_y.shape[1]): #comp
                    data_tmp = {};
                    data_tmp['response_name'] = factors_unique[r];
                    data_tmp['loadings_response'] = loadings_y[r,c]
                    data_tmp['axis'] = c+1;
                    data_tmp['correlations_response'] = cor_y[r,c];
                    data_tmp['pls_model'] = pls_model_I;
                    data_tmp['pls_method'] = method;
                    data_tmp['pls_options'] = {'pls_scale':scale,
                                                'lower':lower,
                                                'upper':upper,
                                                'weights':weights,
                        };
                    #data_tmp['experiment_id'] = experiment_ids[cnt];
                    #data_tmp['time_point'] = time_points[cnt];
                    data_loadings_factors.append(data_tmp);
                    cnt+=1;
        except Exception as e:
            print(e);
            exit(-1);
        return data_scores,data_loadings,data_loadings_factors;
    def extract_mvr_scoresAndLoadings(self,
            mvr_model_I,
            ):
        '''extract out mvr scores and loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('%s' %(mvr_model_I));
            ans = robjects.r(r_statement);
            # get the canonical.correlations
            #canonical_correlations = np.array(ans.rx2('canonical.correlations'));#dim 1 = features, dim 2 = comp
            # get the fit residuals
            residuals = np.array(ans.rx2("residuals")); #dim1=features, dim2=comp, dim3=factors
            #get the scores
            scores_x = np.array(ans.rx2('scores')); #dim 1 = samples, dim 2 = comp
            scores_y = np.array(ans.rx2('Yscores'));
            #get the loadings
            loadings_x = np.array(ans.rx2('loadings')); #dim 1 = features, dim 2 = comp
            loadings_y = np.array(ans.rx2('Yloadings')); #dim 1 = factors, dim 2 = comp
            #get the means
            means_x = np.array(ans.rx2('Xmeans')); #dim 1 = features, dim 2 = comp
            means_y = np.array(ans.rx2('Ymeans')); #dim 1 = factors, dim 2 = comp
            # get the variance of each component
            var_x = np.array(ans.rx2("Xvar")); #dim 1 = comp
            var_x_total = np.array(ans.rx2('Xtotvar'));
            # calculate the correlation matrix
            cor_x = self.calculate_mvr_correlation(
                        mvr_model_I,
                        'correlation_m',
                        comps='1:'+str(loadings_x.shape[1]),
                        );
            cor_y = self.calculate_mvr_correlationResponse(
                        mvr_model_I,
                        'correlation_response_m',
                        comps='1:'+str(loadings_y.shape[1]),
                        );
            # get the explained variance
            var_proportion, var_cumulative = self.calculate_mvr_explainedVariance(
                    mvr_model_I,);
            # extract out scores
            return scores_x,scores_y,var_proportion,var_cumulative,loadings_x,loadings_y,cor_x,cor_y;
        except Exception as e:
            print(e);
            exit(-1);
    
    def calculate_mvr_performance(self,
            mvr_model_I,
            pls_model_I,
            ncomp,
            Y_add,
            scale,
            validation,
            segments,
            method,
            stripped,
            lower,
            upper, 
            trunc_pow, 
            weights,
            p_method,
            nperm):
        '''extract out mvr scores and loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_perf
        '''
        try:
            msep_reduced,rmsep_reduced,r2_reduced,q2_reduced,r2x_reduced=self.extract_mvr_performance(
                    mvr_model_I,
                    validation,)
            # extract out the model performance statistics
            data_perf = [];
            for i in range(len(msep_reduced)): #model
                data_tmp = {};
                data_tmp['pls_model'] = pls_model_I;
                data_tmp['pls_method'] = method;
                data_tmp['pls_options'] = {'pls_scale':scale,
                                            'lower':lower,
                                            'upper':upper,
                                            'weights':weights,
                    };
                data_tmp['pls_msep'] = msep_reduced[i];
                data_tmp['pls_rmsep'] = rmsep_reduced[i];
                data_tmp['pls_r2'] = r2_reduced[i]
                data_tmp['pls_q2'] = q2_reduced[i];
                data_tmp['pls_r2x'] = r2x_reduced[i];
                data_tmp['crossValidation_ncomp'] = i;
                data_tmp['crossValidation_method'] = validation;
                data_tmp['crossValidation_options'] = {'segments':segments,
                    };
                data_tmp['permutation_nperm']=nperm;
                data_tmp['permutation_pvalue']=None;
                #data_tmp['permutation_pvalue_corrected']=pvalue_corrected[i];
                #data_tmp['permutation_pvalue_corrected_description']=pvalue_method[i];
                data_tmp['permutation_pvalue_corrected']=None;
                data_tmp['permutation_pvalue_corrected_description']=None;
                data_tmp['permutation_options']={};
                data_perf.append(data_tmp);
        except Exception as e:
            print(e);
            exit(-1);
        return data_perf;
    def extract_mvr_performance(self,
            mvr_model_I,
            validation,):
        '''extract out mvr performance
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_perf
        '''
        try:
            r_statement = ('%s' %(mvr_model_I));
            ans = robjects.r(r_statement);
            # get model validation
            validation_results = ans.rx2("validation");
            press = np.array(validation_results.rx2('PRESS'));
            adj = np.array(validation_results.rx2('adj'));
            # get the model cross validation statistics
            # dim 1 = "CV", "adjCV"
            # dim 2 = response variables
            #           e.g. for 6 factors/response variables, there should be 6 rows
            # dim 3 = models from ncomps = no comps, ncomps[1]:ncomps[length(ncomps)]
            #           e.g. for 5 comps there should be 6 models
            r_statement = ('RMSEP(%s)'%(mvr_model_I)); # root mean squared error
            ans = robjects.r(r_statement);
            rmsep = np.array(ans.rx2('val'))[1]; #adjCV
            rmsep_reduced = np.zeros_like(rmsep[0,:]);
            for j in range(rmsep.shape[1]):
                rmsep_reduced[j]=rmsep[:,j].sum();
            r_statement = ('MSEP(%s)'%(mvr_model_I));
            ans = robjects.r(r_statement);
            msep = np.array(ans.rx2('val'))[1]; #adjCV
            msep_reduced = np.zeros_like(msep[0,:]);
            for j in range(msep.shape[1]):
                msep_reduced[j]=msep[:,j].sum();
            #calculate the Q2
            r_statement = ('mvrValstats(%s, estimate="%s")'%(mvr_model_I,validation));
            ans = robjects.r(r_statement);
            sse = np.array(ans.rx2('SSE'))[0];
            sst = np.array(ans.rx2('SST'))[0];
            #q2 = 1-(sse/sst);
            sse_reduced = np.zeros_like(sse[0,:]);
            sst_reduced = np.zeros_like(sse[0,:]);
            for j in range(sse.shape[1]):
                sse_reduced[j]=sse[:,j].sum();
                sst_reduced[j]=sse.shape[0]*sst[0];
            q2_reduced = 1-(sse_reduced/sst_reduced);
            #calculate the R2
            r_statement = ('mvrValstats(%s, estimate="train")'%(mvr_model_I));
            ans = robjects.r(r_statement);
            sse = np.array(ans.rx2('SSE'))[0];
            sst = np.array(ans.rx2('SST'))[0];
            #r2 = 1-(sse/sst);
            sse_reduced = np.zeros_like(sse[0,:]);
            sst_reduced = np.zeros_like(sse[0,:]);
            for j in range(sse.shape[1]):
                sse_reduced[j]=sse[:,j].sum();
                sst_reduced[j]=sse.shape[0]*sst[0];
            r2_reduced = 1-(sse_reduced/sst_reduced);
            #calculate the R2X
            var_proportion, var_cumulative = self.calculate_mvr_explainedVariance(
                    mvr_model_I,);
            r2x_reduced = np.insert(var_cumulative,0,0)
            # extract out the model performance statistics
            return msep_reduced,rmsep_reduced,r2_reduced,q2_reduced,r2x_reduced;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_mvr_permutations(self,
                                   ):
        '''perform permutation test on the model
        INPUT:
        OUTPUT:
        TODO:
        '''
        try:
            ## perform permutation test
            #r_statement = ('\
            #    test = MVA.test(%s, %s, cmv = %s,  ncomp = %s,\
            #    kout = %s, kinn = %s,\
            #    model = "%s",\
            #    scale = %s, validation = "%s", segments = %s,\
            #    lower = %s, upper = %s,\
            #    weights = %s,\
            #    nperm = %s)'\
            #            %('concentrations_mt','factors_f',"TRUE",
            #              ncomp,
            #              segments,segments,
            #              pls_model_I,
            #              scale,validation,segments,
            #              #method, , method = "%s"
            #              lower,upper,
            #              weights,
            #              nperm ));
            #ans = robjects.r(r_statement);
            #pvalue_corrected = np.array(ans.rx2('p.value'));
            #pvalue_method = np.array(ans.rx2('p.adjust.method'));
            # extract out the model performance statistics
            data_perf = [];
            for i in range(len(msep_reduced)): #model
                data_tmp = {};
                data_tmp['pls_model'] = pls_model_I;
                data_tmp['pls_method'] = method;
                data_tmp['pls_options'] = {'pls_scale':scale,
                                            'lower':lower,
                                            'upper':upper,
                                            'weights':weights,
                    };
                data_tmp['pls_msep'] = msep_reduced[i];
                data_tmp['pls_rmsep'] = rmsep_reduced[i];
                data_tmp['pls_r2'] = r2_reduced[i]
                data_tmp['pls_q2'] = q2_reduced[i];
                data_tmp['crossValidation_ncomp'] = i;
                data_tmp['crossValidation_method'] = validation;
                data_tmp['crossValidation_options'] = {'segments':segments,
                    };
                data_tmp['permutation_nperm']=nperm;
                data_tmp['permutation_pvalue']=None;
                #data_tmp['permutation_pvalue_corrected']=pvalue_corrected[i];
                #data_tmp['permutation_pvalue_corrected_description']=pvalue_method[i];
                data_tmp['permutation_pvalue_corrected']=None;
                data_tmp['permutation_pvalue_corrected_description']=None;
                data_tmp['permutation_options']={};
                data_perf.append(data_tmp);
        except Exception as e:
            print(e);
            exit(-1);
        return data_perf;
    def calculate_mvr_coefficients(self,
            mvr_model_I,
            pls_model_I,
            factors_unique,
            cn_sorted,
            cgn,
            scale,
            lower,
            upper,
            weights,
            method,
            ):
        '''extract out the coefficients of the mvr model

        requires caret

        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_coefficients

        '''
        try:
            coefficients,coefficients_comps_reduced,coefficients_reduced=self.extract_mvr_coefficients(mvr_model_I);                    
            # extract out coefficients
            data_coefficients = [];
            for r in range(coefficients_comps_reduced.shape[0]):
                for c in range(coefficients_comps_reduced.shape[1]):
                    data_tmp = {};
                    data_tmp['response_name'] = factors_unique[c];
                    data_tmp['component_name'] = cn_sorted[r];
                    data_tmp['component_group_name'] = cgn[r];
                    data_tmp['pls_coefficients'] = coefficients_comps_reduced[r,c];
                    data_tmp['pls_model'] = pls_model_I;
                    data_tmp['pls_method'] = method;
                    data_tmp['pls_options'] = {'pls_scale':scale,
                                                'lower':lower,
                                                'upper':upper,
                                                'weights':weights,
                        };
                    data_coefficients.append(data_tmp);
            for r in range(coefficients_comps_reduced.shape[0]):
                data_tmp = {};
                data_tmp['response_name'] = 'all';
                data_tmp['component_name'] = cn_sorted[r];
                data_tmp['component_group_name'] = cgn[r];
                data_tmp['pls_coefficients'] = coefficients_reduced[r];
                data_tmp['pls_model'] = pls_model_I;
                data_tmp['pls_method'] = method;
                data_tmp['pls_options'] = {'pls_scale':scale,
                                            'lower':lower,
                                            'upper':upper,
                                            'weights':weights,
                    };
                data_coefficients.append(data_tmp);
        except Exception as e:
            print(e);
            exit(-1);
        return data_coefficients;
    def extract_mvr_coefficients(self,
            mvr_model_I,
            ):
        '''extract out the coefficients of the mvr model

        requires caret

        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_coefficients

        '''
        try:
            r_statement = ('%s' %(mvr_model_I));
            ans = robjects.r(r_statement);
            #get the coefficients
            coefficients = np.array(ans.rx2('coefficients')); #c(features, features, comp)
            coefficients_comps_reduced = np.zeros([coefficients.shape[0],coefficients.shape[1]]);
            coefficients_reduced = np.zeros([coefficients.shape[0]]);
            for i in range(coefficients.shape[0]):
                coefficients_reduced[i] = np.abs(coefficients[i,:,:]).sum();
                for j in range(coefficients.shape[1]):
                    coefficients_comps_reduced[i,j] = np.abs(coefficients[i,j,:]).sum();
                    
            # extract out coefficients
            return coefficients,coefficients_comps_reduced,coefficients_reduced;
        except Exception as e:
            print(e);
            exit(-1);
        return data_coefficients;

    def calculate_mvr_correlation(self,
            mvr_model_I,
            correlation_O,
            comps='1:2',
            method = "pearson"
            ):
        '''Calculate the correlation matrix
        INPUT:
        mvr_model_I
        ncomps
        OUTPUT:
        data_O = correlation matrix
        correlation_O = name of the R object containing the correlation matrix
        '''
        try:
            r_statement = ('%s <- cor(model.matrix(%s), %s$scores, method = "%s")' %(
                correlation_O,mvr_model_I,mvr_model_I,method))
            #should this be "%s$loadings"?
            ans = robjects.r(r_statement);
            data_O = np.array(ans); #dim 1 = features, dim 2 = comps
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_mvr_correlationResponse(self,
            mvr_model_I,
            correlation_O,
            comps='1:2',
            ):
        '''Calculate the correlation matrix
        INPUT:
        mvr_model_I
        ncomps
        OUTPUT:
        data_O = correlation matrix
        correlation_O = name of the R object containing the correlation matrix
        '''
        try:
            r_statement = ('%s <- cor(model.matrix(%s), %s$Yscores)' %(correlation_O,mvr_model_I,mvr_model_I,))
            ans = robjects.r(r_statement);
            data_O = np.array(ans); #dim 1 = features, dim 2 = comps
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_mvr_explainedVariance(self,
            mvr_model_I,):
        '''calculate the explained variance in X
        INPUT:
        OUTPUT:
        '''
        try:
            # get the explained variance
            r_statement = ('explvar(%s)'%(mvr_model_I)); #var_ex = (1-(var_x_total[0]-var_x)/var_x_total[0])*100
            ans = robjects.r(r_statement);
            var_proportion = np.array(ans); 
            var_proportion = var_proportion/100.0; #remove the percent
            var_cumulative = np.zeros_like(var_proportion);
            for i in range(len(var_proportion)):
                if i==0:
                    var_cumulative[i] = var_proportion[i];
                else:
                    var_cumulative[i] = var_proportion[i]+var_cumulative[i-1];
            return var_proportion,var_cumulative;
        except Exception as e:
            print(e);
            exit(-1);
    def extract_mvr_scores(self,
            mvr_model_I,
            ):
        '''extract out mvr scores
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('scores(%s)' %(mvr_model_I));
            ans = robjects.r(r_statement);
            scores_x = np.array(ans);
            return scores_x;
        except Exception as e:
            print(e);
            exit(-1);

    def extract_mvr_Yscores(self,
            mvr_model_I,
            ):
        '''extract out mvr Yscores
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('Yscores(%s)' %(mvr_model_I));
            ans = robjects.r(r_statement);
            scores_y = np.array(ans);
            return scores_y;
        except Exception as e:
            print(e);
            exit(-1);

    def extract_mvr_loadings(self,
            mvr_model_I,
            ):
        '''extract out mvr loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('loadings(%s)' %(mvr_model_I));
            ans = robjects.r(r_statement);
            loadings_x = np.array(ans);
            return loadings_x;
        except Exception as e:
            print(e);
            exit(-1);

    def extract_mvr_Yloadings(self,
            mvr_model_I,
            ):
        '''extract out mvr scores
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('Yloadings(%s)' %(mvr_model_I));
            ans = robjects.r(r_statement);
            loadings_y = np.array(ans);
            return loadings_y;
        except Exception as e:
            print(e);
            exit(-1);

            loadings_x = np.array(ans.rx2('loadings')); #dim 1 = features, dim 2 = comp
            loadings_y = np.array(ans.rx2('Yloadings')); #dim 1 = factors, dim 2 = comp
            #get the means
            means_x = np.array(ans.rx2('Xmeans')); #dim 1 = features, dim 2 = comp
            means_y = np.array(ans.rx2('Ymeans')); #dim 1 = factors, dim 2 = comp
            # get the variance of each component
            var_x = np.array(ans.rx2("Xvar")); #dim 1 = comp
            var_x_total = np.array(ans.rx2('Xtotvar'));
            # calculate the correlation matrix
            cor_x = self.calculate_mvr_correlation(
                        mvr_model_I,
                        'correlation_m',
                        comps='1:'+str(loadings_x.shape[1]),
                        );
            cor_y = self.calculate_mvr_correlationResponse(
                        mvr_model_I,
                        'correlation_response_m',
                        comps='1:'+str(loadings_y.shape[1]),
                        );
            # get the explained variance
            var_proportion, var_cumulative = self.calculate_mvr_explainedVariance(
                    mvr_model_I,);
            # extract out scores
            return scores_x,scores_y,var_proportion,var_cumulative,loadings_x,loadings_y,cor_x,cor_y;
        except Exception as e:
            print(e);
            exit(-1);            

    def extract_mvr_coefficientsMeanAndVar(self,
            mvr_model_I,
            ):
        '''extract out mvr scores
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        NOTES:
        requires a cross validated model
        '''
        try:
            r_statement = ('%s' %(mvr_model_I));
            ans = robjects.r(r_statement);
            #get the means
            means_x = np.array(ans.rx2('Xmeans')); #dim 1 = features, dim 2 = comp
            means_y = np.array(ans.rx2('Ymeans')); #dim 1 = factors, dim 2 = comp
            # get the variance of each component
            var_x = np.array(ans.rx2("Xvar")); #dim 1 = comp
            var_x_total = np.array(ans.rx2('Xtotvar'));
            return means_x,means_y,var_x,var_x_total;
        except Exception as e:
            print(e);
            exit(-1);
