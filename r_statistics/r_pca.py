from .r_dependencies import *
from .r_correlation import r_correlation
class r_pca(r_correlation):
    def calculate_pca_princomp(self,data_I,cor_I = "FALSE", scores_I = "TRUE", covmat_I="NULL", na_action_I='na.omit',robust_pca_I=False,
            center_I = "TRUE",
            scale_I="TRUE"):
        '''PCA analysis using princomp from R
        Note: The calculation is done using eigen on the correlation or covariance matrix, as determined by cor

        princomp(x, cor = FALSE, scores = TRUE, covmat = NULL,
         subset = rep_len(TRUE, nrow(as.matrix(x))), ...)
        
        INPUT:
        data	
        an optional data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
        subset	
        an optional vector used to select rows (observations) of the data matrix x.
        na.action	
        a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The factory-fresh default is na.omit.
        x	
        a numeric matrix or data frame which provides the data for the principal components analysis.
        cor	
        a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.)
        scores	
        a logical value indicating whether the score on each principal component should be calculated.
        covmat	
        a covariance matrix, or a covariance list as returned by cov.wt (and cov.mve or cov.mcd from package MASS). If supplied, this is used rather than the covariance matrix of x.
        ...	
        arguments passed to or from other methods. If x is a formula one might specify cor or scores.
        
        perform robust PCA
        pc.rob <- princomp(X, covmat= MASS::cov.rob(X))
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
            column_variables_I=['sample_name_abbreviation'],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        sna = column_variables['sample_name_abbreviation'];
        sna_unique = listdict.get_uniqueValues_list(sna);
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
                # run PCA
                self.princomp_pca(
                        'concentrations_mt',
                        'result',
                        robust_pca_I,
                        na_action_I,cor_I, scores_I,covmat_I,
                        center_I,
                        scale_I);
                # get the scores and loadings
                data_scores,data_loadings = self.extract_pcaMethods_scoresAndLoadings(
                    'result',
                    'concentrations_mt',
                    pca_model_I,
                    pca_method_I,
                    sns_sorted,
                    sna,sna_unique,
                    cn_sorted,
                    cgn,
                    scale,
                    center,
                    );
                # get the validation
                data_perf = [];
            except Exception as e:
                print(e);
                exit(-1);
            return data_scores,data_loadings,data_perf;
        else:
            print('missing values found!');
    def calculate_pca_prcomp(self,data_I,
        retx_I = "TRUE", center_I = "TRUE", na_action_I='na.omit',scale_I="TRUE"):
        '''PCA analysis using prcomp from R
        Note: calculations are made using svd'''
        
        # format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList_pd(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name'],
            column_variables_I=['sample_name_abbreviation'],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        sna = column_variables['sample_name_abbreviation'];
        sna_unique = listdict.get_uniqueValues_list(sna);
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
                # calls for pca analysis
                r_statement = ('pc = prcomp(t(concentrations_m), na.action=%s, retx=%s, center=%s, scale=%s)' %(na_action_I,retx_I,center_I,scale_I));
                #need to send in the transpose
                ans = robjects.r(r_statement);
                r_statement = ('summary(pc)')
                ans = robjects.r(r_statement);
                importance = np.array(ans.rx2('importance'))
                stdev = np.array(ans.rx2('sdev')) # or importance[0,:]
                var_proportion = importance[1,:]
                var_cumulative = importance[2,:]
                loadings = np.array(ans.rx2('rotation')) #nrow = mets x ncol = loadings axis
                scores = np.array(ans.rx2('x')) #nrow = samples x ncol = pc axis
                # extract out scores
                data_scores = [];
                cnt=0;
                for r in range(scores.shape[0]):
                    for c in range(scores.shape[1]):
                        data_tmp = {};
                        data_tmp['sample_name_short'] = sns_sorted[r];
                        data_tmp['score'] = scores[r,c];
                        data_tmp['axis'] = c+1;
                        data_tmp['var_proportion'] = var_proportion[c];
                        data_tmp['var_cumulative'] = var_cumulative[c];
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_scores.append(data_tmp);
                        cnt+=1;
                # extract out loadings
                data_loadings = [];
                cnt=0;
                for r in range(loadings.shape[0]):
                    for c in range(loadings.shape[1]):
                        data_tmp = {};
                        data_tmp['component_name'] = cn_sorted[r]; #need to double check
                        data_tmp['component_group_name'] = cgn[r];
                        data_tmp['loadings'] = loadings[r,c]; #need to double check
                        data_tmp['axis'] = c+1;
                        #data_tmp['experiment_id'] = experiment_ids[cnt];
                        #data_tmp['time_point'] = time_points[cnt];
                        data_loadings.append(data_tmp);
                        cnt+=1;
            except Exception as e:
                print(e);
                exit(-1);
            return data_scores,data_loadings
        else:
            print('missing values found!')
    def calculate_pca_pcaMethods(self,
            data_I,
            pca_model_I="pca",
            pca_method_I="svd",
            imputeMissingValues="TRUE",
            cv="q2",
            ncomps="7",
            scale="uv",
            center="TRUE",
            segments="10",
            nruncv="1",
            crossValidation_type = "krzanowski",
            return_data_imputed_I=False,
            geometric_imputation_I=False):
        '''Calculate pca using various methods from pcaMethods
        Bioinformatics (2007) 23 (9): 1164-1167.
        doi: 10.1093/bioinformatics/btm069
        INPUT:
        data_I = listDict
        pca_model_I = string
        pca_method_I = string, pca method
        OPTIONAL INPUT for pcaMethods:
        imputeMissingValues
        ncomps = nPcs
        scale
        center
        completeObs
        
        OPTIONAL INPUT for missing value imputation
        return_data_imputed_I = boolean, return only imputed data
        geometric_imputation_I = boolean, log normalize values before PCA
            required to ensure only positive values are returned

        OUTPUT:
        scores "matrix", the calculated scores
        loadings "matrix", the calculated loadings
        R2cum "numeric", the cumulative R2 values
        sDev "numeric", the individual standard deviations of the score vectors
        R2 "numeric", the individual R2 values
        cvstat "numeric", cross-validation statistics
        nObs "numeric", number of observations
        nVar "numeric", number of variables
        centered "logical", data was centered or not
        center "numeric", the original variable centers
        scaled "logical", data was scaled or not
        scl "numeric", the original variable scales
        varLimit "numeric", the exceeded variance limit
        nPcs,nP "numeric", the number of calculated PCs
        method "character", the method used to perform PCA
        missing "numeric", the total amount of missing values in original data
        completeObs "matrix", the estimated complete observations
        network "nlpcaNet", the network used by non-linear PCA

        DESCRIPTION:
        pca(object, method=listPcaMethods(), nPcs=2, scale=c("none",
            "pareto", "vector", "uv"), center=TRUE, completeObs=TRUE,
            subset, cv=c("none", "q2"), ...)

        svd: Uses classical prcomp.See documentation for svdPca.
        nipals: An iterative method capable of handling small amounts of missing values. See documentation
        for nipalsPca.
        rnipals: Same as nipals but implemented in R.
        bpca: An iterative method using a Bayesian model to handle missing values. See documentation
        for bpca.
        ppca: An iterative method using a probabilistic model to handle missing values. See documentation
        for ppca.
        svdImpute: Uses expectation maximation to perform SVD PCA on incomplete data. See documentation
        for svdImpute.

        Other methods
        robustPca(Matrix, nPcs=2, verbose=interactive(), ...)

        nlpca(Matrix, nPcs=2, maxSteps=2 * prod(dim(Matrix)), unitsPerLayer,
            functionsPerLayer, weightDecay=0.001, weights,
            verbose=interactive(), ...)


        robustSvd(x)
        computes x = u d v'.
            output: u, d, v

        #scale/center the data
        mdc <- scale(metaboliteDataComplete, center=TRUE, scale=FALSE)
        #calculate the pca using one of several methods
        resPCA <- pca(mdC, method="svd", center=FALSE, nPcs=5)
        resPPCA <- pca(md, method="ppca", center=FALSE, nPcs=5)
        resBPCA <- pca(md, method="bpca", center=FALSE, nPcs=5)
        resSVDI <- pca(md, method="svdImpute", center=FALSE, nPcs=5)
        resNipals <- pca(md, method="nipals", center=FALSE, nPcs=5)
        resNLPCA <- pca(md, method="nlpca", center=FALSE, nPcs=5, maxSteps=300)

        q2SVDI <- Q2(resSVDI, mdC, fold=10)
        R2cum

        errEsti <- kEstimate(md, method = "ppca", evalPcs=1:5, nruncv=1, em="nrmsep")
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
            column_variables_I=['sample_name_abbreviation','experiment_id','time_point'],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        sna = column_variables['sample_name_abbreviation'];
        sna_unique = listdict.get_uniqueValues_list(sna);
        # check if there were any missing values in the data set in the first place
        mv = 0;
        mv = listdict.count_missingValues_pivotTable();
        if mv==0 or (mv!=0 and pca_method_I in ["svdImpute",'ppca','nipals','bpca']):
            # Call to R
            try:
                # clear the workspace
                self.clear_workspace();
                # convert lists to R matrix
                self.make_matrixFromList(concentrations,len(cn_sorted),len(sns_sorted),'concentrations_m');
                # convert to the transpose
                self.transpose_matrix('concentrations_mt','concentrations_m');
                # OPTIONAL: impute in geometric space
                if geometric_imputation_I:
                    r_statement = 'concentrations_mt = log(concentrations_mt)'
                    ans = robjects.r(r_statement);
                #scale/center the data
                #self.pcaMethods_scale('concentrations_mt',
                #        'concentrations_scaled',
                #        center,
                #        scale);
                # call pcaMethods #robjects.r('slplot(result)')
                self.pcaMethods_pca(
                        'concentrations_mt',
                        'result',
                        pca_method_I,
                        ncomps,
                        imputeMissingValues,
                        center,
                        scale,
                        cv,
                        segments,
                        nruncv,
                        crossValidation_type,);
                # get the scores and loadings
                data_scores,data_loadings = self.extract_pcaMethods_scoresAndLoadings(
                    'result',
                    'concentrations_mt',
                    pca_model_I,
                    pca_method_I,
                    sns_sorted,
                    sna,sna_unique,
                    cn_sorted,
                    cgn,
                    scale,
                    center,
                    );
                # get the validation
                if cv != "none":
                    data_perf = self.extract_pcaMethods_crossValidation(
                        'result',
                        'concentrations_mt',
                        pca_model_I,
                        pca_method_I,
                        sns_sorted,
                        sna,sna_unique,
                        cn_sorted,
                        cgn,
                        scale,
                        center,
                        ncomps,
                        cv,
                        segments,
                        nruncv,
                        crossValidation_type,
                        );
                    return data_scores,data_loadings,data_perf;
                if return_data_imputed_I:
                    # get the full data set
                    data_imputed = self.extract_pcaMethods_missingValues_listDict(
                        pcamethods_var_I='result',
                        imputed_O = 'concentrations_imputed',
                        row_variables_I = row_variables,
                        column_variables_I = column_variables,
                        geometric_imputation_I = geometric_imputation_I
                        );
                    return data_imputed;
                else:
                    return data_scores,data_loadings;
            except Exception as e:
                print(e);
                exit(-1);
        else:
            print('missing values found!');
    def pcaMethods_scale(self,
            data_var_I,
            data_var_O,
            center,
            scale,
            ):
        '''Scale and/or center the data
        INPUT:
        data_var_I = name of the R workspace variable to scale/center
        center = string, default, TRUE
        scale = string, default, FALSE
                c("none","pareto", "vector", "uv")
        OUTPUT:
        data_var_O = name of the R object containing the scaled/centered data
        '''
        try:
            r_statement = ('%s <- prep(%s, center=%s, scale="%s")' %(
                data_var_O,data_var_I,center,scale,))
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def pcaMethods_pca(self,
            data_var_I,
            data_var_O,
            pca_method_I,
            ncomps,
            imputeMissingValues,
            #prep arguments
            center="FALSE",
            scale="uv",
            #Q2 arguments
            cv="q2",
            segments="10",
            nruncv="1",
            type = "krzanowski",
            ):
        '''Call pcaMethods pca routine
        INPUT:
        data_var_I = name of the R workspace variable to calculate the pca
        center = string, default, TRUE
        scale = string, default, FALSE
                c("none","pareto", "vector", "uv")
        imputeMissingValues = string, default, TRUE
        cv=c("none", "q2")
        segments = The number of groups to divide the data in.
        nruncv = The number of times to repeat the whole cross-validation
        crossValidation_type = krzanowski or imputation type cross-validation


        OUTPUT:
        data_var_O = name of the resulting R object 
        '''
        try:
            #r_statement = ('%s <- pca(%s, method="%s", center=%s, scale="%s", nPcs=%s, completeObs=%s, cv="%s",\
            #    fold=%s, nruncv=%s, type="%s")' %(
            #    data_var_O,data_var_I,pca_method_I,center,scale,ncomps,imputeMissingValues,cv,
            #    segments,nruncv,type)); #BUG: specifying nruncv results in an error
            r_statement = ('%s <- pca(%s, method="%s", center=%s, scale="%s", nPcs=%s, completeObs=%s, cv="%s",\
                fold=%s, type="%s")' %(
                data_var_O,data_var_I,pca_method_I,center,scale,ncomps,imputeMissingValues,cv,
                segments,type));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def pcaMethod_kEstimate(self,
            data_var_I,
            output_var_O,
            pca_method_I,
            ncomps,
            em="nrmsep",
            segments = "10",
            nruncv="5",
            allVariables = "TRUE"
            ):
        '''errEsti <- kEstimate(md, method = "ppca", evalPcs=1:5, nruncv=1, em="nrmsep")
        INPUT:
        data_var_I = name of the R workspace variable
        pca_method_I = name of the pcaMethod
        ncomps = number of components to use in the validation
        segments = float string, number of segments for cross validation
        nruncv = float string, times the whole cv is repeated
        em = error measure, "nrmsep" or "q2"
        OUTPUT:
        output_var_O = name of the R workspace variable

        bestNPcs number of PCs or k for which the minimal average NRMSEP or the maximal
                Q2 was obtained.
        eError an array of of size length(evalPcs). Contains the average error of the cross validation
            runs for each number of components.
        variableWiseError
            Matrix of size incomplete_variables x length(evalPcs). Contains the NRMSEP
            or Q2 distance for each variable and each number of PCs. This allows to
            easily see for wich variables imputation makes sense and for which one it should
            not be done or mean imputation should be used.
        evalPcs The evaluated numbers of components or number of neighbours (the same as
            the evalPcs input parameter).
        variableIx Index of the incomplete variables. This can be used to map the variable wise
            error to the original data.
        '''
        try:
            r_statement = ('%s <- kEstimate(%s, method="%s", evalPcs=1:%s, segs=%s, em="%s", nruncv=%s, allVariables = TRUE)' %(
                output_var_O,data_var_I,pca_method_I,ncomps,segments,em,nruncv));
            ans = robjects.r(r_statement);
            error = np.array(ans.rx2('eError')) #estimated error NRMSEP
            verror = np.array(ans.rx2('variableWiseError'))
            r_statement = ('print(%s$bestNPcs)'%output_var_O)
            ans = robjects.r(r_statement);
            return error;
        except Exception as e:
            print(e);
            exit(-1);
    def extract_pcaMethods_scoresAndLoadings(self,
            pcaMethod_model_I,
            data_var_I,
            pca_model_I,
            pca_method_I,
            sns_sorted,
            sna,sna_unique,
            cn_sorted,
            cgn,
            scale,
            center,
            ):
        '''extract out mvr scores and loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        data_scores
        data_loadings
        '''
        try:
            r_statement = ('%s' %(pcaMethod_model_I));
            #r_statement = ('summary(%s)' %(pcaMethod_model_I));
            ans = robjects.r(r_statement);
            # get the loadings
            r_statement = ('result_loadings <- loadings(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            loadings_x = np.array(ans); #dim 1 = features, dim 2 = comps
            #get the scores
            r_statement = ('result_scores <- scores(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            scores_x = np.array(ans); #dim 1 = samples, dim 2 = comp
            #get the scores stdev
            var_proportion,var_cumulative = self.calculate_pcaMethods_explainedVariance(pcaMethod_model_I);
            #TODO
            ## calculate the correlation matrix
            #cor_x = self.calculate_mvr_correlation(
            #            pcaMethod_model_I,
            #            data_var_I,
            #            comps='1:'+str(loadings_x.shape[1]),
            #            );
            # extract out scores
            data_scores = [];
            cnt=0;
            for r in range(scores_x.shape[0]):
                for c in range(scores_x.shape[1]):
                    data_tmp = {};
                    data_tmp['sample_name_short'] = sns_sorted[r];
                    data_tmp['sample_name_abbreviation'] = sna[r];
                    data_tmp['score'] = scores_x[r,c];
                    data_tmp['axis'] = c+1;
                    data_tmp['var_proportion'] = var_proportion[c];
                    data_tmp['var_cumulative'] = var_cumulative[c];
                    data_tmp['pca_model'] = pca_model_I;
                    data_tmp['pca_method'] = pca_method_I;
                    data_tmp['pca_options'] = {'pca_scale':scale,
                                               'pca_center':center,
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
                    data_tmp['correlations'] = None; #TODO: cor_x[r,c];
                    data_tmp['pca_model'] = pca_model_I;
                    data_tmp['pca_method'] = pca_method_I;
                    data_tmp['pca_options'] = {'pca_scale':scale,
                                               'pca_center':center,
                        };
                    #data_tmp['experiment_id'] = experiment_ids[cnt];
                    #data_tmp['time_point'] = time_points[cnt];
                    data_loadings.append(data_tmp);
                    cnt+=1;
        except Exception as e:
            print(e);
            exit(-1);
        return data_scores,data_loadings;

    def extract_pcaMethods_R2cum(self,pcaMethod_model_I):
        ''' '''
        try:
            #get the r2 values
            r_statement = ('result_R2cum <- R2cum(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            r2 = np.array(ans); #dim 1 = samples, dim 2 = comp
            return r2;
        except Exception as e:
            print(e);
            exit(-1);

    def extract_pcaMethods_cvstat(self,pcaMethod_model_I):
        ''' '''
        try:
            # get the q2
            r_statement = ('result_cvstat <- cvstat(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            q2 = np.array(ans);
            return q2;
        except Exception as e:
            print(e);
            exit(-1);


    def extract_pcaMethods_crossValidation(self,
            pcaMethod_model_I,
            data_var_I,
            pca_model_I,
            pca_method_I,
            sns_sorted,
            sna,sna_unique,
            cn_sorted,
            cgn,
            scale,
            center,
            ncomps,
            cv="q2",
            segments="10",
            nruncv="1",
            type = "impute",
            ):
        '''extract out mvr scores and loadings
        INPUT:
        mvr_model_I = name of the R mvr model workspace variable
        OUTPUT:
        r2,q2,nrmsep,msep
        '''
        try:
            r_statement = ('%s' %(pcaMethod_model_I));
            ans = robjects.r(r_statement);
            #get the r2 values
            r_statement = ('result_R2cum <- R2cum(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            r2 = np.array(ans); #dim 1 = samples, dim 2 = comp
            # get the q2
            r_statement = ('result_cvstat <- cvstat(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            q2 = np.array(ans);
            # get the error estimates
            nrmsep = self.pcaMethod_kEstimate(
                    data_var_I,
                    'result_errors',
                    pca_method_I,
                    ncomps,
                    "nrmsep",
                    segments,
                    nruncv);
            msep = np.square(nrmsep);
            
            data_perf = [];
            for i in range(len(q2)): #model
                data_tmp = {};
                data_tmp['pca_model'] = pca_model_I;
                data_tmp['pca_method'] = pca_method_I;
                data_tmp['pca_options'] = {'pca_scale':scale,
                                               'pca_center':center,
                    };
                data_tmp['pca_msep'] = msep[i];
                data_tmp['pca_rmsep'] = nrmsep[i];
                data_tmp['pca_r2'] = r2[i]
                data_tmp['pca_q2'] = q2[i];
                data_tmp['crossValidation_ncomp'] = i;
                data_tmp['crossValidation_method'] = "CV";
                data_tmp['crossValidation_options'] = {
                    'segments':segments,
                    'nruncv':nruncv,
                    'type':type,
                    };
                data_tmp['permutation_nperm']=None;
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

    def calculate_pcaMethods_explainedVariance(self,
            pcaMethod_model_I,):
        '''calculate the explained variance in X
        INPUT:
        OUTPUT:
        '''
        try:
            # get the explained variance
            #get the scores stdev
            r_statement = ('result_sDev <- sDev(%s)' %(pcaMethod_model_I))
            ans = robjects.r(r_statement);
            scores_x_stdev = np.array(ans); #dim 1 = comp
            scores_x_var = np.square(scores_x_stdev)
            #scores_x_stdev_total = scores_x_stdev.sum();
            scores_x_var_total = scores_x_var.sum();
            ##var_ex = (1-(var_x_total[0]-var_x)/var_x_total[0])*100
            var_proportion = np.zeros_like(scores_x_stdev)
            var_cumulative = np.zeros_like(scores_x_stdev);
            for i in range(len(var_proportion)):
                #var_ex = 1-(scores_x_stdev_total-scores_x_stdev[i])/scores_x_stdev_total;
                var_ex = 1-(scores_x_var_total-scores_x_var[i])/scores_x_var_total;
                if i==0:
                    var_proportion[i] = var_ex;
                    var_cumulative[i] = var_ex;
                else:
                    var_proportion[i] = var_ex;
                    var_cumulative[i] = var_ex+var_cumulative[i-1];
            return var_proportion,var_cumulative;
        except Exception as e:
            print(e);
            exit(-1);
            
    def prcomp_pca(self,
            data_var_I,
            data_var_O,
            retx_I = "TRUE",
            center_I = "TRUE",
            na_action_I='na.omit',
            scale_I="TRUE"
            ):
        '''Calculate the correlation matrix
        INPUT:
        data_var_I = name of the R workspace variable to scale/center
        OUTPUT:
        data_var_O = name of the R object containing the scaled/centered data
        '''
        try:
            # calls for pca analysis
            r_statement = ('%s = prcomp(%s, na.action=%s, retx=%s, center=%s, scale=%s)' %(
                data_var_O,data_var_I,na_action_I,retx_I,center_I,scale_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
            
    def princomp_pca(self,
            data_var_I,
            data_var_O,
            robust_I=False,
            cor_I = "FALSE",
            scores_I = "TRUE",
            covmat_I="NULL",
            na_action_I='na.omit',
            center_I = "TRUE",
            scale_I="TRUE"
            ):
        '''Calculate the correlation matrix
        INPUT:
        data_var_I = name of the R workspace variable to scale/center
        OUTPUT:
        data_var_O = name of the R object containing the scaled/centered data

        perform robust PCA
        pc.rob <- princomp(X, covmat= MASS::cov.rob(X))

        Error in MASS::cov.rob(concentrations_mt) : at least 91 cases are needed
        '''
        try:
            # calls for pca analysis
            covmat = covmat_I;
            if robust_I:
                covmat = ('MASS::cov.rob(%s)' %(data_var_I));
            r_statement = ('%s = princomp(%s, na.action=%s, cor=%s, scores=%s, covmat=%s, center=%s, scale=%s)' %(
                data_var_O,data_var_I,na_action_I,cor_I, scores_I,covmat,center_I,scale_I));
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
        data_O = np array of the full dataset
        '''
        try:
            r_statement = ('%s <- completeObs(%s)' %(imputed_O,pcamethods_var_I));
            ans = robjects.r(r_statement);
            data_O = np.array(ans);
            return data_O;

        except Exception as e:
            print(e);
            exit(-1);

    def extract_pcaMethods_missingValues_listDict(self,
                pcamethods_var_I,
                imputed_O,
                row_variables_I,
                column_variables_I,
                geometric_imputation_I=False
                ):
        '''Extract out imputed missing values from pcaMethods as a listDict
        '''
        data = self.extract_pcaMethods_missingValues(
                pcamethods_var_I,
                imputed_O);
        if geometric_imputation_I: data = np.exp(data);
        
        #try:
        #    r_statement = ('%s' %('concentrations_mt'));
        #    ans = robjects.r(r_statement);
        #    concentrations = np.array(ans);
        #    dif = data-concentrations
        #    dif_min = min(abs(dif));
        #    dif_max = max(abs(dif));
        #except Exception as e:
        #    print(e);
        #    exit(-1);

        data_listDict = [];
        for r in range(data.shape[0]):
            for c in range(data.shape[1]):
                data_tmp = {};
                data_tmp['sample_name_short'] = column_variables_I['sample_name_short'][r];
                data_tmp['experiment_id'] = column_variables_I['experiment_id'][r]; #
                data_tmp['time_point'] = column_variables_I['time_point'][r]; #
                data_tmp['component_name'] = row_variables_I['component_name'][c]; #
                data_tmp['component_group_name'] = row_variables_I['component_group_name'][c]; #
                data_tmp['calculated_concentration'] = data[r,c]; #
                data_tmp['sample_name_abbreviation'] = column_variables_I['sample_name_abbreviation'][r];
                data_listDict.append(data_tmp);
        return data_listDict;
    def extract_pcaMethods_scores(self,
            pcaMethods_model_I,
            ):
        '''extract out pcaMethods scores
        INPUT:
        mvr_model_I = name of the R pcaMethods model workspace variable
        OUTPUT:
        data_scores
        '''
        try:
            r_statement = ('pcaMethods::scores(%s)' %(pcaMethods_model_I));
            ans = robjects.r(r_statement);
            scores_x = np.array(ans);
            return scores_x;
        except Exception as e:
            print(e);
            exit(-1);
    def extract_pcaMethods_loadings(self,
            pcaMethods_model_I,
            ):
        '''extract out pcaMethods loadings
        INPUT:
        mvr_model_I = name of the R pcaMethods model workspace variable
        OUTPUT:
        data_loadings
        '''
        try:
            r_statement = ('pcaMethods::loadings(%s)' %(pcaMethods_model_I));
            ans = robjects.r(r_statement);
            loadings_x = np.array(ans);
            return loadings_x;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_pcaMethods_correlation(self,
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
        TODO: broken, cannot use model.matrix, need to find a substitute
        '''
        try:
            r_statement = ('%s <- cor(pcaMethods::completeObs(%s), pcaMethods::scores(%s), method = "%s")' %(
                correlation_O,mvr_model_I,mvr_model_I,method))
            #should this be "%s$loadings"?
            ans = robjects.r(r_statement);
            data_O = np.array(ans); #dim 1 = features, dim 2 = comps
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);