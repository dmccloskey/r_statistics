from .r_dependencies import *
from .r_base import r_base
class r_dataNormalization(r_base):
    def calculate_glogNormalization_v1(self,data_I):
        '''normalize the data using a glog transformation using LMGene
        https://www.bioconductor.org/packages/release/bioc/html/LMGene.html
        Citation: Rocke D, Lee GC, Tillinghast J, Durbin-Johnson B and Wu S (2013). LMGene: LMGene Software for Data Transformation and Identification of Differentially Expressed Genes in Gene Expression Arrays. R package version 2.26.0, http://dmrocke.ucdavis.edu/software.html.
        INPUT:
        data_I = listDict
        ...
        OUTPUT:
        data_O = listDict of the transformed data
        concentrations = original data matrix
        concentrations_glog = normalized data matrix
        TODO:
        1. break into individual functions and calls to R
        2. add in optional input for calls to tranest()
        '''

        #make the ExpressionSet
        
        #format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'

        sns = []
        cn = []
        #replicates = [];
        sample_name_abbreviations = [];
        for d in data_I:
                sns.append(d['sample_name_short']);    
                #replicates.append(d['sample_replicate']);  
                sample_name_abbreviations.append(d['sample_name_abbreviation']) 
                cn.append(d['component_name']);
        sns_sorted = sorted(set(sns))
        #replicates_sorted = sorted(set(replicates))
        cn_sorted = sorted(set(cn))
        sample_name_abbreviations_sorted = sorted(set(sample_name_abbreviations))

        # extract out replicates
        replicates_dict = {};
        for sns in sns_sorted:
            replicates_dict[sns]=None;
        cnt_reps = 0;
        for sna_sorted in sample_name_abbreviations_sorted:
            for sns in sns_sorted:
                for d in data_I:
                    if d['sample_name_short'] == sns and d['sample_name_abbreviation'] == sna_sorted:
                        replicates_dict[sns] = cnt_reps;
                        cnt_reps+=1;
                        break;
            cnt_reps = 0;

        concentrations = ['NA' for r in range(len(sns_sorted)*len(cn_sorted))];
        experiment_ids = ['' for r in range(len(sns_sorted)*len(cn_sorted))];
        time_points = ['' for r in range(len(sns_sorted)*len(cn_sorted))];
        component_group_names = ['' for r in range(len(sns_sorted)*len(cn_sorted))];
        analysis_ids = ['' for r in range(len(sns_sorted)*len(cn_sorted))];
        calculated_concentration_units = ['' for r in range(len(sns_sorted)*len(cn_sorted))];
        cnt = 0;
        cnt_bool = True;
        cnt_reps = 0;
        sna = []
        replicates = []
        for c in cn_sorted:
                for s in sns_sorted:
                    for d in data_I:
                        if d['sample_name_short'] == s and d['component_name'] == c:
                            if d['calculated_concentration']:
                                concentrations[cnt] = d['calculated_concentration'];
                                experiment_ids[cnt] = d['experiment_id'];
                                time_points[cnt] = d['time_point'];
                                component_group_names[cnt] = d['component_group_name'];
                                analysis_ids[cnt] = d['analysis_id'];
                                calculated_concentration_units[cnt] = d['calculated_concentration_units'];
                                if cnt_bool:
                                    sna.append(d['sample_name_abbreviation']);
                                    replicates.append(replicates_dict[s]);
                                    #replicates.append(replicates_sorted[cnt_reps]);
                                    #if cnt_reps < len(replicates_sorted)-1:
                                    #    cnt_reps+=1;
                                    #else:
                                    #    cnt_reps=0;
                                break;
                    cnt = cnt+1
                cnt_bool = False;
        # check if there were any missing values in the data set in the first place
        mv = 0;
        for c in concentrations:
            if c=='NA':
                mv += 1;
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
                # convert lists to R list
                sna_r = '';
                for c in sna:
                    sna_r = (sna_r + ',' + '"' + c + '"');
                sna_r = sna_r[1:];
                replicates_r = '';
                for c in replicates:
                    replicates_r = (replicates_r + ',' + str(c));
                replicates_r = replicates_r[1:];
                r_statement = ('sna = c(%s)' % sna_r);
                ans = robjects.r(r_statement);
                r_statement = ('replicates = c(%s)' % replicates_r);
                ans = robjects.r(r_statement);
                r_statement = ('concentrations_l = list(sna=sna,replicates=replicates)');
                ans = robjects.r(r_statement);
                #convert to Expression Set
                r_statement = ('eS = neweS(concentrations_m,concentrations_l)');
                ans = robjects.r(r_statement);
                #estimate the parameters for g-log transformation
                #r_statement = ('tranpar = tranest(eS)');
                #r_statement = ('tranpar = tranest(eS, lowessnorm=TRUE)');
                #r_statement = ('tranpar = tranest(eS, mult=TRUE, lowessnorm=TRUE)');
                r_statement = ('tranpar = tranest(eS, mult=TRUE)'); # Matches metabo-analyst and produces the most uniform distribution
                ans = robjects.r(r_statement);
                r_statement = ('eS_transformed <- transeS(eS, tranpar$lambda, tranpar$alpha)');
                ans = robjects.r(r_statement);
                # extract out data matrices
                r_statement = ('exprs(eS_transformed)');
                ans = robjects.r(r_statement);
                concentrations_glog = numpy.array(ans);
                # convert array back to dict
                data_O = [];
                cnt = 0;
                for c in range(len(cn_sorted)):
                    for s in range(len(sns_sorted)):
                        if isinstance(concentrations_glog[c,s], (int, float, complex)):
                            data_tmp = {};
                            data_tmp['sample_name_short'] = sns_sorted[s]
                            data_tmp['component_name'] = cn_sorted[c]
                            data_tmp['component_group_name'] = component_group_names[cnt];
                            data_tmp['calculated_concentration'] = concentrations_glog[c,s];
                            data_tmp['experiment_id'] = experiment_ids[cnt];
                            data_tmp['time_point'] = time_points[cnt];
                            data_tmp['analysis_id'] = analysis_ids[cnt];
                            data_tmp['calculated_concentration_units'] = calculated_concentration_units[cnt]+ '_glog_normalized';
                            data_tmp['comment_'] = None;
                            data_tmp['used_'] = True;
                            data_O.append(data_tmp);
                            cnt+=1;
                        else:
                            print('concentration value is not a number.');
                #for c in range(len(sns_sorted)):
                #    for r in range(len(cgn_sorted)):
                        #if isinstance(concentrations_glog[r,c], (int, long, float, complex)):
                        #    data_tmp = {};
                        #    data_tmp['sample_name_short'] = sns_sorted[c]
                        #    data_tmp['component_name'] = cgn_sorted[r]
                        #    data_tmp['calculated_concentration'] = concentrations_glog[r,c];
                        #    #sns_O.append(sns_sorted[c]);
                        #    #cn_O.append(cgn_sorted[r]);
                        #    #cc_O.append(ans[c*len(cgn_sorted)+r]);
            except Exception as e:
                print(e);
                exit(-1);

            # reshape original concentrations
            concentrations_original = numpy.array(concentrations);
            concentrations = concentrations_original.reshape(len(cn_sorted),len(sns_sorted));
            return data_O, concentrations, concentrations_glog;
        else:
            print('missing values found in data!');

    def calculate_glogNormalization(self,data_I,
            mult="TRUE",
            lowessnorm="FALSE"
            ):
        '''normalize the data using a glog transformation using LMGene
        https://www.bioconductor.org/packages/release/bioc/html/LMGene.html
        Citation: Rocke D, Lee GC, Tillinghast J, Durbin-Johnson B and Wu S (2013). LMGene: LMGene Software for Data Transformation and Identification of Differentially Expressed Genes in Gene Expression Arrays. R package version 2.26.0, http://dmrocke.ucdavis.edu/software.html.
        INPUT:
        data_I = listDict
        ...
        OUTPUT:
        data_O = listDict of the transformed data
        concentrations = original data matrix
        concentrations_glog = normalized data matrix
        TODO:
        1. break into individual functions and calls to R
        2. add in optional input for calls to tranest()
        '''

        listdict = listDict(data_I);
        concentrations,cn_sorted,sns_sorted,row_variables,column_variables = listdict.convert_listDict2dataMatrixList(
            row_label_I='component_name',
            column_label_I='sample_name_short',
            value_label_I='calculated_concentration',
            row_variables_I=['component_group_name','calculated_concentration_units'],
            column_variables_I=['sample_name_abbreviation','experiment_id','time_point','analysis_id'],
            data_IO=[],
            na_str_I="NA");
        cgn = row_variables['component_group_name'];
        calculated_concentration_units = row_variables['calculated_concentration_units'];
        experiment_ids = column_variables['experiment_id'];
        time_points = column_variables['time_point'];
        analysis_ids = column_variables['analysis_id'];
        sna = column_variables['sample_name_abbreviation'];
        nsna_unique,sna_unique = listdict.get_uniqueValues('sample_name_abbreviation');

        #make replicate numbers for each sample abbreviation
        # extract out replicates
        replicates_dict = {};
        for sns in sns_sorted:
            replicates_dict[sns]=None;
        cnt_reps = 0;
        for sna_sorted in sna_unique:
            for sns in sns_sorted:
                for d in data_I:
                    if d['sample_name_short'] == sns and d['sample_name_abbreviation'] == sna_sorted:
                        replicates_dict[sns] = cnt_reps;
                        cnt_reps+=1;
                        break;
            cnt_reps = 0;
        replicates = [];
        for s in sns_sorted:
            replicates.append(replicates_dict[s]);

        # check if there were any missing values in the data set in the first place
        mv = 0;
        for c in concentrations:
            if c=='NA':
                mv += 1;
        if mv==0:
            # Call to R
            try:
                # clear the R workspace
                self.clear_workspace();

                # convert lists to R matrix
                self.make_matrixFromList(concentrations,len(cn_sorted),len(sns_sorted),'concentrations_m');
                #concentrations_r = '';
                #for c in concentrations:
                #    concentrations_r = (concentrations_r + ',' + str(c));
                #concentrations_r = concentrations_r[1:];
                #r_statement = ('concentrations = c(%s)' % concentrations_r);
                #ans = robjects.r(r_statement);
                #r_statement = ('concentrations_m = matrix(concentrations, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cn_sorted),len(sns_sorted)));
                #ans = robjects.r(r_statement);

                # convert lists to R list
                self.make_vectorFromList(sna,'sna');
                #sna_r = '';
                #for c in sna:
                #    sna_r = (sna_r + ',' + '"' + c + '"');
                #sna_r = sna_r[1:];
                #r_statement = ('sna = c(%s)' % sna_r);
                #ans = robjects.r(r_statement);
                self.make_vectorFromList(replicates,'replicates');
                #replicates_r = '';
                #for c in replicates:
                #    replicates_r = (replicates_r + ',' + str(c));
                #replicates_r = replicates_r[1:];
                #r_statement = ('replicates = c(%s)' % replicates_r);
                #ans = robjects.r(r_statement);

                # make the R factor list
                self.make_factorList('sna','replicates','concentrations_l');
                #r_statement = ('concentrations_l = list(sna=sna,replicates=replicates)');
                #ans = robjects.r(r_statement);

                #convert to Expression Set
                self.convert_matrix2ExpressionSet(matrix_I='concentrations_m',vlist_I='concentrations_l',es_O='eS');
                #r_statement = ('eS = neweS(concentrations_m,concentrations_l)');
                #ans = robjects.r(r_statement);

                # estimate the parameters
                self.call_tranest('eS','tranpar',
                                mult=mult,
                                lowessnorm=lowessnorm
                                );
                #r_statement = ('tranpar = tranest(eS, mult=TRUE)'); # Matches metabo-analyst and produces the most uniform distribution
                #ans = robjects.r(r_statement);

                # transform the expression set
                self.call_transeS('eS','tranpar','eS_transformed');
                #r_statement = ('eS_transformed <- transeS(eS, tranpar$lambda, tranpar$alpha)');
                #ans = robjects.r(r_statement);

                # extract out data matrices
                concentrations_glog = self.extract_expressionSet('eS_transformed');
                #r_statement = ('exprs(eS_transformed)');
                #ans = robjects.r(r_statement);
                #concentrations_glog = numpy.array(ans);

                # convert array back to dict
                data_O = [];
                cnt = 0;
                for c in range(len(cn_sorted)):
                    for s in range(len(sns_sorted)):
                        if isinstance(concentrations_glog[c,s], (int, float, complex)):
                            data_tmp = {};
                            data_tmp['sample_name_short'] = sns_sorted[s]
                            data_tmp['component_group_name'] = cgn[c]
                            data_tmp['component_name'] = cn_sorted[c]
                            data_tmp['calculated_concentration_units'] = calculated_concentration_units[c] + '_glog_normalized';
                            data_tmp['calculated_concentration'] = concentrations_glog[c,s];
                            data_tmp['experiment_id'] = experiment_ids[s];
                            data_tmp['time_point'] = time_points[s];
                            data_tmp['analysis_id'] = analysis_ids[s];
                            data_tmp['imputation_method'] = None;
                            data_tmp['normalization_method'] = 'glog';
                            data_tmp['normalization_ooptions'] = {'mult':"TRUE",'lowessnorm':"FALSE"};
                            data_tmp['comment_'] = None;
                            data_tmp['used_'] = True;
                            data_O.append(data_tmp);
                            cnt+=1;
                        else:
                            print('concentration value is not a number.');
                #for c in range(len(sns_sorted)):
                #    for r in range(len(cgn_sorted)):
                        #if isinstance(concentrations_glog[r,c], (int, long, float, complex)):
                        #    data_tmp = {};
                        #    data_tmp['sample_name_short'] = sns_sorted[c]
                        #    data_tmp['component_name'] = cgn_sorted[r]
                        #    data_tmp['calculated_concentration'] = concentrations_glog[r,c];
                        #    #sns_O.append(sns_sorted[c]);
                        #    #cn_O.append(cgn_sorted[r]);
                        #    #cc_O.append(ans[c*len(cgn_sorted)+r]);
            except Exception as e:
                print(e);
                exit(-1);

            # reshape original concentrations
            concentrations_original = numpy.array(concentrations);
            concentrations = concentrations_original.reshape(len(cn_sorted),len(sns_sorted));
            return data_O, concentrations, concentrations_glog;
        else:
            print('missing values found in data!');

    def convert_matrix2ExpressionSet(self,matrix_I,vlist_I,es_O):
        '''
        Convert a matrix to an expressions set in R
        INPUT:
        matrix_I = string, matrix variable in the R workspace
        vlist_I = string, list variable in the R workspace
        OUTPUT:
        es_O = string, name of the expressionSet variable in the R workspace

        Description
        This function converts a data matrix into an ExpressionSet object.
        Usage
        neweS(mat, vlist, vlabel = as.list(names(vlist)))
        Arguments
        mat A data matrix to be converted.
        vlist A list, each component of which describes a factor in the experimental design.
        vlabel A list of labels for each component of vlist.
        Details
        Each element of a component of vlist corresponds to a column of mat. See vlist for an example.
        Value
        eset An ExpressionSet object.
        '''
        
        try:
            r_statement = ('%s = neweS(%s,%s)' %(es_O,matrix_I,vlist_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def make_factorList(self,sna,replicates,list_O):
        '''
        make factor list for LMGene
        INPUT:
        sna = string, name of the R workspace variable
        replicates = string, name of the R workspace variable
        OUTPUT:
        list_O = string, name of the R workspace variable                
        '''
        try:
            r_statement = ('%s = list(sna=%s,replicates=%s)' %(list_O,sna,replicates));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_expressionSet(self,es_I):
        '''
        Extract out data matrices from an expression set
        INPUT:
        es_I = string, name of the expression set in the R workspace
        OUTPUT:
        data_O = numpy matrix
        '''
        data_O = None;
        try:
            r_statement = ('exprs(eS_transformed)');
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
        except Exception as e:
            print(e);
            exit(-1);
        return data_O;

    def call_tranest(self,es_I,transpar_O,
            mult="TRUE",
            lowessnorm="FALSE"
            ):
        '''
        estimate the glog transformation parameters
        INPUT:
        es_I = string, name of the expression set R workspace variable
        OUTPUT:
        transpar_O = string, name of the R woskpace variable

        NOTES: r_statement = ('tranpar = tranest(eS, mult=TRUE)'); # Matches metabo-analyst and produces the most uniform distribution

        Description
        Estimates parameters for the glog transformation, by maximum likelihood or by minimizing the
        stability score.

        Usage
        tranest(eS, ngenes = -1, starting = FALSE, lambda = 1000, alpha = 0,
        gradtol = 1e-3, lowessnorm = FALSE, method=1, mult=FALSE, model=NULL,
        SD = FALSE, rank = TRUE, model.based = TRUE, rep.arrays = NULL)

        Arguments
        eS An ExpressionSet object
        ngenes Number of genes to be used in parameter estimation. Default is to use all genes
        unless there are more than 100,000, in which case a subset of 50,000 genes is
        selected at random.
        starting If TRUE, user-specified starting values for lambda and alpha are input to the
        optimization routine
        lambda Starting value for parameter lambda. Ignored unless starting = TRUE
        alpha Starting value for parameter alpha. Ignored unless starting = TRUE
        gradtol A positive scalar giving the tolerance at which the scaled gradient is considered
        close enough to zero to terminate the algorithm
        lowessnorm If TRUE, lowess normalization (using lnorm) is used in calculating the likelihood.
        method Determines optimization method. Default is 1, which corresponds to a Newtontype
        method (see nlm and details.)
        mult If TRUE, tranest will use a vector alpha with one (possibly different) entry per
        sample. Default is to use same alpha for every sample. SD and mult may not
        both be TRUE.
        model Specifies model to be used. Default is to use all variables from eS without
        interactions. See details.
        SD If TRUE, transformation parameters are estimated by minimizing the stability
        score rather than by maximum likelihood. See details.
        rank If TRUE, the stability score is calculated by regressing the replicate standard deviations
        on the ranks of the gene/row means (rather than on the means themselves).
        Ignored unless SD = TRUE
        model.based If TRUE, the stability score is calculated using the standard deviations of residuals
        from the linear model in model. Ignored unless SD = TRUE
        rep.arrays List of sets of replicate arrays. Each element of rep.arrays should be a vector
        with entries corresponding to arrays (columns) in exprs(eS) conducted under
        the same experimental conditions, i.e., with identical rows in pData(eS). Ignored
        unless SD = TRUE and model.based = FALSE
        tranest 19

        Details
        If you have data in a matrix and information about experimental design factors, then you can use
        neweS to convert the data into an ExpressionSet object. Please see neweS for more detail.
        The model argument is an optional character string, constructed like the right-hand side of a formula
        for lm. It specifies which of the variables in the ExpressionSet will be used in the model
        and whether interaction terms will be included. If model=NULL, it uses all variables from the
        ExpressionSet without interactions. Be careful of using interaction terms with factors; this often
        leads to overfitting, which will yield an error.
        The default estimation method is maximum likelihood. The likelihood is derived by assuming that
        there exist values for lambda and alpha such that the residuals from the linear model in model, fit
        to glog-transformed data using those values for lambda and alpha, follow a normal distribution.
        See Durbin and Rocke (2003) for details.
        If SD = TRUE, lambda and alpha are estimated by minimizing the stability score rather than by
        maximum likelihood. The stability score is defined as the absolute value of the slope coefficient
        from the regression of the replicate/residual standard deviation on the gene/row means, or on the
        rank of the gene/row means. If model.based = TRUE, the stability score is calculated using the
        standard deviation of residuals from the linear model in model. Otherwise, the stability score is
        calculated using the pooled standard deviation over sets of replicates in rep.arrays. See Wu and
        Rocke (2009) for details.
        Optimization methods in method are as follows:
        1 = Newton-type method, using nlm
        2 = Nelder-Mead, using optim
        3 = BFGS, using optim
        4 = Conjugate gradients, using optim
        5 = Simulated annealing, using optim (may only be used when mult = TRUE)

        Value
        A list with components:
        lambda Estimate of transformation parameter lambda
        alpha Estimate of transformation parameter alpha

                
        '''
        try:
            r_statement = ('%s = tranest(%s, mult=%s, lowessnorm=%s)'
                           %(transpar_O,es_I,mult,lowessnorm));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def call_transeS(self,es_I,transpar_I,es_O):
        '''
        call transeS Function to apply the glog transform to an expression set.
        INPUT:
        es_I = string, name of the expression set R workspace variable
        transpar_I = string, name of the R woskpace variable
        OUTPUT:
        es_I = string, name of the transformed expression set R workspace variable

        Description
        For each element in the array of expression data, this function applies the glog transform y -> glog
        (y-alpha, lambda). If alpha is a vector, it must have one element for each column in exprs(eS).

        Usage
        transeS(eS, lambda, alpha)

        Arguments
        eS An ExpressionSet or AffyBatch object
        lambda The parameter lambda to be used in the glog transform.
        alpha The alpha parameter(s) for the glog transform. May be a single number used for
        all samples, or a vector with one entry per sample.

        Details
        The glog transformation of a variable y is defined as log(y + sqrt(y^2 + lambda)). Using
        lambda = 0 corresponds to the log transformation, up to a scale factor of 2. (Other, equivalent
        expressions exist for the glog transformation. See Durbin et al. (2002) and Huber et al. (2002) for
        futher details.)
        transeS subtracts a (scalar or vector) parameter alpha prior to application of the glog transformation,
        resulting in the expression log(y - alpha + sqrt((y - alpha)^2 + lambda)).
        The parameters lambda and alpha may be estimated using tranest.

        Value
        Returns an ExpressionSet or AffyBatch object with the expression matrix glog-transformed.
        '''
        try:
            r_statement = ('%s = transeS(%s, %s$lambda, %s$lambda)'
                           %(es_O,es_I,transpar_I,transpar_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

