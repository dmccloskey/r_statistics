from .r_dependencies import *
from .r_base import r_base
class r_statistics(r_base):
    def calculate_ave_CV_R(self,data_I):
        # calculate average and CV of data
        # Call to R
        try:
            # convert lists to R objects
            data_R = robjects.FloatVector(data_I);

            data_ave_R = self.stats.ave(data_R);
            data_ave_O = data_ave_R.rx2(1)[0];

            data_var_R = self.stats.var(data_R);
            data_var = data_var_R.rx2(1)[0];
            data_CV_O = sqrt(data_var)/data_ave_O*100;

            return data_ave_O, data_CV_O;
        except Exception as e:
            print(e);
    def calculate_ave_var_R(self,data_I):
        # calculate average and CV of data
        # Call to R
        try:
            # convert lists to R objects
            data_R = robjects.FloatVector(data_I);

            data_ave_R = self.stats.ave(data_R);
            data_ave_O = data_ave_R.rx2(1)[0];

            data_var_R = self.stats.var(data_R);
            data_var_O = data_var_R.rx2(1)[0];

            return data_ave_O, data_var_O;
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_pairwiseTTest(self,data_I,pooled_sd_I = "FALSE", paired_I="TRUE",padjusted_method_I = "bonferroni",alternative_I = "two.sided"):
        '''calculate a pairwise t-test using R's built in Stats package
        padjusted_methods: ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
        alternative_tests: ("greater","less","two.sided")
        Note pooled_sd and paired cannot both be True
        '''

        #make the dataFrame
        
        #format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        sns = []
        cn = []
        for d in data_I:
                sns.append(d['sample_name_short']);    
                cn.append(d['component_name']);
        sns_sorted = sorted(set(sns))
        cn_sorted = sorted(set(cn))
        concentrations = ['NA' for r in range(len(sns_sorted)*len(cn_sorted))];
        cnt = 0;
        cnt_bool = True;
        sna = []
        for c in cn_sorted:
                for s in sns_sorted:
                    for d in data_I:
                        if d['sample_name_short'] == s and d['component_name'] == c:
                            if d['calculated_concentration']:
                                concentrations[cnt] = d['calculated_concentration'];
                                if cnt_bool:
                                    sna.append(d['sample_name_abbreviation']);
                                break;
                    cnt = cnt+1
                cnt_bool = False;
        if len(cn_sorted)>0:
            print('more than one component detected!')
            return None;
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
                r_statement = ('concentrations_v = c(%s)' % concentrations_r);
                ans = robjects.r(r_statement);
                # convert lists to R list
                sna_r = '';
                for c in sna:
                    sna_r = (sna_r + ',' + '"' + c + '"');
                sna_r = sna_r[1:];
                r_statement = ('sna_v = c(%s)' % sna_r);
                ans = robjects.r(r_statement);
                # get basic stats
                mean = None; #same order as sna
                var = None;
                n = None;
                r_statement = ('tapply(concentrations_v,sna_v,mean)'); # calculate the mean
                ans = robjects.r(r_statement);
                mean = numpy.array(ans);
                r_statement = ('tapply(concentrations_v,sna_v,var)'); # calculate the variance
                ans = robjects.r(r_statement);
                var = numpy.array(ans);
                r_statement = ('tapply(concentrations_v,sna_v,length)'); # calculate the # of samples
                ans = robjects.r(r_statement);
                n = numpy.array(ans);
                #convert to Data Frame
                r_statement = ('dF = data.frame(concentrations_v,sna_v)');
                ans = robjects.r(r_statement);
                r_statement = ('names(dF) = c("concentrations","sna")');
                ans = robjects.r(r_statement);
                r_statement = ('attach(dF)');
                ans = robjects.r(r_statement);
                # call paired T-test with without correction
                r_statement = ('pairwise.t.test(concentrations_v, sna_v, p.adjust.method = "none", pool.sd = %s, paired = %s, alternative = "%s")' %(pooled_sd_I, paired_I ,alternative_I));
                ans = robjects.r(r_statement);
                test_description = ans.rx('method')[0][0]
                pvalues = numpy.array(ans.rx('p.value')[0]);
                rownames = numpy.array(ans[2].rownames);
                colnames = numpy.array(ans[2].colnames);
                # call paired T-test with correction
                r_statement = ('pairwise.t.test(concentrations_v, sna_v, p.adjust.method = "%s", pool.sd = %s, paired = %s, alternative = "%s")' %(padjusted_method_I, pooled_sd_I, paired_I ,alternative_I))
                ans = robjects.r(r_statement);
                test_description = ans.rx('method')[0][0]
                pvalues_adjusted = numpy.array(ans.rx('p.value')[0]);
                pvalue_adjusted_description = ans.rx('p.adjust.method')[0][0]
                rownames_adjusted = numpy.array(ans[2].rownames);
                colnames_adjusted = numpy.array(ans[2].colnames);
                # convert array back to dict
                data_pairwise = [];
                # extract out unique sna's in order
                sna_set = [];
                for s in sna:
                    if not(s in sna_set):
                        sna_set.append(s);
                # extract out unique sna's in order
                for c1 in range(len(rownames)):
                    for c2 in range(len(colnames)):
                        if c1 != c2 and pvalues[c1,c2]!='NA':
                            # extract out post hoc results
                            pair = colnames[c2];
                            pvalue = pvalues[c1,c2];
                            pvalue_adjusted = pvalues_adjusted[c1,c2];
                            #foldChange = mean[c2]/mean[c1];
                            for r in range(len(cn_sorted)):
                                data_tmp = {};
                                data_tmp['sample_name_abbreviation_1'] = rownames[c1];
                                data_tmp['sample_name_abbreviation_2'] = pair;
                                data_tmp['component_name'] = cn_sorted[r];
                                #data_tmp['mean'] = mean[c1];
                                #data_tmp['var'] = var[c1];
                                #data_tmp['n'] = n[c1];
                                data_tmp['test_stat'] = None;
                                data_tmp['test_description'] = test_description;
                                data_tmp['pvalue'] = pvalue;
                                data_tmp['pvalue_corrected'] = pvalue_adjusted;
                                data_tmp['pvalue_corrected_description'] = pvalue_adjusted_description;
                                #data_tmp['fold_change'] = foldChange;
                                data_pairwise.append(data_tmp);
            except Exception as e:
                print(e);
                exit(-1);
        return data_anova,data_pairwise;
    def calculate_twoSampleTTest(self,data_1_I, data_2_I, alternative_I = "two.sided", mu_I = 0, paired_I="TRUE", var_equal_I = "TRUE", ci_level_I = 0.95, padjusted_method_I = "bonferroni"):
        '''calculate a two Sample t-test using R's built in Stats package
        padjusted_methods: ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
        alternative_tests: ("greater","less","two.sided")
        '''
        #make the dataFrame

        #format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        concentrations_1 = [];
        for d in data_1_I:
            if d:
                concentrations_1.append(d);
            else:
                concentrations_1.append('NA')
        concentrations_2 = [];
        for d in data_2_I:
            if d:
                concentrations_2.append(d);
            else:
                concentrations_2.append('NA')
        # Call to R
        try:
            # convert lists to R lists
            concentrations_1_r = '';
            for c in concentrations_1:
                concentrations_1_r = (concentrations_1_r + ',' + str(c));
            concentrations_1_r = concentrations_1_r[1:];
            r_statement = ('concentrations_1_v = c(%s)' % concentrations_1_r);
            ans = robjects.r(r_statement);
            concentrations_2_r = '';
            for c in concentrations_2:
                concentrations_2_r = (concentrations_2_r + ',' + str(c));
            concentrations_2_r = concentrations_2_r[1:];
            r_statement = ('concentrations_2_v = c(%s)' % concentrations_2_r);
            ans = robjects.r(r_statement);
            # call paired T-test without correction
            r_statement = ('t.test(concentrations_1_v,concentrations_2_v, alternative = "%s", mu = %s,  paired = %s, var.equal = %s, conf.level = %s)'\
                %(alternative_I,mu_I, paired_I, var_equal_I ,ci_level_I));
            ans = robjects.r(r_statement);
            test_stat = ans.rx2('statistic')[0]
            test_description = ans.rx2('method')[0]
            pvalue = ans.rx2('p.value')[0]
            mean = ans.rx2('estimate')[0]
            ci = numpy.array(ans.rx2('conf.int'))
            # adjust the p-value
            r_statement = ('p.adjust(%s, method = "%s")' %(pvalue,padjusted_method_I));
            ans = robjects.r(r_statement);
            pvalue_adjusted = ans[0]
            pvalue_adjusted_description = padjusted_method_I
            # extract out data
            data_tmp = {};
            data_tmp['mean'] = mean;
            data_tmp['ci_lb'] = ci[0];
            data_tmp['ci_ub'] = ci[1];
            data_tmp['ci_level'] =ci_level_I;
            data_tmp['test_stat'] = test_stat;
            data_tmp['test_description'] = test_description;
            data_tmp['pvalue'] = pvalue;
            data_tmp['pvalue_corrected'] = pvalue_adjusted;
            data_tmp['pvalue_corrected_description'] = pvalue_adjusted_description;
            return data_tmp;
        except Exception as e:
            print(e);
            #exit(-1);
            return None;
    def calculate_oneSampleTTest(self,data_1_I, alternative_I = "two.sided", mu_I = 0, paired_I="TRUE", var_equal_I = "TRUE", ci_level_I = 0.95, padjusted_method_I = "bonferroni"):
        '''calculate a two Sample t-test using R's built in Stats package
        padjusted_methods: ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
        alternative_tests: ("greater","less","two.sided")
        '''
        #make the dataFrame

        #format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        concentrations_1 = [];
        for d in data_1_I:
            if d:
                concentrations_1.append(d);
            else:
                concentrations_1.append('NA')
        # Call to R
        try:
            # convert lists to R lists
            concentrations_1_r = '';
            for c in concentrations_1:
                concentrations_1_r = (concentrations_1_r + ',' + str(c));
            concentrations_1_r = concentrations_1_r[1:];
            r_statement = ('concentrations_1_v = c(%s)' % concentrations_1_r);
            ans = robjects.r(r_statement);
            # call paired T-test without correction
            r_statement = ('t.test(concentrations_1_v,alternative = "%s", mu = %s,  paired = %s, var.equal = %s, conf.level = %s)'\
                %(alternative_I,mu_I, paired_I, var_equal_I ,ci_level_I));
            ans = robjects.r(r_statement);
            test_stat = ans.rx2('statistic')[0]
            test_description = ans.rx2('method')[0]
            pvalue = ans.rx2('p.value')[0]
            #mean = ans.rx2('estimate')[0]
            ci = numpy.array(ans.rx2('conf.int'))
            # adjust the p-value
            r_statement = ('p.adjust(%s, method = "%s")' %(pvalue,padjusted_method_I));
            ans = robjects.r(r_statement);
            pvalue_adjusted = ans[0]
            pvalue_adjusted_description = padjusted_method_I
            # get basic stats
            mean = None; #same order as sna
            var = None;
            n = None;
            r_statement = ('mean(concentrations_1_v)'); # calculate the mean
            ans = robjects.r(r_statement);
            mean = ans[0];
            r_statement = ('var(concentrations_1_v)'); # calculate the variance
            ans = robjects.r(r_statement);
            var = ans[0];
            r_statement = ('length(concentrations_1_v)'); # calculate the # of samples
            ans = robjects.r(r_statement);
            n = ans[0];

            # convert array back to dict
            data_tmp = {};
            data_tmp['mean'] = mean;
            data_tmp['var'] = var;
            data_tmp['cv'] = sqrt(var)/abs(mean)*100 #glog normalization will have negative values
            data_tmp['n'] = n;
            data_tmp['ci_lb'] = ci[0];
            data_tmp['ci_ub'] = ci[1];
            data_tmp['ci_level'] =ci_level_I;
            data_tmp['test_stat'] = test_stat;
            data_tmp['test_description'] = test_description;
            data_tmp['pvalue'] = pvalue;
            data_tmp['pvalue_corrected'] = pvalue_adjusted;
            data_tmp['pvalue_corrected_description'] = pvalue_adjusted_description;
        except Exception as e:
            print(e);
            return None;
            #exit(-1);
        return data_tmp;
    def calculate_twoSampleWilcoxonRankSumTest(self,data_1_I, data_2_I,
            alternative_I = "two.sided", mu_I = 0, paired_I="TRUE",
            exact_I = "NULL",correct_I = "TRUE",
            ci_int_I = "TRUE", ci_level_I = 0.95, padjusted_method_I = "bonferroni"
            ):
        ''' 
        Calculate the Wilcoxon Rank Sum and Signed Rank Tests using R stats
        https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
        INPUT:
        OUTPUT:
        DESCRIPTION:
        wilcox.test(x, ...)

        ## Default S3 method:
        wilcox.test(x, y = NULL,
                    alternative = c("two.sided", "less", "greater"),
                    mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                    conf.int = FALSE, conf.level = 0.95, ...)

        ## S3 method for class 'formula'
        wilcox.test(formula, data, subset, na.action, ...)
        
        x	
        numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.

        y	
        an optional numeric vector of data values: as with x non-finite values will be omitted.

        alternative	
        a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.

        mu	
        a number specifying an optional parameter used to form the null hypothesis. See Details.

        paired	
        a logical indicating whether you want a paired test.

        exact	
        a logical indicating whether an exact p-value should be computed.

        correct	
        a logical indicating whether to apply continuity correction in the normal approximation for the p-value.

        conf.int	
        a logical indicating whether a confidence interval should be computed.

        conf.level	
        confidence level of the interval.

        formula	
        a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs a factor with two levels giving the corresponding groups.

        data	
        an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).

        subset	
        an optional vector specifying a subset of observations to be used.

        na.action	
        a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
        '''
        try:
            # clear the workspace
            self.clear_workspace();
            #make R vectors
            data_1 = 'data1';
            self.make_vectorFromList(data_1_I,data_1);
            data_2 = 'data2';
            self.make_vectorFromList(data_1_I,data_2);
            #calculate wilcox.text
            data_wt = 'datawt';
            try:
                self.calculate_wilcoxonTest(
                    data_1,data_2,data_wt,
                    alternative_I = alternative_I, mu_I = mu_I, paired_I=paired_I,
                    exact_I = exact_I,correct_I = correct_I,
                    ci_int_I = ci_int_I, ci_level_I = ci_level_I);
            except Exception as e:
                print(e);
                #retry using an alternative method
                try:
                    self.calculate_wilcoxonExact(
                        data_1,data_2,data_wt,
                        alternative_I = alternative_I, mu_I = mu_I, paired_I=paired_I,
                        exact_I = exact_I,correct_I = correct_I,
                        ci_int_I = ci_int_I, ci_level_I = ci_level_I);
                except Exception as e:
                    print(e);
                    return None;
            #extract variables from the R workspace
            data = self.extract_wilcoxonTest(data_wt);
            #adjust the p-value
            if data['p.value'] is None:
                pvalue_adjusted=None;
            else:
                pvalue_O = 'p.value.corrected';
                pvalue_adjusted = self.calculate_pValueCorrected(data['p.value'],pvalue_O,method_I = padjusted_method_I);
            #extract out the values into listDicts
            pvalue_adjusted_description = padjusted_method_I
            # extract out data
            data_tmp = {};
            data_tmp['mean'] = data['estimate'];
            data_tmp['ci_lb'] = data['conf.int'][0];
            data_tmp['ci_ub'] = data['conf.int'][1];
            data_tmp['ci_level'] = ci_level_I;
            data_tmp['test_stat'] = data['statistic'];
            data_tmp['test_description'] = data['method'];
            data_tmp['pvalue'] = data['p.value'];
            data_tmp['pvalue_corrected'] = pvalue_adjusted;
            data_tmp['pvalue_corrected_description'] = pvalue_adjusted_description;
            return data_tmp;
        except Exception as e:
            print(e);
            return None;

    def calculate_wilcoxonTest(self,
            data_1_I, data_2_I,data_O,
            alternative_I = "two.sided", mu_I = 0, paired_I="TRUE",
            exact_I = "NULL",correct_I = "TRUE",
            ci_int_I = "TRUE", ci_level_I = 0.95, padjusted_method_I = "bonferroni",
            ):
        '''
        call wilcox.text from R
        INPUT:
        data_1_I = string, r workspace variable
        data_2_I = string, r workspace variable
        ...
        OUTPUT:
        data_O = string, r workspace variable
        DESCRIPTION:
        wilcox.test(x, ...)

        ## Default S3 method:
        wilcox.test(x, y = NULL,
                    alternative = c("two.sided", "less", "greater"),
                    mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                    conf.int = FALSE, conf.level = 0.95, ...)

        ## S3 method for class 'formula'
        wilcox.test(formula, data, subset, na.action, ...)
        
        ARGUMENTS:
        x	
        numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.

        y	
        an optional numeric vector of data values: as with x non-finite values will be omitted.

        alternative	
        a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.

        mu	
        a number specifying an optional parameter used to form the null hypothesis. See Details.

        paired	
        a logical indicating whether you want a paired test.

        exact	
        a logical indicating whether an exact p-value should be computed.

        correct	
        a logical indicating whether to apply continuity correction in the normal approximation for the p-value.

        conf.int	
        a logical indicating whether a confidence interval should be computed.

        conf.level	
        confidence level of the interval.

        formula	
        a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs a factor with two levels giving the corresponding groups.

        data	
        an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).

        subset	
        an optional vector specifying a subset of observations to be used.

        na.action	
        a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").

        '''
        try:
            r_statement = ('%s = wilcox.test(%s, %s, alternative = "%s",mu = %s, paired = %s, exact = %s, correct = %s, conf.int = %s, conf.level = %s)'
                %(data_O,data_1_I,data_2_I,alternative_I,mu_I,paired_I,exact_I,correct_I,ci_int_I,ci_level_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            raise(e);
            ## retry not calculating the confidence interval
            #try:
            #    r_statement = ('%s = wilcox.test(%s, %s, alternative = "%s",mu = %s, paired = %s, exact = %s, correct = %s, conf.int = %s, conf.level = %s)'
            #        %(data_O,data_1_I,data_2_I,alternative_I,mu_I,paired_I,exact_I,correct_I,"FALSE",ci_level_I));
            #    ans = robjects.r(r_statement);
            #except Exception as ae:
            #    print(ae);
            #    exit(-1);

    def extract_wilcoxonTest(self,
            data_I
            ):
        '''
        extract variables from R

        INPUT:
        data_I = string, r workspace variable
        ...
        OUTPUT:
        data_O = {} of values
        
        returned variable is class "htest" with the following values
        VALUES:
        
        statistic	
        the value of the test statistic with a name describing it.

        parameter	
        the parameter(s) for the exact distribution of the test statistic.

        p.value	
        the p-value for the test.

        null.value	
        the location parameter mu.

        alternative	
        a character string describing the alternative hypothesis.

        method	
        the type of test applied.

        data.name	
        a character string giving the names of the data.

        conf.int	
        a confidence interval for the location parameter. (Only present if argument conf.int = TRUE.)

        estimate	
        an estimate of the location parameter. (Only present if argument conf.int = TRUE.)

        '''
        data_O = None;
        try:
            r_statement = ('%s' %(data_I));
            ans = robjects.r(r_statement);
            #extract out values
            statistic = ans.rx2('statistic')[0];
            if str(type(ans.rx2('parameter')))=="<class 'rpy2.rinterface.RNULLType'>":
                parameter = None;
            else:
                parameter = ans.rx2('parameter')[0];
            if numpy.isnan(ans.rx2('p.value')[0]):
                pValue = None;
            else:
                pValue = ans.rx2('p.value')[0];
            nullValue = ans.rx2('null.value')[0];
            alternative = ans.rx2('alternative')[0];
            method = ans.rx2('method')[0];
            dataName = ans.rx2('data.name')[0];
            if str(type(ans.rx2('conf.int')))=="<class 'rpy2.rinterface.RNULLType'>":
                ci = [None,None];
            else:
                ci = numpy.array(ans.rx2('conf.int'))
            if str(type(ans.rx2('estimate')))=="<class 'rpy2.rinterface.RNULLType'>":
                estimate = None;
            else:
                estimate = ans.rx2('estimate')[0]
            #copy to a dictionary
            data_O = {"statistic":statistic,
                "parameter":parameter,
                "p.value":pValue,
                "null.value":nullValue,
                "alternative":alternative,
                "method":method,
                "data.name":dataName,
                "conf.int":ci,
                "estimate":estimate,
                };

        except Exception as e:
            print(e);
            exit(-1);
        return data_O;

    def calculate_kolmogorovSmirnovTest(self,
            data_1_I, data_2_I,data_O,
            alternative_I = "two.sided",
            exact_I = "NULL"):
        '''
        Perform a one- or two-sample Kolmogorov-Smirnov test.
        https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ks.test.html

        Usage

        ks.test(x, y, ...,
                alternative = c("two.sided", "less", "greater"),
                exact = NULL)

        Can be used to test if the sample comes from a normal distribution
        e.g.
        
        require(graphics)
        x <- rnorm(50)
        y <- runif(30)
        # Do x and y come from the same distribution?
        ks.test(x, y)
        # Does x come from a shifted gamma distribution with shape 3 and rate 2?
        ks.test(x+2, "pgamma", 3, 2) # two-sided, exact
        ks.test(x+2, "pgamma", 3, 2, exact = FALSE)
        ks.test(x+2, "pgamma", 3, 2, alternative = "gr")
        '''
        try:
            # clear the workspace
            self.clear_workspace();
            #make R vectors
            data_1 = 'data1';
            self.make_vectorFromList(data_1_I,data_1);
            data_2 = 'data2';
            self.make_vectorFromList(data_1_I,data_2);
            #calculate wilcox.text
            data_wt = 'datakst';
            self.calculate_ksTest(
                data_1,data_2,data_O,
                alternative_I = alternative_I,
                exact_I = exact_I);
            #extract variables from the R workspace
            data = self.extract_ksTest(data_wt);
            #adjust the p-value
            pvalue_O = 'p.value.corrected';
            pvalue_adjusted = self.calculate_pValueCorrected(data['p.value'],pvalue_O,method_I = padjusted_method_I);
            #extract out the values into listDicts
            pvalue_adjusted_description = padjusted_method_I
            # extract out data
            data_tmp = {};
            data_tmp['mean'] = data['estimate'];
            data_tmp['ci_lb'] = None;
            data_tmp['ci_ub'] = None;
            data_tmp['ci_level'] = ci_level_I;
            data_tmp['test_stat'] = None;
            data_tmp['test_description'] = 'kolmogorov_smirnov_test';
            data_tmp['pvalue'] = data['p.value'];
            data_tmp['pvalue_corrected'] = pvalue_adjusted;
            data_tmp['pvalue_corrected_description'] = pvalue_adjusted_description;
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_ksTest(self):
        '''call ks.test in R
        INPUT:
        OUTPUT:
        
        ks.test(x, y, ...,
                alternative = c("two.sided", "less", "greater"),
                exact = NULL)

        Arguments

        x	
        a numeric vector of data values.

        y	
        either a numeric vector of data values, or a character string naming a cumulative distribution function or an actual cumulative distribution function such as pnorm. Only continuous CDFs are valid.

        ...	
        parameters of the distribution specified (as a character string) by y.

        alternative	
        indicates the alternative hypothesis and must be one of "two.sided" (default), "less", or "greater". You can specify just the initial letter of the value, but the argument name must be give in full. See Details for the meanings of the possible values.

        exact	
        NULL or a logical indicating whether an exact p-value should be computed. See Details for the meaning of NULL. Not available in the two-sample case for a one-sided test or if ties are present.

        Details

        If y is numeric, a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution is performed.

        Alternatively, y can be a character string naming a continuous (cumulative) distribution function, or such a function. In this case, a one-sample test is carried out of the null that the distribution function which generated x is distribution y with parameters specified by ....

        The presence of ties always generates a warning, since continuous distributions do not generate them. If the ties arose from rounding the tests may be approximately valid, but even modest amounts of rounding can have a significant effect on the calculated statistic.

        Missing values are silently omitted from x and (in the two-sample case) y.

        The possible values "two.sided", "less" and "greater" of alternative specify the null hypothesis that the true distribution function of x is equal to, not less than or not greater than the hypothesized distribution function (one-sample case) or the distribution function of y (two-sample case), respectively. This is a comparison of cumulative distribution functions, and the test statistic is the maximum difference in value, with the statistic in the "greater" alternative being D^+ = max[F_x(u) - F_y(u)]. Thus in the two-sample case alternative = "greater" includes distributions for which x is stochastically smaller than y (the CDF of x lies above and hence to the left of that for y), in contrast to t.test or wilcox.test.

        Exact p-values are not available for the two-sample case if one-sided or in the presence of ties. If exact = NULL (the default), an exact p-value is computed if the sample size is less than 100 in the one-sample case and there are no ties, and if the product of the sample sizes is less than 10000 in the two-sample case. Otherwise, asymptotic distributions are used whose approximations may be inaccurate in small samples. In the one-sample two-sided case, exact p-values are obtained as described in Marsaglia, Tsang & Wang (2003) (but not using the optional approximation in the right tail, so this can be slow for small p-values). The formula of Birnbaum & Tingey (1951) is used for the one-sample one-sided case.

        If a single-sample test is used, the parameters specified in ... must be pre-specified and not estimated from the data. There is some more refined distribution theory for the KS test with estimated parameters (see Durbin, 1973), but that is not implemented in ks.test.

        '''
        try:
            r_statement = ('%s' %(data_O));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_ksest(self,
            data_I
            ):
        '''
        extract variables from R

        INPUT:
        data_I = string, r workspace variable
        ...
        OUTPUT:
        data_O = {} of values

        Value

        A list with class "htest" containing the following components:

        statistic	
        the value of the test statistic.

        p.value	
        the p-value of the test.

        alternative	
        a character string describing the alternative hypothesis.

        method	
        a character string indicating what type of test was performed.

        data.name	
        a character string giving the name(s) of the data.

        '''
        data_O = None;
        try:
            r_statement = ('%s' %(data_O));
            ans = robjects.r(r_statement);
            #extract out values
            statistic = ans.rx2('statistic')[0];
            pValue = ans.rx2('p.value')[0];
            alternative = ans.rx2('alternative')[0];
            method = ans.rx2('method')[0];
            dataName = ans.rx2('data.name')[0];
            #copy to a dictionary
            data_O = {"statistic":statistic,
                "p.value":pValue,
                "alternative":alternative,
                "method":method,
                "data.name":dataName,
                };

        except Exception as e:
            print(e);
            exit(-1);
        return data_O;

    def calculate_wilcoxonExact(self,
            data_1_I, data_2_I,data_O,
            alternative_I = "two.sided", mu_I = 0, paired_I="TRUE",
            exact_I = "NULL",correct_I = "TRUE",
            ci_int_I = "TRUE", ci_level_I = 0.95, padjusted_method_I = "bonferroni",
            ):
        '''
        call wilcox.exact from R package exactRankTests
        https://cran.r-project.org/web/packages/exactRankTests/exactRankTests.pdf
        INPUT:
        data_1_I = string, r workspace variable
        data_2_I = string, r workspace variable
        ...
        OUTPUT:
        data_O = string, r workspace variable
        DESCRIPTION:
        wilcox.exact(x, ...)

        ## Default S3 method:
        wilcox.exact(x, y = NULL,
                    alternative = c("two.sided", "less", "greater"),
                    mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                    conf.int = FALSE, conf.level = 0.95, ...)

        ## S3 method for class 'formula'
        wilcox.exact(formula, data, subset, na.action, ...)
        
        ARGUMENTS:
        x	
        numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.

        y	
        an optional numeric vector of data values: as with x non-finite values will be omitted.

        alternative	
        a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.

        mu	
        a number specifying an optional parameter used to form the null hypothesis. See Details.

        paired	
        a logical indicating whether you want a paired test.

        exact	
        a logical indicating whether an exact p-value should be computed.

        correct	
        a logical indicating whether to apply continuity correction in the normal approximation for the p-value.

        conf.int	
        a logical indicating whether a confidence interval should be computed.

        conf.level	
        confidence level of the interval.

        formula	
        a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs a factor with two levels giving the corresponding groups.

        data	
        an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).

        subset	
        an optional vector specifying a subset of observations to be used.

        na.action	
        a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").

        '''
        try:
            r_statement = ('%s = wilcox.exact(%s, %s, alternative = "%s",mu = %s, paired = %s, exact = %s, correct = %s, conf.int = %s, conf.level = %s)'
                %(data_O,data_1_I,data_2_I,alternative_I,mu_I,paired_I,exact_I,correct_I,ci_int_I,ci_level_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            raise(e);