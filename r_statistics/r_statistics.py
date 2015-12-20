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
        except:
            print('error in R')
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
        except Exception as e:
            print(e);
            exit(-1);
        return data_tmp;
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