from .r_dependencies import *
from .r_base import r_base 
class r_anova(r_base):
    def calculate_anova(self,data_I):
        '''calculate the 1-way anova using R's built in Stats package
        Note: 1-way anova is equivalent to an independent t-test'''

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
        # check that there is only one component_name:
        if len(cn_sorted)>1:
            print('more than one component detected!')
            return None,None;
        #check if there were any missing values in the data set in the first place
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
                mean = np.array(ans);
                r_statement = ('tapply(concentrations_v,sna_v,var)'); # calculate the variance
                ans = robjects.r(r_statement);
                var = np.array(ans);
                r_statement = ('tapply(concentrations_v,sna_v,length)'); # calculate the # of samples
                ans = robjects.r(r_statement);
                n = np.array(ans);
                #convert to Data Frame
                r_statement = ('dF = data.frame(concentrations_v,sna_v)');
                ans = robjects.r(r_statement);
                r_statement = ('names(dF) = c("concentrations","sna")');
                ans = robjects.r(r_statement);
                r_statement = ('attach(dF)');
                ans = robjects.r(r_statement);
                # call anova
                r_statement = ('aov.out = aov(concentrations ~ sna, data=dF)'); # call anova
                ans = robjects.r(r_statement);
                r_statement = ('summary(aov.out)'); # anova summary
                ans = robjects.r(r_statement);
                f_stat = ans[0].rx2('F value')[0] # f_value
                pvalue = ans[0].rx2('Pr(>F)')[0] # pvalue
                # other attributes available: ['Df', 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)']
                r_statement = ('TukeyHSD(aov.out)'); # TukeyHSD post Hoc
                ans = robjects.r(r_statement);
                postHocTest = [x for x in ans[0].rownames];
                # ans[0].colnames
                # ['diff', 'lwr', 'upr', 'p adj']
                postHocTest_pvalue = [x for x in np.array(ans[0])[:,3]];
                postHocTest_description = 'TukeyHSD'
                # convert array back to dict
                data_anova = [];
                data_pairwise = [];
                # extract out unique sna's in order
                sna_set = [];
                for s in sna:
                    if not(s in sna_set):
                        sna_set.append(s);
                for r in range(len(cn_sorted)):
                    data_tmp = {};
                    data_tmp['sample_name_abbreviation'] = sna_set;
                    data_tmp['component_name'] = cn_sorted[r];
                    data_tmp['test_stat'] = f_stat;
                    data_tmp['test_description'] = '1-way ANOVA; F value';
                    data_tmp['pvalue'] = pvalue;
                    data_tmp['pvalue_corrected'] = None;
                    data_tmp['pvalue_corrected_description'] = None;
                    data_anova.append(data_tmp);
                # extract out unique sna's in order
                for c1 in range(len(sna_set)):
                    for c2 in range(len(sna_set)):
                        if c1 != c2:
                            # extract out post hoc results
                            PostHocTest_tmp = '';
                            PostHocTest_pvalue_tmp = None;
                            foldChange = None;
                            for i,pht in enumerate(postHocTest):
                                if sna_set[c1] in pht and sna_set[c2] in pht:
                                    PostHocTest_tmp = sna_set[c2];
                                    PostHocTest_pvalue_tmp = postHocTest_pvalue[i];
                                    foldChange = mean[c2]/mean[c1];
                            for r in range(len(cn_sorted)):
                                data_tmp = {};
                                data_tmp['sample_name_abbreviation_1'] = sna_set[c1];
                                data_tmp['sample_name_abbreviation_2'] = PostHocTest_tmp;
                                data_tmp['component_name'] = cn_sorted[r];
                                data_tmp['mean'] = None;
                                #data_tmp['mean'] = mean[c1];
                                #data_tmp['var'] = var[c1];
                                #data_tmp['n'] = n[c1];
                                data_tmp['test_stat'] = None;
                                data_tmp['test_description'] = postHocTest_description;
                                data_tmp['pvalue'] = PostHocTest_pvalue_tmp;
                                data_tmp['pvalue_corrected'] = None;
                                data_tmp['pvalue_corrected_description'] = None;
                                data_tmp['ci_lb'] = None;
                                data_tmp['ci_ub'] = None;
                                data_tmp['ci_level'] = None;
                                data_tmp['fold_change'] = foldChange;
                                data_pairwise.append(data_tmp);
            except Exception as e:
                print(e);
                exit(-1);
        return data_anova,data_pairwise;

    def calculate_aov(self,
            function_I = 'concentrations ~ sna',
            dataFrame_I = 'dF',
            aov_O = 'aov.out',
            ):
        '''call R aov
        INPUT:
        function_I = string or R workspace lm variable
        dataFrame_I = string, R workspace dataframe
        aov_O = string, name of the R aov output variable
        '''
        try:
            # call anova
            r_statement = ('%s = aov(%s, data=%s)'%(aov_O,function_I,dataFrame_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extraction_aov(self,
            aov_I = 'aov.out'):
        '''extract aov output
        INPUT:
        aov_I = string, R aov workspace variable
        OUTPUT:
        f_stat = float
        pvalue = float
        '''
        try:
            # anova summary
            r_statement = ('summary(%s)'%(aov_I));
            ans = robjects.r(r_statement);
            # other attributes available: ['Df', 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)']
            f_stat = ans[0].rx2('F value')[0] # f_value
            pvalue = ans[0].rx2('Pr(>F)')[0] # pvalue
            return f_stat,pvalue;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_tukeyHSD(self,
        aov_I='aov.out',
        ci_level_I = 0.95,
        tukeyhsd_O='tukeyhsd'):
        '''call R TukeyHSD
        INPUT:
        aov_I = R workspace aov output
        tukeyhsd_O = string, name of the R TukeyHSD output variable
        OUTPUT:
        '''
        try:
            # call anova
            r_statement = ('%s = TukeyHSD(%s, conf.level = %s)'%(tukeyhsd_O,aov_I,ci_level_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_tukeyHSD(self,
        tukeyhsd_I='tukeyhsd'):
        '''call R TukeyHSD
        INPUT:
        tukeyhsd_I = string, name of the R TukeyHSD output variable
        OUTPUT:
        diff
        lb
        ub
        pvalue
        '''
        try:
            # call anova
            r_statement = ('%s'%(tukeyhsd_I));
            ans = robjects.r(r_statement);

            labels = [x for x in ans[0].rownames];
            # ans[0].colnames
            # ['diff', 'lwr', 'upr', 'p adj']
            
            postHocTest_diff = np.array(ans[0])[:,0];
            postHocTest_lb = np.array(ans[0])[:,1];
            postHocTest_ub = np.array(ans[0])[:,2];
            postHocTest_pvalue = np.array(ans[0])[:,3];
            #postHocTest_description = 'TukeyHSD'
            return postHocTest_diff,postHocTest_lb,postHocTest_ub,postHocTest_pvalue,labels;
        except Exception as e:
            print(e);
            exit(-1);
