from .r_dependencies import *
from .r_base import r_base
class r_dataNormalization(r_base):
    def calculate_glogNormalization(self,data_I):
        '''normalize the data using a glog transformation using LMGene'''

        #make the ExpressionSet
        
        #format into R matrix and list objects
        # convert data dict to matrix filling in missing values
        # with 'NA'
        sns = []
        cgn = []
        #replicates = [];
        sample_name_abbreviations = [];
        for d in data_I:
                sns.append(d['sample_name_short']);    
                #replicates.append(d['sample_replicate']);  
                sample_name_abbreviations.append(d['sample_name_abbreviation']) 
                cgn.append(d['component_name']);
        sns_sorted = sorted(set(sns))
        #replicates_sorted = sorted(set(replicates))
        cgn_sorted = sorted(set(cgn))
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

        concentrations = ['NA' for r in range(len(sns_sorted)*len(cgn_sorted))];
        experiment_ids = ['' for r in range(len(sns_sorted)*len(cgn_sorted))];
        time_points = ['' for r in range(len(sns_sorted)*len(cgn_sorted))];
        cnt = 0;
        cnt_bool = True;
        cnt_reps = 0;
        sna = []
        replicates = []
        for c in cgn_sorted:
                for s in sns_sorted:
                    for d in data_I:
                        if d['sample_name_short'] == s and d['component_name'] == c:
                            if d['calculated_concentration']:
                                concentrations[cnt] = d['calculated_concentration'];
                                experiment_ids[cnt] = d['experiment_id'];
                                time_points[cnt] = d['time_point'];
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
                r_statement = ('concentrations_m = matrix(concentrations, nrow = %s, ncol = %s, byrow = TRUE)' %(len(cgn_sorted),len(sns_sorted)));
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
                for c in range(len(cgn_sorted)):
                    for s in range(len(sns_sorted)):
                        if isinstance(concentrations_glog[c,s], (int, float, complex)):
                            data_tmp = {};
                            data_tmp['sample_name_short'] = sns_sorted[s]
                            data_tmp['component_name'] = cgn_sorted[c]
                            data_tmp['calculated_concentration'] = concentrations_glog[c,s];
                            data_tmp['experiment_id'] = experiment_ids[cnt];
                            data_tmp['time_point'] = time_points[cnt];
                            data_O.append(data_tmp);
                            cnt+=1;
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
            concentrations = concentrations_original.reshape(len(cgn_sorted),len(sns_sorted));
            return data_O, concentrations, concentrations_glog;
        else:
            print('missing values found in data!');