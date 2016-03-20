from .r_dependencies import *
from .r_base import r_base
class r_dataSplit(r_base):
    '''cross validate using a hold-out of the original data
    INPUT:
    Description
    A series of test/training partitions are created using createDataPartition while createResample
    creates one or more bootstrap samples. createFolds splits the data into k groups while createTimeSlices
    creates cross-validation sample information to be used with time series data.
        
    Usage
    createDataPartition(y,
    times = 1,
    p = 0.5,
    list = TRUE,
    groups = min(5, length(y)))
    createResample(y, times = 10, list = TRUE)
    createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
    createMultiFolds(y, k = 10, times = 5)
    createTimeSlices(y, initialWindow, horizon = 1,
    fixedWindow = TRUE, skip = 0)

    Arguments
    y a vector of outcomes. For createTimeSlices, these should be in chronological
    order.
    times the number of partitions to create
    p the percentage of data that goes to training
    list logical - should the results be in a list (TRUE) or a matrix with the number of
    rows equal to floor(p * length(y)) and times columns.
    30 createDataPartition
    groups for numeric y, the number of breaks in the quantiles (see below)
    k an integer for the number of folds.
    returnTrain a logical. When true, the values returned are the sample positions corresponding
    to the data used during training. This argument only works in conjunction with
    list = TRUE
    initialWindow The initial number of consecutive values in each training set sample
    horizon The number of consecutive values in test set sample
    fixedWindow A logical: if FALSE, the training set always start at the first sample.
    skip An integer specifying how many (if any) resamples to skip to thin the total
    amount.

    Value
    A list or matrix of row position integers corresponding to the training data
    '''
    def createDataPartition(self,
        y, data_R_trainIndex,data_R_train,data_R_test,
        times = 1,
        p = 0.5,
        list = "TRUE",
        groups = "NONE"):
        '''A series of test/training partitions are created using createDataPartition
        INPUT:
        OUTPUT
        '''
        try:
            r_statement = ('%s <- createDataPartition(%s,times=%s,p=%s,list=%s,groups=%s)' %(data_R_trainIndex,y,times,p,list,groups));
            ans = robjects.r(r_statement);
            # make the training and test splits
            r_statement = ('%s <- %s[%s])' %(data_R_train,y,data_R_trainIndex));
            ans = robjects.r(r_statement);
            r_statement = ('%s <- %s[-%s])' %(data_R_test,y,data_R_trainIndex));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def createResample(self,y, data_R_O,times = 10, list = "TRUE"):
        '''creates one or more bootstrap samples
        INPUT:
        OUTPUT
        '''
        try:
            r_statement = ('%s <- createResample(%s,times=%s,list=%s)' %(data_R_O,y,times,list));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
    def createFolds(self,y, data_R_O,k = 10, list = "TRUE", returnTrain = "FALSE"):
        '''splits the data into k groups
        INPUT:
        OUTPUT
        '''
        try:
            r_statement = ('%s <- createFolds(%s,k=%s,list=%s,returnTrain=%s)' %(data_R_O,y,k,list,returnTrain));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
    def createMultiFolds(self,y, data_R_O,k = 10, times = 5):
        '''splits the data into k groups
        INPUT:
        OUTPUT
        '''
        try:
            r_statement = ('%s <- createMultiFolds(%s,k=%s,times=%s)' %(data_R_O,y,k,times));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
    def createTimeSlices(self,y, data_R_O,initialWindow, horizon = 1,
                fixedWindow = "TRUE", skip = 0):
        '''splits the data into k groups
        INPUT:
        OUTPUT
        '''
        try:
            r_statement = ('%s <- createTimeSlices(%s,initialWindow=%s,horizon=%s,fixedWindow=%s,skip=%s)' %(data_R_O,y,initialWindow,horizon,fixedWindow,skip));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_confusionMatrix(self,
            data, reference, data_R_O,
            positive = "NULL",
            dnn = 'c("Prediction", "Reference")',
            prevalence = "NULL"):
        '''
        Description
        Calculates a cross-tabulation of observed and predicted classes with associated statistics.

        Usage
        confusionMatrix(data, ...)
        ## Default S3 method:
        confusionMatrix(data, reference, positive = NULL,
        dnn = c("Prediction", "Reference"),
        prevalence = NULL, ...)
        ## S3 method for class 'table'
        confusionMatrix(data, positive = NULL, prevalence = NULL, ...)        Arguments
        data a factor of predicted classes (for the default method) or an object of class table.
        reference a factor of classes to be used as the true results
        positive an optional character string for the factor level that corresponds to a "positive"
        result (if that makes sense for your data). If there are only two factor levels, the
        first level will be used as the "positive" result.
        dnn a character vector of dimnames for the table
        prevalence a numeric value or matrix for the rate of the "positive" class of the data. When
        data has two levels, prevalence should be a single numeric value. Otherwise,
        it should be a vector of numeric values with elements for each class. The vector
        should have names corresponding to the classes.
        ... options to be passed to table. NOTE: do not include dnn here        Value
        a list with elements
        table the results of table on data and reference
        positive the positive result level
        overall a numeric vector with overall accuracy and Kappa statistic values
        byClass the sensitivity, specificity, positive predictive value, negative predictive value,
        prevalence, detection rate, detection prevalence and balanced accuracy for each
        class. For two class systems, this is calculated once using the positive argument
        '''
        try:
            r_statement = ('%s <- confusionMatrix(%s,%s,positive=%s,dnn=%s,prevalence=%s)' %(data_R_O,data,reference,positive,dnn,prevalence));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);