from .r_dependencies import *
from .r_base import r_base
class r_crossValidation(r_base):

    def calculate_trainControl(self,data_I,
            ):
        '''
        INPUT:
        data_I = listDict

        Description
        Control the computational nuances of the train function

        Usage
        trainControl(method = "boot",
        number = ifelse(grepl("cv", method), 10, 25),
        repeats = ifelse(grepl("cv", method), 1, number),
        p = 0.75,
        search = "grid",
        initialWindow = NULL,
        horizon = 1,
        fixedWindow = TRUE,
        verboseIter = FALSE,
        returnData = TRUE,
        returnResamp = "final",
        savePredictions = FALSE,
        classProbs = FALSE,
        summaryFunction = defaultSummary,
        selectionFunction = "best",
        preProcOptions = list(thresh = 0.95, ICAcomp = 3, k = 5),
        sampling = NULL,
        index = NULL,
        indexOut = NULL,
        timingSamps = 0,
        predictionBounds = rep(FALSE, 2),
        seeds = NA,
        adaptive = list(min = 5, alpha = 0.05,
        method = "gls", complete = TRUE),
        trim = FALSE,
        allowParallel = TRUE)

        Arguments
        method The resampling method: "boot", "boot632", "cv", "repeatedcv", "LOOCV",
        "LGOCV" (for repeated training/test splits), "none" (only fits one model to the
        entire training set), "oob" (only for random forest, bagged trees, bagged earth,
        bagged flexible discriminant analysis, or conditional tree forest models), "adaptive_cv",
        "adaptive_boot" or "adaptive_LGOCV"
        number Either the number of folds or number of resampling iterations
        repeats For repeated k-fold cross-validation only: the number of complete sets of folds
        to compute
        verboseIter A logical for printing a training log.
        returnData A logical for saving the data
        returnResamp A character string indicating how much of the resampled summary metrics
        should be saved. Values can be "final", "all" or "none"
        savePredictions
        an indicator of how much of the hold-out predictions for each resample should
        be saved. Values can be either "all", "final", or "none". A logical value
        can also be used that convert to "all" (for true) or "none" (for false). "final"
        saves the predictions for the optimal tuning parameters.
        p For leave-group out cross-validation: the training percentage
        search Either "grid" or "random", describing how the tuning parameter grid is determined.
        See details below.
        initialWindow, horizon, fixedWindow
        possible arguments to createTimeSlices
        classProbs a logical; should class probabilities be computed for classification models (along
        with predicted values) in each resample?
        summaryFunction
        a function to compute performance metrics across resamples. The arguments to
        the function should be the same as those in defaultSummary.
        selectionFunction
        the function used to select the optimal tuning parameter. This can be a name of
        the function or the function itself. See best for details and other options.
        preProcOptions A list of options to pass to preProcess. The type of pre-processing (e.g. center,
        scaling etc) is passed in via the preProc option in train.
        sampling a single character value describing the type of additional sampling that is conducted
        after resampling (usually to resolve class imbalances). Values are "none",
        "down", "up", "smote", or "rose". The latter two values require the DMwR
        and ROSE packages, respectively. This argument can also be a list to facilitate
        custom sampling and these details can be found on the caret package website
        for sampling (link below).
        index a list with elements for each resampling iteration. Each list element is a vector
        of integers corresponding to the rows used for training at that iteration.
        indexOut a list (the same length as index) that dictates which data are held-out for each
        resample (as integers). If NULL, then the unique set of samples not contained in
        index is used.
        timingSamps the number of training set samples that will be used to measure the time for predicting
        samples (zero indicates that the prediction time should not be estimated.
        predictionBounds
        a logical or numeric vector of length 2 (regression only). If logical, the predictions
        can be constrained to be within the limit of the training set outcomes. For
        example, a value of c(TRUE, FALSE) would only constrain the lower end of predictions.
        If numeric, specific bounds can be used. For example, if c(10, NA),
        values below 10 would be predicted as 10 (with no constraint in the upper side).
        seeds an optional set of integers that will be used to set the seed at each resampling
        iteration. This is useful when the models are run in parallel. A value of NA will
        stop the seed from being set within the worker processes while a value of NULL
        will set the seeds using a random set of integers. Alternatively, a list can be used.
        The list should have B+1 elements where B is the number of resamples. The first
        B elements of the list should be vectors of integers of length M where M is the
        number of models being evaluated. The last element of the list only needs to be
        a single integer (for the final model). See the Examples section below and the
        Details section.
        adaptive a list used when method is "adaptive_cv", "adaptive_boot" or "adaptive_LGOCV".
        See Details below.
        trim a logical. If TRUE the final model in object$finalModel may have some components
        of the object removed so reduce the size of the saved object. The
        predict method will still work, but some other features of the model may not
        work. triming will occur only for models where this feature has been implemented.
        allowParallel if a parallel backend is loaded and available, should the function use it?

        Value
        An echo of the parameters specified
        '''
        #handle the input data
        try:
            r_statement = ('%s <- trainControl(%s)' %(data_R_O,y));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);
            
    def calculate_train(self,
            data_R_O,
            x=None, y=None,
            form=None,data=None,
            trControl=None,
            method="rf",tuneGrid="NULL"):
        '''
        INPUT:
        OUTPUT:
        NOTES:
        see http://topepo.github.io/caret/bytag.html for a list of models that can be trained (i.e., method = )

        Description
        This function sets up a grid of tuning parameters for a number of classification and regression
        routines, fits each model and calculates a resampling based performance measure.

        Usage
        train(x, ...)
        ## Default S3 method:
        train(x, y,
        method = "rf",
        preProcess = NULL,
        ...,
        weights = NULL,
        metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
        maximize = ifelse(metric %in% c("RMSE", "logLoss"), FALSE, TRUE),
        trControl = trainControl(),
        tuneGrid = NULL,
        tuneLength = 3)
        ## S3 method for class 'formula'
        train(form, data, ..., weights, subset, na.action, contrasts = NULL)        Arguments
        x an object where samples are in rows and features are in columns. This could
        be a simple matrix, data frame or other type (e.g. sparse matrix). See Details
        below.
        y a numeric or factor vector containing the outcome for each sample.
        form A formula of the form y ~ x1 + x2 + ...
        146 train
        data Data frame from which variables specified in formula are preferentially to be
        taken
        '''
        #handle the input data
        try:
            if not x is None and not y is None:
                r_statement = ('%s <- train(%s,%s,trControl="%s",method=%s,tuneGrid=%s)' %(data_R_O,x,y,trControl,method,tuneGrid));
            elif not form is None and not data is None:
                r_statement = ('%s <- train(%s,data=%s,trControl="%s",method=%s,tuneGrid=%s)' %(data_R_O,form,data,trControl,method,tuneGrid));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_predict(self,object,data_R_O,
        newdata, se_fit = "FALSE", scale = "NULL", df = "Inf",
        interval = 'c("none", "confidence", "prediction")',
        level = 0.95, type = c("response", "terms"),
        terms = "NULL", na_action = "na.pass",
        pred_var = "NULL", weights = 1):
        '''Predict method for Linear Model Fits
        INPUT:
        data_R_I = name of the R workspace variable
        OUTPUT:
        data_R_O = name of the R workspace variable
        data_O = numpy array of data

        Description

        Predicted values based on linear model object.

        Usage

        ## S3 method for class 'lm'
        predict(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                interval = c("none", "confidence", "prediction"),
                level = 0.95, type = c("response", "terms"),
                terms = NULL, na.action = na.pass,
                pred.var = res.var/weights, weights = 1, ...)
        '''
        try:
            r_statement = ('%s <- predict(%s, se.fit = %s ,scale = %s, df = %s, \
                            interval = %s,level = %s, type = %s, terms = %s, \
                            na.action = %s, pred.var = %s, weights = %s)' %(data_R_O,
                                object,se_fit,scale,df,
                                interval,level,type,terms,
                                na_action,pred_var,weights));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);

    def make_expandGrid(self,
            data_str_I,
            data_R_O):
        '''Create a Data Frame from All Combinations of Factor Variables
        INPUT:
        data_str_I = string to pass to expand.grid as input
        OUTPUT:
        data_R_O = name of the R workspace variable
        data_O = numpy array of data

        Description

        Create a data frame from all combinations of the supplied vectors or factors. See the description of the return value for precise details of the way this is done.

        Usage

        expand.grid(..., KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)
        Arguments

        ...	
        vectors, factors or a list containing these.

        KEEP.OUT.ATTRS	
        a logical indicating the "out.attrs" attribute (see below) should be computed and returned.

        stringsAsFactors	
        logical specifying if character vectors are converted to factors.
        '''
        try:
            r_statement = ('%s <- expand.grid(%s)' %(data_R_O,data_str_I));
            ans = robjects.r(r_statement);
            data_O = numpy.array(ans);
            return data_O;
        except Exception as e:
            print(e);
            exit(-1);