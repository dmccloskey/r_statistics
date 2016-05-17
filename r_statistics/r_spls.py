from .r_dependencies import *
from .r_base import r_base

class r_spls(r_base):

    def calculate_sgpls(self,
        spls_O='spls.o',x='x',y='y',K=3,eta=0.8,scale_x="TRUE",
        eps=1e-5, denom_eps=1e-20,
        zero_eps=1e-5,maxsteps=100,
        br="TRUE",ftype='iden'):
        '''
        Description
        Fit a SGPLS classification model.
        Usage
        sgpls( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
        Arguments
        x Matrix of predictors.
        y Vector of class indices.
        K Number of hidden components.
        eta Thresholding parameter. eta should be between 0 and 1.
        scale.x Scale predictors by dividing each predictor variable by its sample standard deviation?
        eps An effective zero for change in estimates. Default is 1e-5.
        denom.eps An effective zero for denominators. Default is 1e-20.
        zero.eps An effective zero for success probabilities. Default is 1e-5.
        maxstep Maximum number of Newton-Raphson iterations. Default is 100.
        br Apply Firths bias reduction procedure?
        ftype Type of Firths bias reduction procedure. Alternatives are "iden" (the approximated
        version) or "hat" (the original version). Default is "iden".
        Details
        The SGPLS method is described in detail in Chung and Keles (2010). SGPLS provides PLS-based
        classification with variable selection, by incorporating sparse partial least squares (SPLS) proposed
        in Chun and Keles (2010) into a generalized linear model (GLM) framework. y is assumed to have
        numerical values, 0, 1, ..., G, where G is the number of classes subtracted by one.
        Value
        A sgpls object is returned. print, predict, coef methods use this object
        '''
        try:
            r_statement = ('%s <- sgpls(x, y, K, eta, \
                            scale.x=TRUE, eps=1e-5, \
                            denom.eps=1e-20, zero.eps=1e-5, \
                            maxstep=100, br=TRUE, ftype="%s")'
                % (spls_O,x,y,K,eta,scale_x,
                    eps, denom_eps,
                    zero_eps,maxsteps,
                    br,ftype));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_splsda(self,
        spls_O='spls.o',x='x',y='y',K=3,eta=0.8,
        kappa=0.5,classifier='lda',scale_x="TRUE",
        ):
        '''
        Package spls:
        Description
        Fit a SPLSDA classification model.
        Usage
        splsda( x, y, K, eta, kappa=0.5,
        classifier=c('lda','logistic'), scale.x=TRUE, ... )
        Arguments
        x Matrix of predictors.
        y Vector of class indices.
        K Number of hidden components.
        eta Thresholding parameter. eta should be between 0 and 1.
        kappa Parameter to control the effect of the concavity of the objective function and the
        closeness of original and surrogate direction vectors. kappa is relevant only for
        multicategory classification. kappa should be between 0 and 0.5. Default is 0.5.
        classifier Classifier used in the second step of SPLSDA. Alternatives are "logistic" or
        "lda". Default is "lda".
        scale.x Scale predictors by dividing each predictor variable by its sample standard deviation?
        ... Other parameters to be passed through to spls.
        Value
        A splsda object is returned. print, predict, coef methods use this object.

        Package caret:
        probMethod = "softmax", "Bayes"

        '''
        try:
            #r_statement = ('%s <- splsda(%s, %s, K=%s, eta=%s, kappa=%s, \
            #        classifier="%s", scale.x=%s)'
            #    % (spls_O,x,y,K,eta,
            #        kappa,classifier,scale_x));
            #BUG IN PACKAGE SPLS:  Error in spls::spls(x, y, ...) : unused argument (classifier = "lda")
            #lda will be used by default
            #r_statement = ('%s <- spls::splsda(%s, %s, K=%s, eta=%s, kappa=%s, \
            #        scale.x=%s)'
            #    % (spls_O,x,y,K,eta,
            #        kappa,scale_x));
            r_statement = ('%s <- caret::splsda(%s, %s, K=%s, eta=%s, kappa=%s, \
                    scale.x=%s)'
                % (spls_O,x,y,K,eta,
                    kappa,scale_x));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_splsda_mixOmics(self,
        spls_O='spls.o',x='x',y='y',
        ncomp=5,max_iter=500,
        tol=1e-06,near_zero_var="TRUE"
        ):
        '''Description
        Function to perform sparse Partial Least Squares to classify samples (supervised analysis) and select variables.

        Usage
        splsda(X, Y, ncomp = 2, keepX = rep(ncol(X), ncomp),
               max.iter = 500, tol = 1e-06, near.zero.var = TRUE)
        Arguments
        X
        numeric matrix of predictors. NAs are allowed.
        Y
        a factor or a class vector for the discrete outcome.
        ncomp
        the number of components to include in the model (see Details).
        keepX
        numeric vector of length ncomp, the number of variables to keep in X-loadings. By default all variables are kept in the model.
        max.iter
        integer, the maximum number of iterations.
        tol
        a positive real, the tolerance used in the iterative algorithm.
        near.zero.var
        boolean, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations.
        Details
        splsda function fit sPLS models with 1, ... ,ncomp components to the factor or class vector Y. The appropriate indicator (dummy) matrix is created.

        '''
        try:
            r_statement = ('%s <- mixOmics::splsda(%s, %s, ncomp = %s, \
                max.iter = %s, tol = %s, near.zero.var = %s)'
                % (spls_O,x,y,ncomp,
                    max_iter,tol,near_zero_var));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_scoresAndLoadings_splsda(
        self,spls_I='spls.o'):
        '''Extract out the scores and loadings
        from splsda object
        INPUT:
        spls_I = string, name of the R workspace variable
        '''
        try:
            r_statement = ('%s'
                % (spls_I));
            ans = robjects.r(r_statement);
            # extract out the objects in ans
            #list(ans.names)
            #   ['x', 'y', 'betahat', 'A', 'betamat',
            #   'new2As', 'mu', 'meanx', 'normx',
            #   'normy', 'eta', 'K', 'kappa', 'select',
            #   'fit', 'projection', 'obsLevels', 'probMethod']
            betahat = np.array(ans.rx2('betahat')); #dim1 = features, dim2 = unique responses
            a = np.array(ans.rx2('A')); #dim1 = feature indices
            betamat = [np.array(new2as) for new2as in ans.rx2('betamat')]; #dim1 = K, dim2 = features, dim3 = responses
            new2As = [np.array(new2as) for new2as in ans.rx2('new2As')]
            mu = np.array(ans.rx2('mu')); #dim1 = K
            meanx = np.array(ans.rx2('meanx')); #dim1 = features
            normx = np.array(ans.rx2('normx')); #dim1 = features
            normy = np.array(ans.rx2('normy')); #dim1 = features
            fit = str(ans.rx2('fit')[0])
            obsLevels = np.array(ans.rx2('obsLevels')); #dim1 = unique responses
            projection = np.array(ans.rx2('projection')); #dim1 = reduced reatures, dim2 = K
            probMethod = str(ans.rx2('probMethod')[0])
            return data;
        except Exception as e:
            print(e);
            exit(-1);  

    def extract_scoresAndLoadings_splsda_mixOmics(
        self,spls_I='spls.o'):
        '''Extract out the scores and loadings
        from splsda object
        INPUT:
        spls_I = string, name of the R workspace variable

        Values
        splsda returns an object of class "splsda", a list that contains the following components:

        X
        the centered and standardized original predictor matrix.
        Y
        the centered and standardized indicator response vector or matrix.
        ind.mat
        the indicator matrix.
        ncomp
        the number of components included in the model.
        keepX
        number of X variables kept in the model on each component.
        mat.c
        matrix of coefficients to be used internally by predict.
        variates
        list containing the variates.
        loadings
        list containing the estimated loadings for the X and     Y variates.
        names
        list containing the names to be used for individuals and variables.
        nzv
        list containing the zero- or near-zero predictors information.
        tol
        the tolerance used in the iterative algorithm, used for subsequent S3 methods
        max.iter
        the maximum number of iterations, used for subsequent S3 methods
        iter
        Number of iterations of the algorthm for each component
        '''
        try:
            r_statement = ('%s'
                % (spls_I));
            ans = robjects.r(r_statement);
            # extract out the objects in ans
            mat_c = np.array(ans.rx2('mat.c')); #dim1 = features, dim2 = unique responses
            variates = [np.array(variate) for variate in ans.rx2('variates')]; 
                #dim1 = , dim2 = responses, dim3 = components
            loadings = [np.array(loading) for loading in ans.rx2('loadings')]; 
                #dim1 = , dim2 = features, dim3 = components
            unique_factors = np.array(ans.rx2('names')[1]); 
            factors_index = np.array(ans.rx2('names')[2]); #dim1 = features
            return mat_c,variates,loadings,names;
        except Exception as e:
            print(e);
            exit(-1);  

    def cv_sgpls(self,):
        '''
        Description
        Draw heatmap of v-fold cross-validated misclassification rates and return optimal eta (thresholding
        parameter) and K (number of hidden components).
        Usage
        cv.sgpls( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
        br=TRUE, ftype='iden', n.core=8 )
        Arguments
        x Matrix of predictors.
        y Vector of class indices.
        fold Number of cross-validation folds. Default is 10-folds.
        K Number of hidden components.
        eta Thresholding parameter. eta should be between 0 and 1.
        scale.x Scale predictors by dividing each predictor variable by its sample standard deviation?
        plot.it Draw the heatmap of cross-validated misclassification rates?
        br Apply Firth’s bias reduction procedure?
        ftype Type of Firth’s bias reduction procedure. Alternatives are "iden" (the approximated
        version) or "hat" (the original version). Default is "iden".
        n.core Number of CPUs to be used when parallel computing is utilized.
        Details
        Parallel computing can be utilized for faster computation. Users can change the number of CPUs to
        be used by changing the argument n.core.
        Value
        Invisibly returns a list with components:
        err.mat Matrix of cross-validated misclassification rates. Rows correspond to eta and
        columns correspond to number of components (K).
        eta.opt Optimal eta.
        K.opt Optimal K.
        '''
        try:
            r_statement = ('%s <- sgpls(%s)'
                % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def cv_splsda(self,
        cvspls_O='cvspls.o',x='x',y='y',
        fold = 10, K="K",eta="eta",
        kappa=0.5,classifier='lda',scale_x="TRUE",
        plot_it="FALSE",n_core=1):
        '''
        Description
        Draw heatmap of v-fold cross-validated misclassification rates and return optimal eta (thresholding
        parameter) and K (number of hidden components).
        Usage
        cv.splsda( x, y, fold=10, K, eta, kappa=0.5,
        classifier=c('lda','logistic'), scale.x=TRUE, plot.it=TRUE, n.core=8 )
        Arguments
        x Matrix of predictors.
        y Vector of class indices.
        fold Number of cross-validation folds. Default is 10-folds.
        K Number of hidden components.
        eta Thresholding parameter. eta should be between 0 and 1.
        kappa Parameter to control the effect of the concavity of the objective function and the
        closeness of original and surrogate direction vectors. kappa is relevant only for
        multicategory classification. kappa should be between 0 and 0.5. Default is 0.5.
        classifier Classifier used in the second step of SPLSDA. Alternatives are "logistic" or
        "lda". Default is "lda".
        scale.x Scale predictors by dividing each predictor variable by its sample standard deviation?
        plot.it Draw the heatmap of the cross-validated misclassification rates?
        n.core Number of CPUs to be used when parallel computing is utilized.
        Details
        Parallel computing can be utilized for faster computation. Users can change the number of CPUs to
        be used by changing the argument n.core.
        Value
        Invisibly returns a list with components:
        err.mat Matrix of cross-validated misclassification rates. Rows correspond to eta and
        columns correspond to number of components (K).
        eta.opt Optimal eta.
        K.opt Optimal K.
        '''
        try:
            r_statement = ('%s <- cv.splsda( %s, %s, fold=%s, K=%s, eta=%s, kappa=%s, \
                        classifier="%s", scale.x=%s, plot.it=%s, n.core=%s )'
                    % (cvspls_O,x,y,fold,K,eta,
                    kappa,classifier,scale_x,
                    plot_it,n_core));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_cv_splsda(self,
            cvspls_I="spls.O"):
        '''
        extract out cv grid
        INPUT:
        OUTPUT:
        '''
        try:
            r_statement = ('%s'
                % (cvspls_I));
            ans = robjects.r(r_statement);
            # extract out the objects in ans
            #ans.names: ['mspemat', 'eta.opt', 'K.opt']
            data = np.array(ans.rx2('mspemat')); #dim1 = eta, dim2 = K
            return data;
        except Exception as e:
            print(e);
            exit(-1);            

    def calculate_spls(self,):
        '''
        '''
        try:
            r_statement = ('%s <- spls(%s)'
                % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def cv_spls(self,):
        '''
        Description
        Draw heatmap of v-fold cross-validated mean squared prediction error and return optimal eta (thresholding
        parameter) and K (number of hidden components).
        Usage
        cv.spls( x, y, fold=10, K, eta, kappa=0.5,
        select="pls2", fit="simpls",
        scale.x=TRUE, scale.y=FALSE, plot.it=TRUE )
        Arguments
        x Matrix of predictors.
        y Vector or matrix of responses.
        fold Number of cross-validation folds. Default is 10-folds.
        K Number of hidden components.
        eta Thresholding parameter. eta should be between 0 and 1.
        kappa Parameter to control the effect of the concavity of the objective function and
        the closeness of original and surrogate direction vectors. kappa is relevant only
        when responses are multivariate. kappa should be between 0 and 0.5. Default is
        0.5.
        select PLS algorithm for variable selection. Alternatives are "pls2" or "simpls".
        Default is "pls2".
        fit PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls",
        "simpls", or "oscorespls". Default is "simpls".
        scale.x Scale predictors by dividing each predictor variable by its sample standard deviation?
        scale.y Scale responses by dividing each response variable by its sample standard deviation?
        plot.it Draw heatmap of cross-validated mean squared prediction error?
        Value
        Invisibly returns a list with components:
        mspemat Matrix of cross-validated mean squared prediction error. Rows correspond to
        eta and columns correspond to the number of components (K).
        eta.opt Optimal eta.
        K.opt Optimal K.
        '''
        try:
            r_statement = ('%s <- spls(%s)'
                % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def calculate_splsCI(self,):
        '''
        Description
        Calculate bootstrapped confidence intervals of coefficients of the selected predictors and generate
        confidence interval plots.
        Usage
        ci.spls( object, coverage=0.95, B=1000,
        plot.it=FALSE, plot.fix="y",
        plot.var=NA, K=object$K, fit=object$fit )
        Arguments
        object A fitted SPLS object.
        coverage Coverage of confidence intervals. coverage should have a number between 0
        and 1. Default is 0.95 (95% confidence interval).
        B Number of bootstrap iterations. Default is 1000.
        plot.it Plot confidence intervals of coefficients?
        plot.fix If plot.fix="y", then plot confidence intervals of the predictors for a given
        response. If plot.fix="x", then plot confidence intervals of a given predictor
        across all the responses. Relevant only when plot.it=TRUE.
        plot.var Index vector of responses (if plot.fix="y") or predictors (if plot.fix="x")
        to be fixed in plot.fix. The indices of predictors are defined among the set of
        the selected predictors. Relevant only when plot.it=TRUE.
        K Number of hidden components. Default is to use the same K as in the original
        SPLS fit.
        fit PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls",
        "simpls", or "osco
        Value
        Invisibly returns a list with components:
        cibeta A list with as many matrix elements as the number of responses. Each matrix
        element is p by 2, where i-th row of the matrix lists the upper and lower bounds
        of the bootstrapped confidence interval of the i-th predictor.
        betahat Matrix of original coefficients of the SPLS fit.
        lbmat Matrix of lower bounds of confidence intervals (for internal use).
        ubmat Matrix of upper bounds of confidence intervals (for internal use).

        '''
        try:
            r_statement = ('%s <- ci.spls(%s)'
                % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
    def correct_splsCI(self,):
        '''
        Description
        Correct initial SPLS coefficient estimates of the selected predictors based on bootstrapped confi-
        dence intervals and draw heatmap of original and corrected coefficient estimates.
        Usage
        correct.spls( object, plot.it=TRUE )
        Arguments
        object An object obtained from the function ci.spls.
        plot.it Draw the heatmap of original coefficient estimates and corrected coefficient estimates?
        Details
        The set of the selected variables is updated by setting the coefficients with zero-containing confi-
        dence intervals to zero.
        Value
        Invisibly returns a matrix of corrected coefficient estimates

        EXAMPLE:
        data(mice)
        # SPLS with eta=0.6 & 1 latent components
        f <- spls( mice$x, mice$y, K=1, eta=0.6 )
        # Calculate confidence intervals of coefficients
        ci.f <- ci.spls(f)
        # Corrected coefficient estimates
        cf <- correct.spls( ci.f )
        cf[20,1:5]
        '''
        try:
            r_statement = ('%s <- correct.spls(%s)'
                % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);