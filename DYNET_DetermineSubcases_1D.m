%% This function determines which is the most likely data case at hand
function [Result] = DYNET_DetermineSubcases_1D(X,X_ref,n_KLD,alpha)

    n_subjects = length(X);

    % Number of parameters in both cases (Gaussian and GMM), used then for
    % the BIC computations
    D = 1;
    N = length(X);
    N_ref = length(X_ref);
    
    P_Gauss = 2*D + D*(D-1)/2;
    
    % Fits a Gaussian for the reference group and for the other (ML sense)
    [Mu_ref,Sigma_ref] = ET_FitGaussian_1D(X_ref);
    
    % Here, we want to loop and perform leave-one-out cross-validation
    for s = 1:n_subjects
        
        s
        
        % Segmenting into test and train data
        X_train = X;
        X_test = X(s);
        X_train(s) = [];
        
        % The model is fit on the training data
        [Mu_Gauss{s},Sigma_Gauss{s}] = ET_FitGaussian_1D(X_train');
        
        % The test data point is evaluated
        [~,LL_Gauss_test(s)] = ET_EvaluateGaussian_1D(X_test',Mu_Gauss{s},Sigma_Gauss{s});
    end

    % Total log-likelihood
    LL_Gauss = sum(LL_Gauss_test);
    
    % Determines the Bayesian Information Criterion for both cases, and
    % also the Akaike Information Criterion
    AIC_Gauss = 2*P_Gauss - 2*LL_Gauss;
    
    BIC_Gauss = P_Gauss*log(N) - 2*LL_Gauss;
    
    Result.Mean_Gauss = Mu_Gauss;
    Result.Sigma_Gauss = Sigma_Gauss;
    
    % We save the evaluation metrics as well
    Result.LL_Gauss = LL_Gauss;
   
    Result.BIC_Gauss = BIC_Gauss;
    
    Result.AIC_Gauss = AIC_Gauss;
    
    [Mu_Gauss_all,Sigma_Gauss_all] = ET_FitGaussian_1D(X);
    
    % Kullback-Leibler divergence (actual)
    [KLD,DM,DS] = ET_ComputeExtendedGaussianStatistics_1D(Mu_ref,...
        Mu_Gauss_all,Sigma_ref,Sigma_Gauss_all);
        
    % We run many null settings to generate a null KLD distribution,
    % and compare it to the actual value. If the actual distance is
    % larger, then we know our second group differs from the reference
    for n = 1:n_KLD

        tmp_data = [X';X_ref'];
        id = randperm(N+N_ref);
        tmp_data = tmp_data(id);
        tmp_X = tmp_data(1:N);
        tmp_X_ref = tmp_data(N+1:end);

        [tmp_mu,tmp_sigma] = ET_FitGaussian_1D(tmp_X);
        [tmp_mu_ref,tmp_sigma_ref] = ET_FitGaussian_1D(tmp_X_ref);

        % Is there a significant group difference or not?
        [KLD_null(n),DM_null(n),DS_null(n)] = ...
            ET_ComputeExtendedGaussianStatistics_1D(tmp_mu_ref,tmp_mu,...
            tmp_sigma_ref,tmp_sigma);
    end
        
    KLD_thresh = prctile(KLD_null,100-alpha);
        
    % Saves the distance results
    Result.KLD = KLD;
    Result.KLDN = KLD_null;
    
    Result.DM = DM;
    Result.DMN = DM_null;
    
    Result.DS = DS;
    Result.DSN = DS_null;
        
    if KLD > KLD_thresh
        disp('Group difference (KLD)!');
        Result.GD = 1;
    else
        disp('No group difference!');
        Result.GD = 0;
    end
end





