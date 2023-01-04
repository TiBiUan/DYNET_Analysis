%% This function determines which is the most likely data case at hand
function [Result] = DYNET_DetermineSubcases(X,X_ref,n_KLD,alpha)

    n_subjects = size(X,1);

    % Number of parameters in both cases (Gaussian and GMM), used then for
    % the BIC computations
    D = size(X,2);
    N = size(X,1);
    N_ref = size(X_ref,1);
    
    P_Gauss = 2*D + D*(D-1)/2;
    
    % Fits a Gaussian for the reference group and for the other (ML sense)
    [Mu_ref,Sigma_ref] = ET_FitGaussian(X_ref);
    
    % Here, we want to loop and perform leave-one-out cross-validation
    for s = 1:n_subjects
        
        % Segmenting into test and train data
        X_train = X;
        X_test = X(s,:);
        X_train(s,:) = [];
        
        % The model is fit on the training data
        [Mu_Gauss{s},Sigma_Gauss{s}] = ET_FitGaussian(X_train);
        
        % The test data point is evaluated
        [~,LL_Gauss_test(s)] = ET_EvaluateGaussian(X_test',Mu_Gauss{s},Sigma_Gauss{s});
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
    
    [Mu_Gauss_all,Sigma_Gauss_all] = ET_FitGaussian(X);
    
    % Kullback-Leibler divergence and other properties
    [DeltaMu,DeltaLambda,DeltaTheta,KLD,DM1,DM2,DM3,DS11,DS12,DS13,DS22,DS23,DS33] = ...
        ET_ComputeExtendedGaussianStatistics(Mu_ref,Mu_Gauss_all,Sigma_ref,Sigma_Gauss_all);

    % We run many null settings to generate a null KLD distribution,
    % and compare it to the actual value. If the actual distance is
    % larger, then we know our second group differs from the reference
    for n = 1:n_KLD

        tmp_data = [X;X_ref];
        id = randperm(N+N_ref);
        tmp_data = tmp_data(id,:);
        tmp_X = tmp_data(1:N,:);
        tmp_X_ref = tmp_data(N+1:end,:);

        [tmp_mu,tmp_sigma] = ET_FitGaussian(tmp_X);
        [tmp_mu_ref,tmp_sigma_ref] = ET_FitGaussian(tmp_X_ref);

        % Is there a significant group difference or not?
        [DeltaMu_null(n),DeltaLambda_null(n),DeltaTheta_null(n),KLD_null(n),...
            DM1_null(n),DM2_null(n),DM3_null(n),DS11_null(n),DS12_null(n),...
            DS13_null(n),DS22_null(n),DS23_null(n),DS33_null(n)] = ...
            ET_ComputeExtendedGaussianStatistics(tmp_mu_ref,tmp_mu,tmp_sigma_ref,tmp_sigma);
    end
        
    DeltaMu_thresh = prctile(DeltaMu_null,100-alpha);
    DeltaLambda_thresh = prctile(DeltaLambda_null,100-alpha);
    DeltaTheta_thresh = prctile(DeltaTheta_null,100-alpha);
    KLD_thresh = prctile(KLD_null,100-alpha);
        
    % Saves the distance results
    Result.KLD = KLD;
    Result.KLDN = KLD_null;
    
    Result.DM = DeltaMu;
    Result.DMN = DeltaMu_null;
    
    Result.DL = DeltaLambda;
    Result.DLN = DeltaLambda_null;
    
    Result.DT = DeltaTheta;
    Result.DTN = DeltaTheta_null;
    
    Result.DM1 = DM1;
    Result.DM1N = DM1_null;
    
    Result.DM2 = DM2;
    Result.DM2N = DM2_null;
    
    Result.DM3 = DM3;
    Result.DM3N = DM3_null;
    
    Result.DS11 = DS11;
    Result.DS11N = DS11_null;
    
    Result.DS12 = DS12;
    Result.DS12N = DS12_null;
   
    Result.DS13 = DS13;
    Result.DS13N = DS13_null;
    
    Result.DS22 = DS22;
    Result.DS22N = DS22_null;
    
    Result.DS23 = DS23;
    Result.DS23N = DS23_null;
    
    Result.DS33 = DS33;
    Result.DS33N = DS33_null;
        
    if KLD > KLD_thresh
        disp('Group difference (KLD)!');
        Result.GD = 1;
    else
        disp('No group difference!');
        Result.GD = 0;
    end
end





