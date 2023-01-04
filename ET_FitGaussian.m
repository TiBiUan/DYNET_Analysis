%% X is the data at hand (n_samples x n_dims)
% idx contains the indexing information we will use for the GMM case (can
% be several)
function [Mu,Sigma,L_HC,LL_HC,PD,Lambda,Weights] = ET_FitGaussian(X)
    
    % Parameters
    Mu = mean(X)';
    Sigma = cov(X);
    
    % Extraction of the principal directions of the data and weights along
    % them
    [PD,Lambda] = eig(Sigma);
    
    % Sorts to have the first direction first (max variance expressed)
    [~,idx] = sort(diag(Lambda),'descend');
    Lambda = diag(Lambda);
    Lambda = Lambda(idx);
    PD = PD(:,idx);
    
    % Gives the weights (3 x 29) along the three principal directions,
    % having remembered to subtract the mean to the data points prior to
    % projection
    Weights = PD'*(X'-repmat(Mu,1,size(X,1)));

    % For control subjects, likelihood is computed for each region
    [L_HC,LL_HC] = ET_EvaluateGaussian(X',Mu,Sigma);

end