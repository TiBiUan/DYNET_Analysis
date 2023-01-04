function [Likelihood,LLikelihood] = ET_EvaluateGaussian_1D(x,Mu,Sigma)

    % Number of subjects to evaluate for
    n_subjects = length(x);
    
    % For each...
    for s = 1:n_subjects
        
        % Likelihood
        Likelihood(s) = 1/(realmin+sqrt(2*pi)*Sigma)   *   ...
            exp( -1/(realmin+2*Sigma^2)  *  (x(s) - Mu)^2);
    end
    
    LLikelihood = log(Likelihood);
end