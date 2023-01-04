function [Likelihood,LLikelihood] = ET_EvaluateGaussian(x,Mu,Sigma)

    % Number of subjects to evaluate for
    n_subjects = size(x,2);
    
    % For each...
    for s = 1:n_subjects
        
        % Likelihood
        Likelihood(s) = 1/(realmin+sqrt(det(2*pi*Sigma)))   *   ...
            exp( -1/2  *  (x(:,s) - Mu)'    *   (Sigma\(x(:,s) - Mu))      );
    end
    
    LLikelihood = log(Likelihood);
end