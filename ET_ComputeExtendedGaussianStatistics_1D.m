function [KLD,DM,DS] = ET_ComputeExtendedGaussianStatistics_1D(mu1,mu2,sigma1,sigma2)

    % Kullback-Leibler divergence
    KLD = 1/2*(trace(inv(sigma2)*sigma1) - size(sigma1,1) + ...
        (mu2-mu1)'*inv(sigma2)*(mu2-mu1) + log(det(sigma2)/det(sigma1)));

    
    
    DM = mu1 - mu2;
    
    
    DS = sigma1 - sigma2;
   
end