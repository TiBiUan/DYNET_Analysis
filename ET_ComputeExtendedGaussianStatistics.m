function [DeltaMu,DeltaLambda,DeltaTheta,KLD,DM1,DM2,DM3,DS11,DS12,DS13,DS22,DS23,DS33] = ET_ComputeExtendedGaussianStatistics(mu1,mu2,sigma1,sigma2)

    % Kullback-Leibler divergence
    KLD = 1/2*(trace(inv(sigma2)*sigma1) - size(sigma1,1) + ...
        (mu2-mu1)'*inv(sigma2)*(mu2-mu1) + log(det(sigma2)/det(sigma1)));

    % Difference in mean (distance between both centers)
    DeltaMu = sqrt(sum((mu1 - mu2).^2));
    
    [u1,s1] = svd(sigma1);
    [u2,s2] = svd(sigma2);
    
    s1 = s1(1);
    s2 = s2(1);
    
    u1 = u1(:,1);
    u2 = u2(:,1);
    
    u1 = u1/norm(u1);
    u2 = u2/norm(u2);
    
    DeltaLambda = s1^2 - s2^2;
    
    DeltaTheta = u1'*u2;
    
    DM1 = mu1(1) - mu2(1);
    DM2 = mu1(2) - mu2(2);
    DM3 = mu1(3) - mu2(3);
    
    DS11 = sigma1(1,1) - sigma2(1,1);
    DS12 = sigma1(1,2) - sigma2(1,2);
    DS13 = sigma1(1,3) - sigma2(1,3);
    DS22 = sigma1(2,2) - sigma2(2,2);
    DS23 = sigma1(2,3) - sigma2(2,3);
    DS33 = sigma1(3,3) - sigma2(3,3);
end