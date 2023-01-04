function [Sign_DM,Sign_DS,p_DS] = DYNET_AssessSignificance_SC(R,N)

    Sign_DM = 0;

    a = prctile(R.DMN,[2.5/N,100-2.5/N]);
    
    if R.DM > a(2) || R.DM < a(1)
        Sign_DM = sign(R.DM);
    end
    
    
    
    Sign_DS = 0;

    a = prctile(R.DSN,[2.5/N,100-2.5/N]);
    
    p_DS = 0;
    
    if R.DS > a(2) || R.DS < a(1)
        Sign_DS = sign(R.DS);
        p_DS(1,1) = N*min([sum(R.DSN > R.DS),sum(R.DSN < R.DS)]);
    end
    
    
end