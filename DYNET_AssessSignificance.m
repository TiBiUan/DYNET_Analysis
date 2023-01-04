function [Sign_DM,Sign_DS,p_DS] = DYNET_AssessSignificance(R,N)

    Sign_DM = zeros(3,1);

    a = prctile(R.DM1N,[2.5/N,100-2.5/N]);
    
    if R.DM1 > a(2) || R.DM1 < a(1)
        Sign_DM(1) = sign(R.DM1);
    end
    
    a = prctile(R.DM2N,[2.5/N,100-2.5/N]);
    
    if R.DM2 > a(2) || R.DM2 < a(1)
        Sign_DM(2) = sign(R.DM2);
    end

    
    a = prctile(R.DM3N,[2.5/N,100-2.5/N]);
    
    if R.DM3 > a(2) || R.DM3 < a(1)
        Sign_DM(3) = sign(R.DM3);
    end
    
    Sign_DS = zeros(3,3);

    a = prctile(R.DS11N,[2.5/N,100-2.5/N]);
    
    p_DS = zeros(3,3);
    
    if R.DS11 > a(2) || R.DS11 < a(1)
        Sign_DS(1,1) = sign(R.DS11);
        p_DS(1,1) = N*min([sum(R.DS11N > R.DS11),sum(R.DS11N < R.DS11)]);
    end
    
    a = prctile(R.DS22N,[2.5/N,100-2.5/N]);
    
    if R.DS22 > a(2) || R.DS22 < a(1)
        Sign_DS(2,2) = sign(R.DS22);
        p_DS(2,2) = N*min([sum(R.DS22N > R.DS22),sum(R.DS22N < R.DS22)]);
    end
    
    a = prctile(R.DS33N,[2.5/N,100-2.5/N]);
    
    if R.DS33 > a(2) || R.DS33 < a(1)
        Sign_DS(3,3) = sign(R.DS33);
        p_DS(3,3) = N*min([sum(R.DS33N > R.DS33),sum(R.DS33N < R.DS33)]);
    end
    
    a = prctile(R.DS12N,[2.5/N,100-2.5/N]);
    
    if R.DS12 > a(2) || R.DS12 < a(1)
        Sign_DS(1,2) = sign(R.DS12);
        p_DS(1,2) = N*min([sum(R.DS12N > R.DS12),sum(R.DS12N < R.DS12)]);
    end
    
    a = prctile(R.DS13N,[2.5/N,100-2.5/N]);
    
    if R.DS13 > a(2) || R.DS13 < a(1)
        Sign_DS(1,3) = sign(R.DS13);
        p_DS(1,3) = N*min([sum(R.DS13N > R.DS13),sum(R.DS13N < R.DS13)]);
    end
    
    a = prctile(R.DS23N,[2.5/N,100-2.5/N]);
    
    if R.DS23 > a(2) || R.DS23 < a(1)
        Sign_DS(2,3) = sign(R.DS23);
        p_DS(2,3) = N*min([sum(R.DS23N > R.DS23),sum(R.DS23N < R.DS23)]);
    end
end