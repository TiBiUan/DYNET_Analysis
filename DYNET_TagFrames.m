function [Tags] = DYNET_TagFrames(FD,T)


    Tags = zeros(size(FD));
    
    for s = 1:size(FD,2)
        
        idx_scrub = [];
        
        for t = 1:size(FD,1)
            
            if FD(t,s) > T
                idx_scrub = unique([idx_scrub,t,t-1,t+1,t+2]);
            end
            
        end
        
        idx_scrub(idx_scrub < 1) = [];
        idx_scrub(idx_scrub > size(FD,1)) = [];
        
        Tags(idx_scrub,s) = 1;
    end
end