%% This function prepares the (d)FC data when the contralateral Vim is used
% as seed region
function [FC_mat,FC_vec,CV_vec,dFC,FD,Labels_scrubbing,artif_FC_vc,idx_TC_toremove,idx_toscrub,W_subj] = ...
    DYNET_PrepareData_VimSeed(Main_Folder,Subjects,W,Delta,Shape,Measure,...
    n_regions,idx_seeds,idx_targets,IDX_atlas,RM_scans,T_FD,P_FD,n_min,nW_min)

    cd(Main_Folder);

    % In case we have time courses to remove because they have no values
    idx_TC_toremove = [];
    
    idx_toscrub = [];
    
    % Running across all subjects
    for s = 1:length(Subjects)

        s

        % Loading and preparing the time courses
        cd(Subjects{s});
        load('RESULTS','TC3_MotReg');
        
        % Computation of framewise displacement from the text file
        % summarizing motion parameters, where rotational movement is
        % converted into rad (pi/180) and multiplied by the assumed brain
        % radius (50 mm)
        tmp_mot = textread('rp_fMRI.txt');
        tmp_mot(1:RM_scans,:) = [];
        tmp_mot(:,4:6) = tmp_mot(:,4:6)*50*pi/180;
        
        FD(:,s) = [0;sum(abs(diff(tmp_mot)),2)];

        % First, we determine the frames not to consider using FD, that is,
        % the frames with FD > T_FD
        Tags(:,s) = DYNET_TagFrames(FD(:,s),T_FD);
        
        % We compute the percentage of such frames, a criterion to exclude
        % some of the subjects
        P_scrub(s) = sum(Tags(:,s))/size(FD(:,s),1)*100;
        
        % Concatenates the desired atlases together to obtain size n_time x
        % n_regions
        TC = [];
        for i = 1:length(IDX_atlas)
            TC = [TC,TC3_MotReg{IDX_atlas(i)}];
        end

        % Fills in the "to-remove TCs" vector if needed, which specifies
        % whether there are "null" time courses that should be excluded
        if(sum(isnan(TC(:))) > 0)
            idx_TC_toremove = unique([idx_TC_toremove,find(isnan(sum(TC)))]);
        end

        % Performs a sliding window analysis to generate the dynamic functional
        % connectivity time courses for all subjects
        [~,dFC{s},Labels_scrubbing{s}] = ET_SlidingWindow_Scrub(TC(:,1:n_regions)',W,Delta,Shape,Measure,Tags(:,s),n_min);
        
        if length(Labels_scrubbing{s})-sum(Labels_scrubbing{s}) < nW_min
            idx_toscrub = [idx_toscrub,s];
        end
        
        
        % Pearson's correlation coefficient for all pairs of regions
        FC_mat(:,:,s) = corr(TC(Tags(:,s)==0,1:n_regions));

        % Makes it a vector
        FC_vec(:,s) = jUpperTriMatToVec(squeeze(FC_mat(:,:,s)));
    end

    
    idx_toscrub = unique([idx_toscrub,find(P_scrub > P_FD)]);
    
    cd(Main_Folder);

    clear tmp

    % Only retains the seed-to-target pairs
    artif_FC = zeros(n_regions,n_regions);
    artif_FC(idx_seeds,[idx_targets]) = 1;
    artif_FC([idx_targets],idx_seeds,:) = 1;
    artif_FC_vc = logical(jUpperTriMatToVec(artif_FC));

    FC_vec = FC_vec(artif_FC_vc,:);
    
    for s = 1:length(dFC)
        dFC{s} = dFC{s}(artif_FC_vc,:);
    end
    
    for s = 1:length(dFC)
        CV_vec(:,s) = nanstd(dFC{s},[],2);
    end

    W_subj = size(dFC{1},3);
end