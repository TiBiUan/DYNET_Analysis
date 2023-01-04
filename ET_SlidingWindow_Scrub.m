%% This function computes the sliding window functional connectivity values
% for a given subject data, with specified parameters
%
% Inputs:
% - TC has size n_ROI x n_TP
function [dFC,dFC_vec,Labels_scrubbing] = ET_SlidingWindow_Scrub(TC,W,Delta,Shape,Measure,Tags,n_min)

    n_ROI = size(TC,1);
    n_TP = size(TC,2);

    % First time span (for first window)
    t_start = 1;
    t_end = W;

    % We compute connectivity as long as we can still move the temporal
    % subwindow forward...
    ind_window = 1;
    while t_end <= n_TP

        % Data from the presently considered temporal window, as well as
        % scrubbing-related tags
        tmp_data = TC(:,t_start:t_end);
        tmp_tags = Tags(t_start:t_end);

        % Computation of connectivity for all region pairs according to the 
        % chosen method
        switch Shape

            % Rectangular window (uniform weights of 1)
            case 'Rectangular'
                
                % If we have enough samples, we compute
                if sum(tmp_tags==0) >= n_min
                    dFC(:,:,ind_window) = corr(tmp_data(:,tmp_tags==0)');
                    Labels_scrubbing(ind_window) = 0;
                
                % Else, we return NaN
                else
                    dFC(:,:,ind_window) = NaN(n_ROI,n_ROI);
                    Labels_scrubbing(ind_window) = 1;
                    % disp('!!!');
                end
        end

        % Shift of the temporal window by a step size
        t_end = t_end + Delta;
        t_start = t_start + Delta;
        ind_window = ind_window + 1;
    end

    % Vectorization
    for t = 1:size(dFC,3)
        tmp = squeeze(dFC(:,:,t));
        dFC_vec(:,t) = jUpperTriMatToVec(tmp);
    end
end