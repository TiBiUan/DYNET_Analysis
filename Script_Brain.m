%% This script runs a joint analysis of surface-based morphometry and dynamic
% functional connectivity features on a population of patients with severe
% essential tremor. These patients were scanned before and 1 year after
% Gamma Knife stereotactic radiosurgical thalamotomy of the
% ventro-intermediate nucleus of the thalamus.
%
% The analyses are architectured into three parts:
% 1. Dynamic functional connectivity (dFC) on resting-state functional magnetic
% resonance imaging data: (A) final preprocessing stages, (B) extraction of 
% dFC states on healthy controls matched for age, (C) assignment of ET dFC 
% samples to the HC states, (D) quantification of temporal occurrences and 
% spatial similarity to the states, (E) assessment of pre versus post 
% thalamotomy group difference, (F) assessment of associations between 
% features pre-thalamotomy and clinical recovery
% 2. Surface-based morphometry derived from T1 images through Freesurfer:
% (A) final preprocessing stages, (B) computation of log-likelihood to be
% issued from the HC distribution, (C) coefficient comparison for
% significant regions, (D) assessment of associations between 
% features pre-thalamotomy and clinical recovery
% 3. SBM/dFC correlational analysis
%
% Written by Thomas A.W. Bolton 



%% Part 1. Parameter definitions and data loading

%%%%%%%%%% First, loading of the data relevant for resting-state fMRI
%%%%%%%%%% analyses. In this case, we consider regional time courses of
%%%%%%%%%% activity in the Schaefer 400 (17 networks), Tian S3 and AAL
%%%%%%%%%% atlases

% The paths to 27 BASE and 23 YEAR recordings to consider
load('BSYS_new.mat');

% Same for 18 healthy controls (HCs)
load('./Data/Subjects_HC.mat');

% Note: some of these subjects are discarded, in what follows, because of 
% unclean data. This is why the numbers reported in the manuscript are
% lower (14 HCs, 18 in each ET group)

% The 566 labels of all the regions from the three atlases that I decided
% to use (Schaefer 400 17 networks, Tian S3 and AAL)
load('./Data/Labels_DYNET_finalized.mat');

% Note: this is trimmed to a lower number later (e.g., removal of cortical
% AAL regions because better parcels are in the Schaefer atlas)

% This contains the names of the region pairs involved in the connections
% to the Vim; there are 475 such regions at first
load('./Data/Labels_DYNET1.mat');
load('./Data/Labels_DYNET2.mat');

% This contains the TSTH scores for 41 subjects (all that survived visual
% inspection)
load('./Data/TSTH_DYNET.mat');

% Same for age and the other variables
load('./Data/Age_DYNET.mat');
load('./Data/SymptomsDuration_DYNET.mat');
load('./Data/LesionVolume_DYNET.mat');
load('./Data/TremorStop_DYNET.mat');
load('./Data/TSTHP_DYNET.mat');
load('./Data/Family_DYNET.mat');
load('./Data/Gender_DYNET.mat');

% 462 values that tag the network to which a given region belongs; there
% are 462 values because here, we have already removed the "bad areas"
% (excluded below for other variables)
load('./Data/Network.mat');

% Same for hemisphere index (462 values)
load('./Data/Hemisphere.mat');

% Labels of the 19 networks
load('./Data/Netlabs.mat');

% 462 values, this time indexing when a network has a different index in
% the left and right hemispheres
load('./Data/Nethem.mat');

% Paths to required functions and utilities
addpath('/Users/bolton/Desktop/HD_content/ET_functional');
addpath(genpath('/Users/bolton/Desktop/Utilities'));
addpath('/Users/bolton/Desktop/HD_content');
rmpath(genpath('/Users/bolton/Desktop/Utilities/spm12'));
rmpath(genpath('/Users/bolton/Desktop/Utilities/PLS_TOOLBOX'));

warning('off');



%%%%%%%%%% Second, loading of the data required for morphometric analyses.
%%%%%%%%%% We essentially only need the values for the three properties
%%%%%%%%%% studied therein (cortical thickness, surface area and mean
%%%%%%%%%% curvature). Note that this data has 29 HCs and 34 ET
%%%%%%%%%% patients, more than for resting-state fMRI
load('./Data/ET_DATA','Morpho_BASE','Morpho_BASE_SC','Morpho_YEAR',...
    'Morpho_YEAR_SC','Morpho_HC','Morpho_HC_SC');
load('./Data/ET_Clinical','TSTH_PercAllev');

% These are the matching indices of the morphometric ET patients for whom we
% also have resting-state fMRI scans. Note that here, the set of indices
% has already been trimmed to relate to only the scans with optimal quality
% both upon visual inspection and following head movement assessment. This
% yields 18 scans both pre- and post-thalamotomy
match_BASE = [2,7,11,12,15,16,17,19,20,22,23,24,25,26,27,30,32,34];
match_YEAR = [1,6,7,11,12,13,14,17,19,20,22,23,25,26,29,30,32,34];



%%%%%%%%%% Third, we need to process some variables and declare a few
%%%%%%%%%% parameters

% Which subjects do we want to run the analysis on? All baseline and 1 year
% recordings! "Subjects" contains the respective paths
Subjects = [BS,YS];

% This tags the sessions that we wish to keep because following visual
% inspection, the data was good enough. There are 41 such sessions, which
% fits with the dimensionality of the loaded data from above!
Label = [1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,1,...
    2,0,0,2,2,2,2,2,2,0,2,0,2,2,2,2,2,2,2,2,2,2,2];

% Which atlases do we use for cortex, subcortex and cerebellum? We try
% Schaefer 400 (17 net), Tian S3 and AAL
IDX_atlas = [6,9,1];

% Total number of regions that will be loaded (not the number considered)
% for the RS fMRI analysis
n_regions = 400+50+116;

% TR (in seconds)
TR = 3.3;

% Window size (in TRs); selected because TR=3.3 s and this way, we have the
% inverse of the lowest frequency present in the data (0.01 Hz)
W = 30;

% Step size (in TRs), i.e., a new computation every 6.6 seconds
Delta = 2;

% Window shape
Shape = 'Rectangular';

% Measure used to compute "functional connectivity"
Measure = 'Correlation';

% Accordingly defines the number of sessions to consider
n_base = sum(Label==1);
n_year = sum(Label==2);
n_HC = length(Subjects_HC);

% Number of different patients included in the final analyses following all
% cleaning and removal steps
n_different_patients = 23;

% Removes the subjects that were not deemed acceptable upon visual
% inspection
Subjects(Label == 0) = [];

% Main folder
Main_Folder = pwd;

% Parameters for consensus clustering
CC_distance = 'cosine';
CC_folds = 200;
CC_percdata = 0.8;
CC_type = 'items';

% Parameters for segmentation into states
n_folds_states = 1000;
max_iter_states = 200;

% Values that we explore for K (to get states)
K_range = [2:20];

% Will we plot the results?
is_plot = 1;

% Our seed is the contralateral Vim (i.e., the left Vim since our patients
% were all right-handed). Left regions follow right regions in the Tian
% atlas, and the region we select is the "THA-VPl" (ventral lateral posterior 
% nucleus of the thalamus), as its ventral subpart is known to include the
% Vim
idx_seeds = [400+31];

% Our targets are all the rest of the brain, including cortical, subcortical 
% and cerebellum. This gives a total of 475 regions
idx_targets = [1:400,...
    400+[1:30,32:50],...
    450+[91:116]];

% Note: this number is slightly trimmed down below owing to some incomplete
% fields of view (mostly in cerebellar subparts)

% How many scans to remove at the start of the recordings for magnetization
% equilibration
RM_scans = 3;

% Parameters for the removal of corrupted scans: we define the minimal
% number of samples that should be clean (FD < T_FD) for the computation of
% a windowed FC (n_min), the threshold used for this purpose (T_FD = 0.5
% mm, a standard value), the percentage of scrubbed frames above which a
% subject is excluded (P_FD = 30 %), and the minimal number of windowed FC
% values that must be available (i.e., clean) for a scan to be retained
n_min = 20;
T_FD = 0.5;
P_FD = 30;
nW_min = 85;

% These are the indices of the subjects for whom both the pre- and
% post-thalamotomy scans were retained in dFC analyses. They are used for
% plotting (Figure 2 from the manuscript)
idx_prepost(1,:) = [2,3,4,7,8,9,10,11,13,14,16,17,18];
idx_prepost(2,:) = [21,22,23,26,27,28,29,30,31,32,34,35,36];

% Numbers of regions (cortex and subcortex) considered in morphometry
% analyses
n_regions_morpho = 68;
n_regions_morpho_SC = 19;

% Number of null folds for nonparametric permutation-based significance
% assessment of individual mean and (co)variance coefficients between the
% ETpre and ETpost groups
n_KLD = 100000;

% Numbers of subjects for the morphometry data
n_HCM = 29;
n_BASE = 34;
n_YEAR = 34;

% Various colormaps used for display purposes
CM_RdBu = flipud(cbrewer('div','RdBu',1000));
CM_RdBu(CM_RdBu<0) = 0;
CM_RdBu(CM_RdBu>1) = 1;

CM_Reds = cbrewer('seq','Reds',1000);
CM_YlGn = cbrewer('seq','YlGn',1000);
CM_YlGn(CM_YlGn<0)=0;
CM_YlGn(CM_YlGn>1)=1;
CM_Greys = cbrewer('seq','Greys',1000);

CM_Groups = [[70,148,73]/255;...
    [56,61,150]/255;...
    [175,54,60]/255];

CM_Networks = cbrewer('qual','Accent',19);
CM_Networks(CM_Networks<0)=0;
CM_Networks(CM_Networks>1)=1;

CM_Patients = cbrewer('qual','Accent',24);
CM_Patients(CM_Patients<0)=0;
CM_Patients(CM_Patients>1)=1;

CM_Sim = cbrewer('seq','YlOrRd',1000);
CM_Sim(CM_Sim<0)=0;
CM_Sim(CM_Sim>1)=1;



%% Part 2. Resting-state fMRI analyses
% In this part, we wish to explore the network of regions that interact
% with the Vim seed. We propose to complement previous findings by applying 
% a dynamic functional connectivity, state-based analysis

%%%%%%%%%% A. Preprocessing: Here, (dynamic) functional connectivity is
%%%%%%%%%% quantified within the network of regions that interact with the
%%%%%%%%%% Vim seed. We consider cortical, subcortical and cerebellar
%%%%%%%%%% areas. We make sure that the outputs from this part are as clean
%%%%%%%%%% as possible when it comes to head movement, by (1) discarding
%%%%%%%%%% scans if too many frames are scrubbed, (2) discarding scans if
%%%%%%%%%% too many windows have too few clean samples for proper dFC
%%%%%%%%%% estimation. We also remove the brain regions for whom the data
%%%%%%%%%% was not consistent across all subjects (i.e., because there may
%%%%%%%%%% have been a trimmed field of view in one or the other case)

% Dynamic functional connectivity and functional connectivity are computed
% with respect to the seed region, for the ET patients and HCs
[FC_MOT,FCv_MOT,CVv_MOT,dFC_MOT,FD_MOT,Labels_scrubbing_MOT,...
    artif_FC_vc_MOT,idx_TC_toremove_MOT,idx_toscrub_MOT,W_subj_MOT] = ...
    DYNET_PrepareData_VimSeed(Main_Folder,Subjects,W,Delta,Shape,...
    Measure,n_regions,idx_seeds,idx_targets,IDX_atlas,RM_scans,T_FD,...
    P_FD,n_min,nW_min);

[FC_HC,FCv_HC,CVv_HC,dFC_HC,FD_HC,Labels_scrubbing_HC,artif_FC_vc_HC,...
    idx_TC_toremove_HC,idx_toscrub_HC,W_subj_HC] = ...
    DYNET_PrepareData_VimSeed(Main_Folder,Subjects_HC,W,Delta,Shape,...
    Measure,n_regions,idx_seeds,idx_targets,IDX_atlas,RM_scans,T_FD,...
    P_FD,n_min,1);

% Note: the final argument is set to 1 for the HC case, because we wish not
% to remove scans with less than 85 clean dFC estimates; indeed, HC samples
% are concatenated alltogether, so this criterion is then not required

% This contains the indices of the areas that we wish to discard from the
% analyses because they had NaN values in at least one subject
idx_toremove = unique([idx_TC_toremove_HC,idx_TC_toremove_MOT]);

% Note: there are 13 such areas, including the LH_LimbicA_TempPole_1, Left
% Caudate Body and Caudate Tail, and 10 cerebellar areas

% We separate the baseline and post-thalamotomy recordings, and remove the
% subjects with excessive movement (> 30% of corrupted frames or too many 
% corrupted windowed FC)
n_base_final = n_base - sum(idx_toscrub_MOT <=n_base);
n_year_final = n_year - sum(idx_toscrub_MOT >n_base);
n_HC_final = n_HC - length(idx_toscrub_HC);

% Note: the numbers of subjects in each group, following preprocessing, are
% not exactly the ones reported in the manuscript

FC_HC(:,:,idx_toscrub_HC) = [];
FC_MOT(:,:,idx_toscrub_MOT) = [];

FCv_HC(:,idx_toscrub_HC) = [];
FCv_MOT(:,idx_toscrub_MOT) = [];

CVv_HC(:,idx_toscrub_HC) = [];
CVv_MOT(:,idx_toscrub_MOT) = [];

dFC_HC(idx_toscrub_HC) = [];
dFC_MOT(idx_toscrub_MOT) = [];

Labels_scrubbing_HC(idx_toscrub_HC) = [];
Labels_scrubbing_MOT(idx_toscrub_MOT) = [];

FC_BASE = FC_MOT(:,:,1:n_base_final);
FC_YEAR = FC_MOT(:,:,n_base_final+1:end);

FCv_BASE = FCv_MOT(:,1:n_base_final);
FCv_YEAR = FCv_MOT(:,n_base_final+1:end);

CVv_BASE = CVv_MOT(:,1:n_base_final);
CVv_YEAR = CVv_MOT(:,n_base_final+1:end);

dFC_BASE = dFC_MOT(1:n_base_final);
dFC_YEAR = dFC_MOT(n_base_final+1:end);

Labels_scrubbing_BASE = Labels_scrubbing_MOT(1:n_base_final);
Labels_scrubbing_YEAR = Labels_scrubbing_MOT(n_base_final+1:end);

% All these variables have size 36 (18 pre and 18 post values)
TSTH_DYNET(idx_toscrub_MOT) = [];
TSTHP_DYNET(idx_toscrub_MOT) = [];
Age_DYNET(idx_toscrub_MOT) = [];
Family_DYNET(idx_toscrub_MOT) = [];
Gender_DYNET(idx_toscrub_MOT) = [];
LesionVolume_DYNET(idx_toscrub_MOT) = [];
SymptomsDuration_DYNET(idx_toscrub_MOT) = [];
TremorStop_DYNET(idx_toscrub_MOT) = [];

% We discard the time courses that should be discarded
M_torem = zeros(n_regions,n_regions);
M_torem(idx_toremove,:) = 1;
M_torem(:,idx_toremove) = 1;

tmp = jUpperTriMatToVec(M_torem);
idx_toremove_vec = logical(tmp(artif_FC_vc_MOT));

% We obtain the final FC matrices
FC_HC(idx_toremove,:,:) = [];
FC_HC(:,idx_toremove,:) = [];

FC_BASE(idx_toremove,:,:) = [];
FC_BASE(:,idx_toremove,:) = [];

FC_YEAR(idx_toremove,:,:) = [];
FC_YEAR(:,idx_toremove,:) = [];

FCv_HC(idx_toremove_vec,:) = [];
FCv_BASE(idx_toremove_vec,:) = [];
FCv_YEAR(idx_toremove_vec,:) = [];

CVv_HC(idx_toremove_vec,:) = [];
CVv_BASE(idx_toremove_vec,:) = [];
CVv_YEAR(idx_toremove_vec,:) = [];

% Similarly, we remove the excess time courses for dFC data
for s = 1:length(dFC_HC)
    dFC_HC{s}(idx_toremove_vec,:) = [];
end

for s = 1:length(dFC_BASE)
    dFC_BASE{s}(idx_toremove_vec,:) = [];
end

for s = 1:length(dFC_YEAR)
    dFC_YEAR{s}(idx_toremove_vec,:) = [];
end

% Also trims the atlas labels
Labels_DYNET_finalized(idx_toremove) = [];
Labels_DYNET1(idx_toremove_vec) = [];
Labels_DYNET2(idx_toremove_vec) = [];

% Trims the list of subjects to reflect only the final, clean ones
Subjects(idx_toscrub_MOT) = [];



%%%%%%%%%% B. Data concatenation: Here, we create the data of interest for
%%%%%%%%%% our dynamic analyses. On the one hand, we add all clean windowed
%%%%%%%%%% FC estimates available across all HCs, in order to subsequently
%%%%%%%%%% estimates dFC states. On the other hand, we also concatenate the
%%%%%%%%%% patients data in a separate array.

% This will contain the data of interest
Data_BASE = [];
Data_HC = [];
Data_YEAR = [];

Labels_HC_time = [];
Labels_HC_subjects = [];

% First, we consider the HCs recordings. We will consider all the windows
% that are "clear", even if it results in some subjects contributing more
% than others. This is because we make the assumption that compared to the
% patients, the pool of HCs is somehow "homogeneous"
for s = 1:n_HC_final

    s

    for t = 1:size(dFC_HC{s},2)

        % Creates the equivalent compatible with the morphometric atlas
        % results
        if Labels_scrubbing_HC{s}(t) == 0
            Data_HC = [Data_HC,dFC_HC{s}(:,t)];
            Labels_HC_time = [Labels_HC_time,t];
            Labels_HC_subjects = [Labels_HC_subjects,s];
        end
    end
end

% Note: Data_HC contains a total of 1136 samples, each of dimension 462

% Then, we consider the baseline and post-intervention recordings; we only
% take the samples that are clean, and if a subject does not have enough,
% it is excluded
for s = 1:n_base_final

    s

    tmp = dFC_BASE{s}(:,Labels_scrubbing_BASE{s}==0);
    Data_BASE = [Data_BASE,tmp(:,1:nW_min)];
end

for s = 1:n_year_final

    s

    tmp = dFC_YEAR{s}(:,Labels_scrubbing_YEAR{s}==0);
    Data_YEAR = [Data_YEAR,tmp(:,1:nW_min)];
end

% Note: In a case where there are more than nW_min=85 samples estimated for
% a given scan, we always select the first available ones. By this mean, we
% avoid effects due to drowsiness in the scanner as much as possible



%%%%%%%%%% C. Extraction of dynamic functional connectivity metrics: First,
%%%%%%%%%% we determine the optimal number of states in healthy
%%%%%%%%%% individuals. Second, we assess the similarity of the ET patients
%%%%%%%%%% frames to the HC states.

% First, we consider the HC data and determine how many states there are.
[Consensus_HC] = CAP_ConsensusClustering({Data_HC},K_range,CC_type,CC_percdata,CC_folds,CC_distance);
[~,PAC_HC] = ComputeClusteringQuality(Consensus_HC,K_range);

if is_plot
   figure;
   bar(K_range,PAC_HC);
end

% Our results point towards K=3 different patterns
K_opt_HC = 3;

% Using this optimum, we derive the states themselves
idx_HC = kmeans(Data_HC',K_opt_HC,'distance',CC_distance,'replicates',...
    n_folds_states,'empty','drop','maxiter',max_iter_states);

% Computes the centroids as the average of frames clustered together
for st = 1:K_opt_HC
    Centroids_HC(st,:) = mean(Data_HC(:,idx_HC==st),2);
end

% Computes the distribution of similarity values of frames to their
% centroids
for f = 1:size(Data_HC,2)
    Sim_toCentroids(f) = corr(Data_HC(:,f),Centroids_HC(idx_HC(f),:)');
end

% Transitions and summary of counts
for s = 1:n_HC_final
    n_frames_HC(s) = sum(Labels_HC_subjects==s);
    
    for st = 1:K_opt_HC
        tmp = idx_HC(Labels_HC_subjects==s);
        Counts_HC(s,st) = sum(tmp==st)/n_frames_HC(s);
    end
end

% Computation of similarity values for each scan
idx_start = 1;

% Note that we can run one loop only because we have the same number of
% BASE and YEAR scans, coincidentally
for s = 1:n_base_final
    
    % We sample the data from the subject at hand
    tmp_data = Data_BASE(:,idx_start:idx_start+nW_min-1);
    tmp_data2 = Data_YEAR(:,idx_start:idx_start+nW_min-1);
    
    % Similarity to the HC states, and indices of the most similar state to
    % the considered frame, for both pre and post-thalamotomy scans
    Sim_BASE(:,:,s) = corr(tmp_data,Centroids_HC');
    [SimMax_BASE(:,s),IndMax_BASE(:,s)] = max(squeeze(Sim_BASE(:,:,s)),[],2);
    
    Sim_YEAR(:,:,s) = corr(tmp_data2,Centroids_HC');
    [SimMax_YEAR(:,s),IndMax_YEAR(:,s)] = max(squeeze(Sim_YEAR(:,:,s)),[],2);
    
    % We also compute the similarity separately for each state. We extract
    % three features: mean similarity, std of similarity, and coefficient
    % of variation (that is, std/mean)
    for st = 1:K_opt_HC
        MeanSimStates_BASE(st,s) = mean(SimMax_BASE(IndMax_BASE(:,s)==st,s));
        StdSimStates_BASE(st,s) = std(SimMax_BASE(IndMax_BASE(:,s)==st,s));
        
        MeanSimStates_YEAR(st,s) = mean(SimMax_YEAR(IndMax_YEAR(:,s)==st,s));
        StdSimStates_YEAR(st,s) = std(SimMax_YEAR(IndMax_YEAR(:,s)==st,s));
        
        CoefVarStates_BASE(st,s) = StdSimStates_BASE(st,s)/MeanSimStates_BASE(st,s);
        CoefVarStates_YEAR(st,s) = StdSimStates_YEAR(st,s)/MeanSimStates_YEAR(st,s);
    end
    
    idx_start = idx_start + nW_min;
end

% Computation of the counts
for st = 1:K_opt_HC
    Counts_BASE(st,:) = sum(IndMax_BASE == st)/nW_min;
    Counts_YEAR(st,:) = sum(IndMax_YEAR == st)/nW_min;
end

% Here, we compute similarity regardless of the state
MeanSim_BASE = mean(SimMax_BASE);
StdSim_BASE = std(SimMax_BASE);

MeanSim_YEAR = mean(SimMax_YEAR);
StdSim_YEAR = std(SimMax_YEAR);

CoefVar_BASE = StdSim_BASE./MeanSim_BASE;
CoefVar_YEAR = StdSim_YEAR./MeanSim_YEAR;



%%%%%%%%%% D. Statistical analysis: to investigate potential group
%%%%%%%%%% differences, we resort to mixed modelling. To assess clinical
%%%%%%%%%% predictive potential, we use a GLM design with lesion volume as
%%%%%%%%%% an interacting factor

% Data for mixed modelling, that is, subject and group variables
% (respectively featuring as random and fixed effect terms in the model)
ind_MM = [(1:n_base_final),19,20,2,3,4,21,22,7,8,9,10,11,13,14,23,16,17,18]';
Subject_MM = nominal([((n_different_patients+1:n_different_patients+n_HC_final)');ind_MM]);
Group_MM = nominal([ones(n_HC_final,1);2*ones(n_base_final,1);3*ones(n_year_final,1)]);

% Note: the subject indices were manually checked. When there is a new
% index after the first 18 indices, it means that there is a new subject in
% the ETpost group that was not included in the ETpre group

%%% We conduct the statistical analysis on dynamic functional connectivity 
%%% metrics, that is, counts, mean similarity and std of similarity to a 
%%% given state, plus coefficient of variation
for st = 1:K_opt_HC
    p_States(st,:) = anovan([Counts_BASE(st,:),Counts_YEAR(st,:)]',{(Group_MM(n_HC_final+1:end)),Subject_MM(n_HC_final+1:end)},'random',2,'display','off');
    [p_MeanSimSt(st,:),ttttt{st}] = anovan([MeanSimStates_BASE(st,:),MeanSimStates_YEAR(st,:)]',{(Group_MM(n_HC_final+1:end)),Subject_MM(n_HC_final+1:end)},'random',2,'display','off');
    [p_StdSimSt(st,:),ttt{st}] = anovan([StdSimStates_BASE(st,:),StdSimStates_YEAR(st,:)]',{(Group_MM(n_HC_final+1:end)),Subject_MM(n_HC_final+1:end)},'random',2,'display','off');
    [p_CoefVarSt(st,:),tt{st}] = anovan([CoefVarStates_BASE(st,:),CoefVarStates_YEAR(st,:)]',{(Group_MM(n_HC_final+1:end)),Subject_MM(n_HC_final+1:end)},'random',2,'display','off');
end

% Note: although we included HCs in the original data fed to anovan, we
% only perform an ETpost versus ETpre comparison in the end, because the
% HCs were already included in the analysis since we consider similarity to
% HC states

% Correction for multiple comparison
p_States = p_States * K_opt_HC;
p_MeanSimSt = p_MeanSimSt * K_opt_HC;
p_StdSimSt = p_StdSimSt * K_opt_HC;
p_CoefVarSt = p_CoefVarSt * K_opt_HC;



%%% Now, we look at links to clinical recovery

% Definition of the required quantities
Age_MM = Age_DYNET(1:n_base_final);
Gender_MM = nominal(Gender_DYNET(1:n_base_final));
Subject_MM = nominal((1:n_base_final)');
LesionVolume_MM = LesionVolume_DYNET(1:n_base_final);
Outcome_MM = TSTH_DYNET(1:n_base_final)-TSTHP_DYNET(1:n_base_final);

for st = 1:K_opt_HC
    
    % Counts
    Metric_MM = Counts_BASE(st,:)';
    TBL = table(Metric_MM,Age_MM,Gender_MM,Subject_MM,LesionVolume_MM,Outcome_MM);
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM +  Age_MM + Gender_MM');
    ClinSign_Counts{st} = glm.Coefficients;
    
    % MeanSim
    Metric_MM = MeanSimStates_BASE(st,:)';
    TBL = table(Metric_MM,Age_MM,Gender_MM,Subject_MM,LesionVolume_MM,Outcome_MM);
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM +  Age_MM + Gender_MM');
    ClinSign_MeanSim{st} = glm.Coefficients;
    
    % StdSim
    Metric_MM = StdSimStates_BASE(st,:)';
    TBL = table(Metric_MM,Age_MM,Gender_MM,Subject_MM,LesionVolume_MM,Outcome_MM);
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM +  Age_MM + Gender_MM');
    ClinSign_StdSim{st} = glm.Coefficients;
    
    % CoefVar
    Metric_MM = CoefVarStates_BASE(st,:)';
    TBL = table(Metric_MM,Age_MM,Gender_MM,Subject_MM,LesionVolume_MM,Outcome_MM);
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM +  Age_MM + Gender_MM');
    ClinSign_CoefVar{st} = glm.Coefficients;
end

% Balance between both states with significance
Metric_MM = Counts_BASE(3,:)'-Counts_BASE(1,:)';
TBL = table(Metric_MM,Age_MM,Gender_MM,Subject_MM,LesionVolume_MM,Outcome_MM);
glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM +  Age_MM + Gender_MM');
ClinSign_StateBalance = glm.Coefficients;



%%%%%%%%%% E. Visualizations for the above results


%%%%%% V1: Visualization of the three dFC states from HCs
DYNET_PlotStateGraph_VimSeed(Centroids_HC(1,:)',0,[],CM_Networks,Network,Hemisphere);
DYNET_PlotStateGraph_VimSeed(Centroids_HC(2,:)',0,[],CM_Networks,Network,Hemisphere);
DYNET_PlotStateGraph_VimSeed(Centroids_HC(3,:)',0,[],CM_Networks,Network,Hemisphere);

% Summary of connectivity in each of the three states across networks
% in terms of mean, std and coefficient of variation
for n = 1:max(Nethem)
    MeanCentroids_HC(:,n) = mean(Centroids_HC(:,Nethem==n),2);
    StdCentroids_HC(:,n) = std(Centroids_HC(:,Nethem==n),[],2);
end

ZScore_Centroids_HC = MeanCentroids_HC./StdCentroids_HC;

figure; 
imagesc(MeanCentroids_HC);
colormap(CM_RdBu);
caxis([-0.4,0.4]);

figure; 
imagesc(StdCentroids_HC);
colormap(cbrewer('seq','Reds',1000));
caxis([0,0.4]);

figure; 
imagesc(ZScore_Centroids_HC);
colormap(CM_RdBu);
caxis([-5,5]);


%%%%%% V3: Visualization as surface plots
CMR = cbrewer('seq','Reds',10000);
CMB = cbrewer('seq','Blues',10000);
CMR(CMR<0)=0;
CMR(CMR>1)=1;
CMB(CMB<0)=0;
CMB(CMB>1)=1;

Weights_Signs = sign(Centroids_HC(3,:)');
Weights = abs(Centroids_HC(3,:)');
Weights = (Weights)/max(Weights);

for e = 1:462
    if Weights_Signs(e) >= 0
        tmp_col2(e,:) = CMR(floor(Weights(e)*10000),:);
    elseif Weights_Signs(e) < 0
        tmp_col2(e,:) = CMB(floor(Weights(e)*10000),:);
    end
end

addpath(genpath('/Users/bolton/Desktop/Utilities/spm12'));

brain_info = spm_vol('Schaefer2018_17N_400R.nii');
SCH_values = spm_read_vols(brain_info);

Weights_surface(1:113) = Weights(1:113).*Weights_Signs(1:113);

% I know that only area 114 was excluded originally owing to FoV
Weights_surface(114) = 0;
Weights_surface(115:400) = Weights(114:399).*Weights_Signs(114:399);

voxel_size = diag(brain_info.mat);
voxel_size = voxel_size(1:end-1)';

voxel_shift = brain_info.mat(:,4);
voxel_shift = voxel_shift(1:end-1)';

V_surface = zeros(brain_info.dim);

for r = 1:400
    V_surface(SCH_values==r) = Weights_surface(r);
end

tmp_NIFTI = make_nii(V_surface,voxel_size,-voxel_shift./voxel_size);
tmp_NIFTI.hdr.dime.datatype=64;
tmp_NIFTI.hdr.dime.bitpix=64;

save_nii(tmp_NIFTI,'Figure_Surf_State3.nii');

rmpath(genpath('/Users/bolton/Desktop/Utilities/spm12'));

%%%%%% V5: Visualization of pre- and post-thalamotomy results for each of
%%%%%% the metrics of interest: counts, mean similarity, std(similarity)
%%%%%% and coefficient of variation. The plots show bars and SEM, as well
%%%%%% as individual data points, connected when belonging to the same
%%%%%% subject before and after intervention

%%%% A. Counts for state 1
figure;
hold on

Cat1 = Counts_BASE(1,:);
Cat2 = Counts_YEAR(1,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.9]);


%%%% B. Counts for state 2
figure;
hold on

% Plot of pre and post values for similarity to state 3
Cat1 = Counts_BASE(2,:);
Cat2 = Counts_YEAR(2,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.9]);


%%%% C. Counts for state 3
figure;
hold on

% Same for std similarity for state 3
Cat1 = Counts_BASE(3,:);
Cat2 = Counts_YEAR(3,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.9]);




%%%% D. Mean similarity for state 1
figure;
hold on

Cat1 = MeanSimStates_BASE(1,:);
Cat2 = MeanSimStates_YEAR(1,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.6]);


%%%% E. Mean similarity for state 2
figure;
hold on

% Plot of pre and post values for similarity to state 3
Cat1 = MeanSimStates_BASE(2,:);
Cat2 = MeanSimStates_YEAR(2,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.6]);


%%%% F. Mean similarity for state 3
figure;
hold on

Cat1 = MeanSimStates_BASE(3,:);
Cat2 = MeanSimStates_YEAR(3,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.6]);




%%%% G. Std similarity (state 1)
figure;
hold on

Cat1 = StdSimStates_BASE(1,:);
Cat2 = StdSimStates_YEAR(1,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.15]);

%%%% H. Std similarity (state 2)
figure;
hold on

Cat1 = StdSimStates_BASE(2,:);
Cat2 = StdSimStates_YEAR(2,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.15]);

%%%% I. Std similarity (state 3)
figure;
hold on

Cat1 = StdSimStates_BASE(3,:);
Cat2 = StdSimStates_YEAR(3,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.15]);

%%%% J. Coef var (state 1)
figure;
hold on

Cat1 = CoefVarStates_BASE(1,:);
Cat2 = CoefVarStates_YEAR(1,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.5]);


%%%% K. Coef var (state 2)
figure;
hold on

Cat1 = CoefVarStates_BASE(2,:);
Cat2 = CoefVarStates_YEAR(2,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.5]);


%%%% L. Coef var (state 3)
figure;
hold on

Cat1 = CoefVarStates_BASE(3,:);
Cat2 = CoefVarStates_YEAR(3,:);

bar([nanmean(Cat1),nanmean(Cat2)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
errorbar([nanmean(Cat1),nanmean(Cat2)],[nanstd(Cat1),nanstd(Cat2)]/sqrt(length(Cat1)),'LineStyle','none','LineWidth',2,'Color','k');

for s = 1:18
    if s == 9 || s == 11 || s == 13
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled');
    else
        scatter(1.1,Cat1(s),60,CM_Patients(ind_MM(s),:),'filled','s');
    end
    
    if s == 28-18 || s == 30-18 || s == 31-18
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled');
    else
        scatter(2.1,Cat2(s),60,CM_Patients(ind_MM(s+18),:),'filled','s');
    end
end

imm_pre = ind_MM(1:18);
imm_post = ind_MM(19:end);

for s = 1:24
    if sum(ind_MM == s) > 1
        plot([1.1,2.1],[Cat1(imm_pre==s),Cat2(imm_post==s)],'LineStyle','--','LineWidth',0.5,'Color',CM_Patients(s,:));
    end
end

xlim([0.5,2.5]);
ylim([0,0.5]);


%%%%%% V6: Visualization of similarity as a function of time, and state
%%%%%% identity across subjects
figure;
subplot(2,1,1);
imagesc(SimMax_BASE);
caxis([-0.1,0.7]);
colormap(CM_Sim);

subplot(2,1,2);
imagesc(SimMax_YEAR);
caxis([-0.1,0.7]);
colormap(CM_Sim);

figure;
subplot(2,1,1);
imagesc(IndMax_BASE);
colormap(CM_Groups);

subplot(2,1,2);
imagesc(IndMax_YEAR);
colormap(CM_Groups);

for st = 1:K_opt_HC
    
    tmp2 = SimMax_BASE(IndMax_BASE == st);
    tmp3 = SimMax_YEAR(IndMax_YEAR == st);
    tmp2 = tmp2(:);
    tmp3 = tmp3(:);

    figure;
    hold on

    h2 = histogram(tmp2,40);
    set(h2,'FaceColor',[56,61,150]/255,'EdgeColor','None','FaceAlpha',0.4);

    h3 = histogram(tmp3,40);
    set(h3,'FaceColor',[175,54,60]/255,'EdgeColor','None','FaceAlpha',0.4);
    
    xlim([-0.1,0.7]);
    ylim([0,60]);
end


%%%%%% V7: Visualization of plots for links to clinical improvement
figure;
scatter(CoefVarStates_BASE(2,:),TSTH_DYNET(1:18)-TSTHP_DYNET(1:18),40,'filled','k');

figure;
scatter(Counts_BASE(3,:)'-Counts_BASE(1,:)',TSTH_DYNET(1:18)-TSTHP_DYNET(1:18),LesionVolume_DYNET(1:18)*300,'filled','k');



%% Part 3. Morphometry analyses
% In this part, we explore to what extent the morphometric data from ET
% patients changes upon thalamotomy, leveraging a recently introduced
% analytical approach of ours.

%%%%%%%%%% A. Assessment of similarity to HCs via log-likelihood: Here, the
%%%%%%%%%% similarity to the HC distribution is assessed

for r = 1:n_regions_morpho

    % Data from the three groups
    tmp_data = squeeze(Morpho_HC(r,:,:));
    tmp_data2 = squeeze(Morpho_BASE(r,:,:));
    tmp_data3 = squeeze(Morpho_YEAR(r,:,:));
    
    % Covariance and mean estimates in the HC population (maximum
    % likelihood estimates)
    COV_HC(r,:,:) = cov(tmp_data);
    MU_HC(r,:) = mean(tmp_data);
    
    % For control subjects, likelihood is computed for each region
    [L_HC_toHC(r,:),LL_HC_toHC(r,:)] = ET_EvaluateGaussian(tmp_data',MU_HC(r,:)',squeeze(COV_HC(r,:,:)));
    
    % Likelihood is also computed for the ET subjects, before and
    % after intervention, when assessing the probability that they could be
    % obtained from the HC distribution
    [L_BASE_toHC(r,:),LL_BASE_toHC(r,:)] = ET_EvaluateGaussian(tmp_data2',MU_HC(r,:)',squeeze(COV_HC(r,:,:)));
    [L_YEAR_toHC(r,:),LL_YEAR_toHC(r,:)] = ET_EvaluateGaussian(tmp_data3',MU_HC(r,:)',squeeze(COV_HC(r,:,:)));
end

% Same on the non-cortical areas (only univariate)
for r = 1:n_regions_morpho_SC
    
    % HC data (size S x M)
    tmp_data = squeeze(Morpho_HC_SC(r,:));
    tmp_data2 = squeeze(Morpho_BASE_SC(r,:));
    tmp_data3 = squeeze(Morpho_YEAR_SC(r,:));

    % Covariance and mean estimates in the control population (maximum
    % likelihood estimates)
    COV_SC_HC(r) = cov(tmp_data);
    MU_SC_HC(r) = mean(tmp_data);
    
    % For control subjects, likelihood is computed for each region
    [L_SC_HC_toHC(r,:),LL_SC_HC_toHC(r,:)] = ET_EvaluateGaussian(tmp_data,MU_SC_HC(r),squeeze(COV_SC_HC(r)));
    [L_SC_BASE_toHC(r,:),LL_SC_BASE_toHC(r,:)] = ET_EvaluateGaussian(tmp_data2,MU_SC_HC(r),squeeze(COV_SC_HC(r)));
    [L_SC_YEAR_toHC(r,:),LL_SC_YEAR_toHC(r,:)] = ET_EvaluateGaussian(tmp_data3,MU_SC_HC(r),squeeze(COV_SC_HC(r)));
end

% Finds the regions for which likelihood to the HC distribution differs
% between the baseline and 1-year ET populations
for r = 1:n_regions_morpho
    [p_LL(r),~,stats_LL{r}] = ranksum(LL_BASE_toHC(r,:)',LL_YEAR_toHC(r,:)');
end

for r = 1:n_regions_morpho_SC
    [p_LL_SC(r),~,stats_LL_SC{r}] = ranksum(LL_SC_BASE_toHC(r,:)',LL_SC_YEAR_toHC(r,:)');
end

% Bonferroni correction is applied for the P=87 tests run in parallel
p_LL = p_LL*(n_regions_morpho+n_regions_morpho_SC);
p_LL_SC = p_LL_SC*(n_regions_morpho+n_regions_morpho_SC);

idx_LL_toHC = find(p_LL < 0.05);
idx_LL_toHC_SC = find(p_LL_SC < 0.05);



%%%%%%%%%% B. Assessment of post- vs pre-thalamotomy group differences in
%%%%%%%%%% terms of morphometric properties: Here, we consider cortical
%%%%%%%%%% thickness, surface area and mean curvature together to yield,
%%%%%%%%%% for cortical regions, a three-dimensional data point per
%%%%%%%%%% subject. We compare the multivariate Gaussian representations of
%%%%%%%%%% patients before and after intervention, for each parameter, to
%%%%%%%%%% assess whether there are significant differences

% Actual computations of group deltas for cortical and non-cortical regions
for r = idx_LL_toHC

    r

    % Our reference population is the set of baseline ET patients
    X_ref = squeeze(Morpho_BASE(r,:,:));
    
    % Our other population is the set of 1-year ET patients
    X = squeeze(Morpho_YEAR(r,:,:));
     
    Results_Morpho_GD{r} = DYNET_DetermineSubcases(X,X_ref,n_KLD,NaN);
end

for r = idx_LL_toHC_SC
    
    % Our reference population is the set of baseline ET patients
    X_ref = squeeze(Morpho_BASE_SC(r,:));
    
    % Our other population is the set of 1-year ET patients
    X = squeeze(Morpho_YEAR_SC(r,:));
    
    [Results_Morpho_SC_GD{r}] = DYNET_DetermineSubcases_1D(X,X_ref,n_KLD,NaN);
end

% Establishes the significant regions for each coefficient
for r = idx_LL_toHC
    [DM(:,r),DS(:,:,r),p_DS(:,:,r)] = DYNET_AssessSignificance(Results_Morpho_GD{r},103);
end

for r = idx_LL_toHC_SC
    [DM_SC(r),DS_SC(r),p_DS_SC(r)] = DYNET_AssessSignificance_SC(Results_Morpho_SC_GD{r},103);
end

% Note: "103" is a hard-coded value that stands for the number of parallel
% tests done (based on the regions that survived significance above, and
% the dimensionality of their morphometric feature space)



%%%%%%%%%% C. Statistical analyses: Here, we perform the same statistical
%%%%%%%%%% analyses as in the previous section in terms of possible
%%%%%%%%%% clinical utility of the pre-intervention similarity to HC

% From the morphometric side, we sample the pre-thalamotomy scans
% for which we also have the RS-fMRI data
LL_pre_toHC = [LL_BASE_toHC(:,match_BASE)];
LL_pre_toHC_SC = [LL_SC_BASE_toHC(:,match_BASE)];

% Stats per se
for r = idx_LL_toHC
    
    Subject_MM = nominal((1:n_base_final)');
    LesionVolume_MM = LesionVolume_DYNET(1:n_base_final);
    Outcome_MM = TSTH_DYNET(1:n_base_final)-TSTHP_DYNET(1:n_base_final);

    Metric_MM = LL_pre_toHC(r,:)';

    % All summarized in a table to enable mixed modelling...
    TBL = table(Metric_MM,Subject_MM,LesionVolume_MM,Outcome_MM);

    % Mixed model itself
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM');

    ClinSign_Morpho{r} = glm.Coefficients;
end

for r = idx_LL_toHC_SC
    
    Subject_MM = nominal((1:n_base_final)');
    LesionVolume_MM = LesionVolume_DYNET(1:n_base_final);
    Outcome_MM = TSTH_DYNET(1:n_base_final)-TSTHP_DYNET(1:n_base_final);

    Metric_MM = LL_pre_toHC_SC(r,:)';

    % All summarized in a table to enable mixed modelling...
    TBL = table(Metric_MM,Subject_MM,LesionVolume_MM,Outcome_MM);

    % Mixed model itself
    glm = fitglm(TBL,'Outcome_MM ~ 1 + LesionVolume_MM * Metric_MM');

    ClinSign_Morpho_SC{r} = glm.Coefficients;
end

% Note: in the above assessments, age and gender were not included because
% they had already been regressed out during processing!



%%%%%%%%%% D. Visualization of the results

if is_plot
    
    %%% A. Examination of the conformity to HCs before and after
    %%% intervention
    
    % Summary of the data: log-likelihood with respect to the HC
    % distribution, and "distance" to the baseline ET distribution between
    % the individual baseline and 1-year measurements
    figure;
    plot(mean(LL_HC_toHC,2),'color',[70,148,73]/255,'LineWidth',2);
    hold on
    plot(mean(LL_BASE_toHC,2),'color',[56,61,150]/255,'LineWidth',2);
    plot(mean(LL_YEAR_toHC,2),'color',[175,54,60]/255,'LineWidth',2);
    errorbar(mean(LL_HC_toHC,2),std(LL_HC_toHC,[],2)/sqrt(n_HCM),'color','k');
    errorbar(mean(LL_BASE_toHC,2),std(LL_BASE_toHC,[],2)/sqrt(n_BASE),'color','k','LineWidth',0.5);
    errorbar(mean(LL_YEAR_toHC,2),std(LL_YEAR_toHC,[],2)/sqrt(n_YEAR),'color','k','LineWidth',0.5);
    set(gca,'Box','off');
    xlabel('Region');
    ylabel('Log-likelihood to HC');
    
    
    % The same subcortically
    figure;
    plot(mean(LL_SC_HC_toHC,2),'color',[70,148,73]/255,'LineWidth',2);
    hold on
    plot(mean(LL_SC_BASE_toHC,2),'color',[56,61,150]/255,'LineWidth',2);
    plot(mean(LL_SC_YEAR_toHC,2),'color',[175,54,60]/255,'LineWidth',2);
    errorbar(mean(LL_SC_HC_toHC,2),std(LL_SC_HC_toHC,[],2)/sqrt(n_HCM),'color','k','LineWidth',0.5);
    errorbar(mean(LL_SC_BASE_toHC,2),std(LL_SC_BASE_toHC,[],2)/sqrt(n_BASE),'color','k','LineWidth',0.5);
    errorbar(mean(LL_SC_YEAR_toHC,2),std(LL_SC_YEAR_toHC,[],2)/sqrt(n_YEAR),'color','k','LineWidth',0.5);
    set(gca,'Box','off');
    xlabel('Region');
    ylabel('Log-likelihood to HC');
    
    % Individual plots of the distributions for the regions showing group
    % differences
    
%     for r = idx_LL_toHC
%         ET_PlotEllipses(Morpho_HC,Morpho_BASE,Morpho_YEAR,r,0.682,TSTH_ForPlot);
%     end
end



%% Part 4. Links between morphometric and resting-state data
% Are the features computed from RS-fMRI data somehow related to the
% morphometric ones?

clear Sim_BASE
clear Sim_YEAR

%%%%%%%%%% A. Computations

LL_before = [LL_BASE_toHC(:,match_BASE);LL_SC_BASE_toHC(:,match_BASE)];
LL_after = [LL_YEAR_toHC(:,match_BASE);LL_SC_YEAR_toHC(:,match_BASE)];

% Actual similarity and null distribution built upon permutation assessment
for k = 1:K_opt_HC
    
    % We determine which entries are NaN for the pre- and post- cases (when
    % a state is not expressed by a subject)
    tmp1 = MeanSimStates_BASE(k,:);
    tmp2 = MeanSimStates_YEAR(k,:);
    
    idx_nan1 = (isnan(tmp1));
    idx_nan2 = (isnan(tmp2));
    
    % We compute correlation with the morpho LLs in both cases
    Sim_BASE(:,k) = corr(LL_before(:,~idx_nan1)',MeanSimStates_BASE(k,~idx_nan1)');
    Sim_YEAR(:,k) = corr(LL_after(:,~idx_nan2)',MeanSimStates_YEAR(k,~idx_nan2)');
    
    % Our statistic of interest is the difference
    DeltaSim(:,k) = Sim_YEAR(:,k) - Sim_BASE(:,k);
    
    % Same for null data
    for n = 1:50000
    
        tmp3 = LL_before(:,~idx_nan1)';
        tmp4 = LL_after(:,~idx_nan2)';
        
        tmp5 = MeanSimStates_BASE(k,~idx_nan1)';
        tmp6 = MeanSimStates_YEAR(k,~idx_nan2)';
        
        idx_1 = randperm(size(tmp3,1));
        idx_2 = randperm(size(tmp4,1));

        Sim_BASE_null(:,k,n) = corr(tmp3,tmp5(idx_1));
        Sim_YEAR_null(:,k,n) = corr(tmp4,tmp6(idx_2));
        
        DeltaSim_null(:,k,n) = Sim_YEAR_null(:,k,n) - Sim_BASE_null(:,k,n);
    end
end

for k = 1:K_opt_HC
    [p_Sim(k),~,stats_sim{k}] = ranksum(Sim_BASE(:,k),Sim_YEAR(:,k));
end



%%%%%%%%%% B. Visualizations

figure;
subplot(1,3,1);
imagesc(Sim_BASE(:,[3,1,2]));
caxis([-0.7,0.7]);
axis off
colormap(CM_RdBu);

subplot(1,3,2);
imagesc(Sim_YEAR(:,[3,1,2]));
caxis([-0.7,0.7]);
axis off
colormap(CM_RdBu);

subplot(1,3,3);
imagesc(DeltaSim(:,[3,1,2]));
caxis([-1,1]);
axis off
colormap(CM_RdBu);

UpperT = squeeze(prctile(DeltaSim_null,100-2.5,3));
LowerT = squeeze(prctile(DeltaSim_null,2.5,3));

Sign_DeltaSim = (DeltaSim > UpperT | DeltaSim < LowerT);

for k = 1:K_opt_HC
    
    tmp2 = Sim_BASE(:,k);
    tmp3 = Sim_YEAR(:,k);

    figure;
    hold on

    h2 = histogram(tmp2,40);
    set(h2,'FaceColor',[56,61,150]/255,'EdgeColor','None','FaceAlpha',0.4);

    h3 = histogram(tmp3,40);
    set(h3,'FaceColor',[175,54,60]/255,'EdgeColor','None','FaceAlpha',0.4);
    
    xlim([-0.7,0.7]);
    ylim([0,9]);
end