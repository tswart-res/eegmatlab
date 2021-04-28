%% Run within comparions for each group

% clear all; clc; close all

if ~exist ( 'ft_EEG', 'var' ), load('C:\MATLAB\EEGData\ASMR-R\FT\ft_EEG_RS_EO.mat'), end

grouplist = {'ASMRR', 'Controls'}

for i = 1:length(grouplist)
%% Frequency analysis
Group1 = grouplist{i}
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
% cfg.foi = [1:80];
cfg.foilim = [1 80];
cfg.taper = 'dpss';
cfg.tapsmofrq = [2];
cfg.keeptrials = 'no'; % grand average wants this as no.

statelist = fieldnamesr(ft_EEG.(Group1), 1)         % genrates a list of structure files with master structure (state structure files here)

for i = 1:length(statelist)
    structstate = statelist{i}
    structpptlist = fieldnamesr(ft_EEG.(Group1).(structstate), 1) % generates list of substructure files in statestructure (ppt)
    if ~exist('ft_tf_EEG','var') || ~isfield(ft_tf_EEG,Group1) || ~isfield(ft_tf_EEG.(Group1),structstate)
        for o = 1:length(structpptlist)
            structppt = structpptlist{o}
            ft_tf_EEG.(Group1).(structstate).(structppt) = ft_freqanalysis(cfg,ft_EEG.(Group1).(structstate).(structppt))
        end
    else
    end
end


%% Select only participants who have data for both experimental and baseline conditions

% clearvars -except ft_EEG ft_tf_EEG
% close all
    statelist = fields(ft_EEG.(Group1))
%     ba_vs_exp_tf = []

% statelist =
% 
%   2×1 cell array
% 
%     {'EC_start'}
%     {'EC_end'  }


    BaselineState = 1
    ExperimentalState = 2
    
    statedescript = {'PostASMR\_Baseline', 'preASMR\_Baseline', 'postASMR\_Relaxed', 'preASMR\_Relaxed', 'WeakASMR', 'StrongASMR', 'preASMR', 'ASMR', 'postASMR', 'EyesOpen\_Start', 'EyesOpen\_End'}
    
    PPTStateA = fields(ft_tf_EEG.(Group1).(statelist{BaselineState}))
    PPTStateB = fields(ft_tf_EEG.(Group1).(statelist{ExperimentalState}))

    %Merge and order lists.

    PPTState_All_merged = [PPTStateA; PPTStateB]
    PPTState_All_merged = unique(PPTState_All_merged, 'rows');
    numberall = numel(PPTState_All_merged)


    % Identify which list has the most amount of common participants
    StateA_IDX = ismember(PPTState_All_merged, PPTStateA)
    StateB_IDX = ismember(PPTState_All_merged, PPTStateB)

    if sum(StateA_IDX) > sum(StateB_IDX)
        big = BaselineState
        small = ExperimentalState
    elseif sum(StateB_IDX) > sum(StateA_IDX)
        big = ExperimentalState
        small = BaselineState
    elseif sum(StateB_IDX) == sum(StateA_IDX)
        big = BaselineState
        small = ExperimentalState
    end
    %% generate a new list with a shared pool of participants
    for i = 1:length(fields(ft_tf_EEG.(Group1).(statelist{big})))
        statepptall = fields(ft_tf_EEG.(Group1).(statelist{big})) % generate particpant list from preASMR_Relax (which has no missing data)
        stateppt = statepptall{i} 
        if isfield(ft_tf_EEG.(Group1).(statelist{small}),(stateppt)) == 1 % checks for ppt in StrongASMR
            ba_vs_exp_tf.(Group1).(statelist{small}).(stateppt) = ft_tf_EEG.(Group1).(statelist{small}).(stateppt) % saves elsewhere if so
            % the same process repeats for state1pre (preASMR_relax)
            ba_vs_exp_tf.(Group1).(statelist{big}).(stateppt) = ft_tf_EEG.(Group1).(statelist{big}).(stateppt)       
        end
    end
% 

for i = 1:length(fields(ba_vs_exp_tf.(Group1))) % ba_vs_exp_tf contains data for all states of interest
        stateall = fields(ba_vs_exp_tf.(Group1)) % fetch the fieldnames of of these states
        state = stateall{i}    % select one state at a time
        stateppt = fields(ba_vs_exp_tf.(Group1).(state)) % generate participant list for this selected state
        for m = 1:length(stateppt) % loop through participants for this statelist
            ppt = stateppt{m}
            % frequency range used to normalize the spectrum
            if length(ba_vs_exp_tf.(Group1).(state).(ppt).freq) == 159
                normalisation_freq = [1 159];
            elseif length(ba_vs_exp_tf.(Group1).(state).(ppt).freq) == 80
                normalisation_freq = [1 80];
            end
            % Normalising preASMR
            ba_vs_exp_tf.(Group1).(state).(ppt).baseline = sum(ba_vs_exp_tf.(Group1).(state).(ppt).powspctrm(:,normalisation_freq(1):normalisation_freq(2)),2);
%             ba_vs_exp_tf.(state).(ppt).baseline = mean(ba_vs_exp_tf.(state).(ppt).powspctrm(:,normalisation_freq(1):normalisation_freq(2)),2);
            ba_vs_exp_tf.(Group1).(state).(ppt).powspctrm_b = bsxfun(@rdivide,ba_vs_exp_tf.(Group1).(state).(ppt).powspctrm, ba_vs_exp_tf.(Group1).(state).(ppt).baseline);
            % Normalising ASMR
        end
    end    

%%  Grand average all participants within each state (freq x freq)
%     ba_vs_exp_tf_ga = []
    cfg = []
    cfg.keepindividual = 'yes'
    cfg.parameter      = {'powspctrm_b', 'powspctrm'} % the name of the normalised power data
    statelist1 = fieldnamesr(ba_vs_exp_tf.(Group1), 1)         % generates a list of structure files with master structure (state structure files here)

    for i = 1:length(statelist1)
        structstate = statelist1{i}
        tempcell = struct2cell(ba_vs_exp_tf.(Group1).(structstate)) % converts structure into temporary cell array which is required for grandaverage function
        ba_vs_exp_tf_ga.(Group1).(structstate) = ft_freqgrandaverage(cfg, tempcell{1:end})
    end

%%  Clustering (freq x freq)
       
%         stat = []
    for j = 1:5

        freqband = [1 4 8 13 30 80]
        freqbandlabel = {'delta' 'theta' 'alpha' 'beta' 'gamma'}
        FOIlower = freqband(j)
        FOIhigher = freqband(j+1)

        cfg = [];
        cfg.latency          = 'all'; % time interval over which the experimental
                                         % conditions must be compared (in seconds)
        cfg.frequency        = [FOIlower FOIhigher];
        cfg.parameter = 'powspctrm_b'
        cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
        cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
                                       % evaluate the effect at the sample level
        cfg.correctm         = 'cluster'; % uses cluster method eg compared to Holms/Bonferonni etc
        cfg.clusteralpha     = 0.05; % alpha level of the sample-specific test statistic that
                                        % will be used for thresholding
        cfg.computeprob = 'yes'
        cfg.computecritval = 'yes'

        cfg.clusterstatistic = 'maxsum' %'maxsum'; % test statistic that will be evaluated under the
        %                                % permutation distribution.  how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
        %                           the option 'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity; 
        %                           see Hayasaka & Nichols (2004) NeuroImage for details
        cfg.clusterthreshold = 'nonparametric_common' %method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
        cfg.minnbchan        = 2; % minimum number of neighborhood channels that is
                                       % required for a selected sample to be included
                                       % in the clustering algorithm (default=0).
        cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
        cfg.clustertail      = cfg.tail;
        cfg.correcttail = 'alpha'
        cfg.alpha            = 0.05;  % alpha level of the permutation test
        cfg.numrandomization = 5000; %  In general you should set this to a higher value (1000 or up).
        % If it turns out that estimated p-value is very close to the the critical alpha-level (0.05 or 0.01), you should increase this number.

        cfg_neighb = [];
        cfg_neighb.method    = 'distance' %  'distance'; 'triangulation'%
        cfg.neighb.layout = 'biosemi64.lay'
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, ft_tf_EEG.(Group1).(statelist{BaselineState}).(ppt)); % the neighbours specify for each sensor with
                                         % which other sensors it can form clusters

        subj = size( ba_vs_exp_tf_ga.(Group1).(statelist{ExperimentalState}).powspctrm_b, 1);

        design = zeros(2,2*subj);
        for i = 1:subj
          design(1,i) = i;
        end
        for i = 1:subj
          design(1,subj+i) = i;
        end

        design(2,1:subj)        = 1;
        design(2,subj+1:2*subj) = 2;

        cfg.design           = design;  % design matrix
        cfg.ivar             = 2;  % number or list with indices indicating the independent variable(s)
                                    % We use cfg.ivar to indicate the row of the design matrix that contains the independent variable. 
                                    % For a between-trials statistical analysis, the minimum design matrix contains only a single row,
                                    % and cfg.ivar is superfluous. However, in other statistical analyses (e.g., those for a within-UO design)
                                    % the design matrix must contain more than one row, and specifying cfg.ivar is essential.
        cfg.uvar     = 1;

        [stat.(Group1).(freqbandlabel{j})] = ft_freqstatistics(cfg, ba_vs_exp_tf_ga.(Group1).(statelist{ExperimentalState}), ba_vs_exp_tf_ga.(Group1).(statelist{BaselineState})); % experimental goes first here
    end
    

 
%     stat = []
%     for j = 1
% 
%         freqband = [13 18]
%         freqbandlabel = {'lowbeta'}
%         FOIlower = freqband(1)
%         FOIhigher = freqband(1+1)
% 
%         cfg = [];
%         cfg.latency          = 'all'; % time interval over which the experimental
%                                          % conditions must be compared (in seconds)
%         cfg.frequency        = [FOIlower FOIhigher];
%         cfg.parameter = 'powspctrm_b'
%         cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
%         cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
%                                        % evaluate the effect at the sample level
%         cfg.correctm         = 'cluster'; % uses cluster method eg compared to Holms/Bonferonni etc
%         cfg.clusteralpha     = 0.05; % alpha level of the sample-specific test statistic that
%                                         % will be used for thresholding
%         cfg.computeprob = 'yes'
%         cfg.computecritval = 'yes'
% 
%         cfg.clusterstatistic = 'maxsum' %'maxsum'; % test statistic that will be evaluated under the
%         %                                % permutation distribution.  how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%         %                           the option 'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity; 
%         %                           see Hayasaka & Nichols (2004) NeuroImage for details
%         cfg.clusterthreshold = 'nonparametric_common' %method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%         cfg.minnbchan        = 2; % minimum number of neighborhood channels that is
%                                        % required for a selected sample to be included
%                                        % in the clustering algorithm (default=0).
%         cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
%         cfg.clustertail      = cfg.tail;
%         cfg.correcttail = 'alpha'
%         cfg.alpha            = 0.05;  % alpha level of the permutation test
%         cfg.numrandomization = 5000; %  In general you should set this to a higher value (1000 or up).
%         % If it turns out that estimated p-value is very close to the the critical alpha-level (0.05 or 0.01), you should increase this number.
% 
%         cfg_neighb = [];
%         cfg_neighb.method    = 'distance' %  'distance'; 'triangulation'%
%         cfg.neighb.layout = 'biosemi64.lay'
%         cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, ft_EEG.(statelist{BaselineState}).(ppt)); % the neighbours specify for each sensor with
%                                          % which other sensors it can form clusters
% 
%         subj = size( ba_vs_exp_tf_ga.(statelist{ExperimentalState}).powspctrm_b, 1);
% 
%         design = zeros(2,2*subj);
%         for i = 1:subj
%           design(1,i) = i;
%         end
%         for i = 1:subj
%           design(1,subj+i) = i;
%         end
% 
%         design(2,1:subj)        = 1;
%         design(2,subj+1:2*subj) = 2;
% 
%         cfg.design           = design;  % design matrix
%         cfg.ivar             = 2;  % number or list with indices indicating the independent variable(s)
%                                     % We use cfg.ivar to indicate the row of the design matrix that contains the independent variable. 
%                                     % For a between-trials statistical analysis, the minimum design matrix contains only a single row,
%                                     % and cfg.ivar is superfluous. However, in other statistical analyses (e.g., those for a within-UO design)
%                                     % the design matrix must contain more than one row, and specifying cfg.ivar is essential.
%         cfg.uvar     = 1;
% 
%         [stat.(freqbandlabel{j})] = ft_freqstatistics(cfg, ba_vs_exp_tf_ga.(statelist{ExperimentalState}), ba_vs_exp_tf_ga.(statelist{BaselineState})); % experimental goes first here
%     end
%         
    %% generate summary table 
    ClusterSummary.(Group1) = {'Stat' 'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma'}
    ClusterSummary.(Group1)(2,1) = {'Pos Prob'} 
    ClusterSummary.(Group1)(3,1) = {'Neg Prob'}
    try 
    ClusterSummary.(Group1)(2,2) = {stat.(Group1).delta.posclusters(1).prob(1)}
    end
    try
    ClusterSummary.(Group1)(3,2) = {stat.(Group1).delta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummary.(Group1)(2,3) = {stat.(Group1).theta.posclusters(1).prob(1)}
    end
    try
    ClusterSummary.(Group1)(3,3) = {stat.(Group1).theta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummary.(Group1)(2,4) = {stat.(Group1).alpha.posclusters(1).prob(1)}
    end
    try
    ClusterSummary.(Group1)(3,4) = {stat.(Group1).alpha.negclusters(1).prob(1)}
    end
    try 
    ClusterSummary.(Group1)(2,5) = {stat.(Group1).beta.posclusters(1).prob(1)}
    end
    try
    ClusterSummary.(Group1)(3,5) = {stat.(Group1).beta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummary.(Group1)(2,6) = {stat.(Group1).gamma.posclusters(1).prob(1)}
    end
    try
    ClusterSummary.(Group1)(3,6) = {stat.(Group1).gamma.negclusters(1).prob(1)}
    end

    
    
%% Plot raw data for baseline and experimental accross classical bands
baselinecondition = ba_vs_exp_tf_ga.(Group1).(statelist{BaselineState})
experimentalcondition = ba_vs_exp_tf_ga.(Group1).(statelist{ExperimentalState})
cfg = [];
cfg.avgoverrpt = 'yes';
cfg.parameter  = {'powspctrm_b'};

base_avg.total = ft_selectdata(cfg, baselinecondition);
exp_avg.total = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [1 4]
cfg.avgoverfreq = 'no'
base_avg.delta = ft_selectdata(cfg, baselinecondition);
exp_avg.delta = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [4 8]
base_avg.theta = ft_selectdata(cfg, baselinecondition);
exp_avg.theta = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [8 13]
base_avg.alpha = ft_selectdata(cfg, baselinecondition);
exp_avg.alpha = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [13 18]
base_avg.lowbeta = ft_selectdata(cfg, baselinecondition);
exp_avg.lowbeta = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [21 30]
base_avg.highbeta = ft_selectdata(cfg, baselinecondition);
exp_avg.highbeta = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [13 30]
base_avg.beta = ft_selectdata(cfg, baselinecondition);
exp_avg.beta = ft_selectdata(cfg, experimentalcondition);

cfg.frequency = [30 80]
base_avg.gamma = ft_selectdata(cfg, baselinecondition);
exp_avg.gamma = ft_selectdata(cfg, experimentalcondition);

%% freq x freq (difference)


% stat.delta.raweffect = ((exp_avg.delta.powspctrm_b - base_avg.delta.powspctrm_b));
% stat.beta.raweffect = ((exp_avg.beta.powspctrm_b - base_avg.beta.powspctrm_b));
% stat.theta.raweffect = ((exp_avg.theta.powspctrm_b - base_avg.theta.powspctrm_b));
% stat.alpha.raweffect = ((exp_avg.alpha.powspctrm_b - base_avg.alpha.powspctrm_b));
% stat.gamma.raweffect = ((exp_avg.gamma.powspctrm_b - base_avg.gamma.powspctrm_b));

%% freq x freq (% difference)

stat.(Group1).delta.raweffect = ((exp_avg.delta.powspctrm_b - base_avg.delta.powspctrm_b)./ abs(base_avg.delta.powspctrm_b))*100;
stat.(Group1).theta.raweffect = ((exp_avg.theta.powspctrm_b - base_avg.theta.powspctrm_b)./ abs(base_avg.theta.powspctrm_b))*100;
stat.(Group1).alpha.raweffect = ((exp_avg.alpha.powspctrm_b - base_avg.alpha.powspctrm_b)./ abs(base_avg.alpha.powspctrm_b))*100;
stat.(Group1).lowbeta.raweffect = ((exp_avg.lowbeta.powspctrm_b - base_avg.lowbeta.powspctrm_b)./ abs(base_avg.lowbeta.powspctrm_b))*100;
stat.(Group1).highbeta.raweffect = ((exp_avg.highbeta.powspctrm_b - base_avg.highbeta.powspctrm_b)./ abs(base_avg.highbeta.powspctrm_b))*100;
stat.(Group1).beta.raweffect = ((exp_avg.beta.powspctrm_b - base_avg.beta.powspctrm_b)./ abs(base_avg.beta.powspctrm_b))*100;
stat.(Group1).gamma.raweffect = ((exp_avg.gamma.powspctrm_b - base_avg.gamma.powspctrm_b)./ abs(base_avg.gamma.powspctrm_b))*100;
ClusterSummary.(Group1)

ClusterSig = ClusterSummary.(Group1)(2:3,2:end)
empties = cellfun('isempty',ClusterSig)
ClusterSig(empties) = {NaN}
ClusterSig = cell2mat(ClusterSig)

statsalpha = 0.025 % 0.025 for 1-tail etc

ClusterSigIDX = ClusterSig < statsalpha
ClusterSigIDX_any = (any ( ClusterSigIDX, 1 ))

%% plot the this difference with the pos/neg cluster mask overlayed.


freqband = [1 4 8 13 30 80]
freqbandlabel = {'delta' 'theta' 'alpha' 'beta' 'gamma'}


%   1×5 cell array
%
%     {'delta'}    {'theta'}    {'alpha'}    {'beta'}    {'gamma'}

for i = 1:length(freqbandlabel)
    if ClusterSigIDX_any(i) == 1
        bandindex = i
        bandlabel = freqbandlabel{bandindex}
        % bandlabel = 'lowbeta'
        freqbandplot = []
        if bandlabel == "delta"
            freqbandplot = [1 4]
        elseif bandlabel == "theta"
            freqbandplot = [4 8]
        elseif bandlabel == "alpha"
            freqbandplot = [8 13]
        elseif bandlabel == "beta"
            freqbandplot = [13 30]
        elseif bandlabel == "gamma"
            freqbandplot = [30 80]
        end
        
        
        cfg = [];
        cfg.alpha  = statsalpha
        cfg.parameter = 'raweffect'; % be sure that the dimensions of this matches the dimensions of stat. eg number of frequencies involved
        cfg.zlim   = 'maxabs'%[0 2]; % set this only after running without initially to find out appropriate range to standardise all maps to.
        cfg.layout = 'biosemi64.lay';
        cfg.gridscale = 200
        cfg.subplotsize = [3 3]
        cfg.colorbartext       =  'difference'
        cfg.colorbar           = 'EastOutside'
        % cfg.marker = 'on'
        % cfg.markersize = 10
        
        % ft_clusterplot(cfg, stat.(Group1).delta);
        % ft_clusterplot(cfg, stat.(Group1).theta);
%         ft_clusterplot(cfg, stat.(Group1).(bandlabel));
        % ft_clusterplot(cfg, stat.(Group1).gamma);
        % ft_clusterplot(cfg, stat.(Group1).beta);
        
        %% lets squeeze down the freqs and sig electrodes across the band
        statbackup = stat
        % stat = statbackup
        try
            if  stat.(Group1).(bandlabel).posclusters(1).prob < statsalpha
                if size(stat.(Group1).(bandlabel).posclusters, 2) == 1
                    stat.(Group1).(bandlabel).siglabelmat1 = stat.(Group1).(bandlabel).mask
                    stat.(Group1).(bandlabel).SigHzIDX = unique(ceil((find(stat.(Group1).(bandlabel).posclusterslabelmat==1) / 64)))'
                elseif size(stat.(Group1).(bandlabel).posclusters, 2)  > 1
                    stat.(Group1).(bandlabel).siglabelmat1 = zeros(size(stat.(Group1).(bandlabel).posclusterslabelmat))
                    stat.(Group1).(bandlabel).siglabelmat1(find(stat.(Group1).(bandlabel).posclusterslabelmat==1)) = 1
                    stat.(Group1).(bandlabel).posclusters = stat.(Group1).(bandlabel).posclusters(1)
                    stat.(Group1).(bandlabel).SigHzIDX = unique(ceil((find(stat.(Group1).(bandlabel).posclusterslabelmat==1) / 64)))'
                    %     elseif stat.(Group1).(bandlabel).posclusters(1).prob > 0.025
                    %         stat.(Group1).(bandlabel) = rmfield(stat.(bandlabel), {'posclusterslabelmat' 'posdistribution' 'posclusters'})
                end
            end
        end
        try
            if  stat.(Group1).(bandlabel).negclusters(1).prob < statsalpha
                if size(stat.(Group1).(bandlabel).negclusters, 2) == 1
                    stat.(Group1).(bandlabel).siglabelmat1 = stat.(Group1).(bandlabel).mask
                    stat.(Group1).(bandlabel).SigHzIDX = unique(ceil((find(stat.(Group1).(bandlabel).negclusterslabelmat==1) / 64)))'
                end
                if size(stat.(Group1).(bandlabel).negclusters, 2) > 1
                    stat.(Group1).(bandlabel).siglabelmat1 = zeros(size(stat.(Group1).(bandlabel).negclusterslabelmat))
                    stat.(Group1).(bandlabel).siglabelmat1(find(stat.(Group1).(bandlabel).negclusterslabelmat==1)) = 1
                    stat.(Group1).(bandlabel).negclusters = stat.(Group1).(bandlabel).negclusters(1)
                    stat.(Group1).(bandlabel).SigHzIDX = unique(ceil((find(stat.(Group1).(bandlabel).negclusterslabelmat==1) / 64)))'
                end
                %     elseif stat.(Group1).(bandlabel).negclusters(1).prob > 0.025
                %         stat.(Group1).(bandlabel) = rmfield(stat.(bandlabel), {'negclusterslabelmat' 'negdistribution' 'negclusters'})
            end
        end
        stat.(Group1).(bandlabel).SigHz = stat.(Group1).(bandlabel).freq(any(stat.(Group1).(bandlabel).siglabelmat1,1))
        %
        %     %% effect size
        %
        %        for j = bandindex
        %
        %         freqband = [1 4 8 13 30 80]
        %         freqbandlabel = {'delta' 'theta' 'alpha' 'beta' 'gamma'}
        %         FOIlower = freqband(j)
        %         FOIhigher = freqband(j+1)
        %
        %
        %         cfg = [];
        %         cfg.latency          = 'all'; % time interval over which the experimental
        %                                          % conditions must be compared (in seconds)
        % %         cfg.channel = stat.(freqbandlabel{j}).label(any(stat.(freqbandlabel{j}).mask,2))
        %         cfg.frequency        = [FOIlower FOIhigher];
        %         cfg.parameter = 'powspctrm_b'
        %         cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
        %         cfg.statistic        = 'cohensd'; % use the independent samples T-statistic as a measure to
        %                                        % evaluate the effect at the sample level
        %         cfg.correctm         = 'cluster'; % uses cluster method eg compared to Holms/Bonferonni etc
        %         cfg.clusteralpha     = 0.05; % alpha level of the sample-specific test statistic that
        %                                         % will be used for thresholding
        %         cfg.computeprob = 'yes'
        %         cfg.computecritval = 'yes'
        %
        %         cfg.clusterstatistic = 'maxsum' %'maxsum'; % test statistic that will be evaluated under the
        %         %                                % permutation distribution.  how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
        %         %                           the option 'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity;
        %         %                           see Hayasaka & Nichols (2004) NeuroImage for details
        %         cfg.clusterthreshold = 'nonparametric_common' %method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
        %         cfg.minnbchan        = 2; % minimum number of neighborhood channels that is
        %                                        % required for a selected sample to be included
        %                                        % in the clustering algorithm (default=0).
        %         cfg.tail             = -1; % -1, 1 or 0 (default = 0); one-sided or two-sided test
        %         cfg.clustertail      = cfg.tail;
        %         cfg.correcttail = 'alpha'
        %         cfg.alpha            = 0.025;  % alpha level of the permutation test
        %         cfg.numrandomization = 5000; %  In general you should set this to a higher value (1000 or up).
        %         % If it turns out that estimated p-value is very close to the the critical alpha-level (0.05 or 0.01), you should increase this number.
        %
        %         cfg_neighb = [];
        %         cfg_neighb.method    = 'distance' %  'distance'; 'triangulation'%
        %         cfg.neighb.layout = 'biosemi64.lay'
        %         cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, ft_EEG.(statelist{BaselineState}).(ppt)); % the neighbours specify for each sensor with
        %                                          % which other sensors it can form clusters
        %
        %         subj = size( ba_vs_exp_tf_ga.(statelist{ExperimentalState}).powspctrm_b, 1);
        %
        %         design = zeros(2,2*subj);
        %         for i = 1:subj
        %           design(1,i) = i;
        %         end
        %         for i = 1:subj
        %           design(1,subj+i) = i;
        %         end
        %
        %         design(2,1:subj)        = 1;
        %         design(2,subj+1:2*subj) = 2;
        %
        %         cfg.design           = design;  % design matrix
        %         cfg.ivar             = 2;  % number or list with indices indicating the independent variable(s)
        %                                     % We use cfg.ivar to indicate the row of the design matrix that contains the independent variable.
        %                                     % For a between-trials statistical analysis, the minimum design matrix contains only a single row,
        %                                     % and cfg.ivar is superfluous. However, in other statistical analyses (e.g., those for a within-UO design)
        %                                     % the design matrix must contain more than one row, and specifying cfg.ivar is essential.
        %         cfg.uvar     = 1;
        %
        %         [cohensdstat.(freqbandlabel{j})] = ft_freqstatistics(cfg, ba_vs_exp_tf_ga.(statelist{ExperimentalState}), ba_vs_exp_tf_ga.(statelist{BaselineState})); % experimental goes first here
        %        end
        % sigvoxels = []
        % sigvoxels = false(size(stat.(freqbandlabel{j}).siglabelmat1))
        % sigvoxels = any(stat.(freqbandlabel{j}).siglabelmat1,3)
        % stat.(freqbandlabel{j}).meancohensd = mean(cohensdstat.(freqbandlabel{j}).cohensd(sigvoxels))
        % sprintf('mean cohens d = %.4g', stat.(freqbandlabel{j}).meancohensd)
        
        %%
        
        statplot_mean = stat.(Group1).(bandlabel)
        statplot_mean.freq = mean ( stat.(Group1).(bandlabel).freq(:, statplot_mean.SigHzIDX) );
        statplot_mean.raweffect = mean( stat.(Group1).(bandlabel).raweffect(:, statplot_mean.SigHzIDX), 2 );
        statplot_mean.prob = mean ( stat.(Group1).(bandlabel).prob(:, statplot_mean.SigHzIDX), 2 );
        try
            statplot_mean.posclusterslabelmat = (mean(statplot_mean.siglabelmat1,2))
            statplot_mean.posclusterslabelmat = ceil(mean(statplot_mean.mask,2))
            statplot_mean.figuredescript = single(round(statplot_mean.posclusters(1).prob,3,'significant'));
            
        end
        try
            statplot_mean.negclusterslabelmat = (mean(statplot_mean.siglabelmat1,2))
            statplot_mean.negclusterslabelmat = ceil(mean(statplot_mean.mask,2))
            statplot_mean.figuredescript = single(round(statplot_mean.negclusters(1).prob,3,'significant'));
        end
        
        if isempty(statplot_mean.negclusters) == 1
            statplot_mean = rmfield(statplot_mean, {'negclusterslabelmat' 'negdistribution' 'negclusters'})
        end
        if isempty(statplot_mean.posclusters) == 1
            statplot_mean = rmfield(statplot_mean, {'posclusterslabelmat' 'posdistribution' 'posclusters'})
        end
        
        
        cfg = [];
        cfg.alpha  = statsalpha
        cfg.parameter = 'raweffect'; % be sure that the dimensions of this matches the dimensions of stat. eg number of frequencies involved
        cfg.zlim   = 'maxabs'%[0 2]; % set this only after running without initially to find out appropriate range to standardise all maps to.
        % cfg.zlim   = [-15 15]
        cfg.layout = 'biosemi64.lay';
        cfg.gridscale = 200
        cfg.subplotsize = [1  1]
        cfg.colorbartext       =  'difference'
        cfg.colorbar           = 'EastOutside'
        cfg.maskparameter      = 'mask'
        sigmarkersize1 = 7
        sigmarkersize5 = 10
        cfg.highlightsizeseries = [sigmarkersize1 sigmarkersize5 sigmarkersize5 sigmarkersize5 sigmarkersize5]
        cfg.comment = 'no'
        
        Testfig = figure;
        Testfig = ft_clusterplot(cfg, statplot_mean);
        % title(sprintf('%s vs %s in the %s band (sig %d-%d Hz; n = %d)',...
        % (statedescript{ExperimentalState}), (statedescript{BaselineState}), bandlabel, (statplot_mean.cfg.frequency(1)+(statplot_mean.SigHzIDX(1))-1) , (statplot_mean.cfg.frequency(1)+(statplot_mean.SigHzIDX(end))-1), subj), 'FontSize', 12)
        cbarHandle = colorbar;
        % cbarHandle.Label.String = sprintf('%% change from %s', (statedescript{BaselineState}))
        cbarHandle.Label.String = '% change'
        % cbarHandle.Label.FontSize = 12;
        cbarHandle.Label.FontSize = 14;
        cbarHandle.FontSize = 14
        
        ClusterSummary
        
        statelistlab = strrep(statelist,'_','\_')
        
        colormap jet
        sprintf('%s %s vs %s in the %s band (sig %g-%g Hz; n = %d; p = %g)',...
            Group1, statelistlab{ExperimentalState}, (statelistlab{BaselineState}), bandlabel, (stat.(Group1).(bandlabel).SigHz(1)) ,...
            stat.(Group1).(bandlabel).SigHz(end), subj, statplot_mean.figuredescript)
       
        % colorbar off % useful for multiplots
        % title('') % useful for multiplots
        title(sprintf('%s %s vs %s in the %s band (sig %g-%g Hz; n = %d; p = %g)',...
            Group1, statelistlab{ExperimentalState}, (statelistlab{BaselineState}), bandlabel, (stat.(Group1).(bandlabel).SigHz(1)) ,...
            stat.(Group1).(bandlabel).SigHz(end), subj, statplot_mean.figuredescript))
        
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        saveas(gcf,sprintf('C:/MATLAB/EEGData/ASMR-R/FT/%s_%s-vs-%s_%s.png',Group1, statelist{ExperimentalState}, (statelist{BaselineState}), bandlabel))

        % x0=100;
        % y0=100;
        % width=450;  % 300 is good for 5 row table
        % height=250 % 200 is good for 5 row table
        % set(gcf,'position',[x0,y0,width,height])
        
        %% plot power spectra across all electrodes, sig freqs (x axis) masked with a box (y = power)
        
        % sigposmask = (stat.alpha.posclusterslabelmat==1) & stat.alpha.mask;
        % signegmask = (stat.alpha.negclusterslabelmat==1) & stat.alpha.mask;
        %
        %
        % base_avg.alpha.mask = sigposmask;
        % exp_avg.alpha.mask = sigposmask;
        %
        % cfg = [];
        % % cfg.elec          = elec;
        % cfg.layout = 'biosemi64.lay'
        % cfg.colorbar      = 'no';
        % cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
        % % cfg.maskstyle     = 'thickness';
        % cfg.parameter     = 'powspctrm_b';
        % cfg.maskfacealpha = 0.5;
        % cfg.showlabels = 'yes'
        % cfg.showscale     = 'yes'
        % % cfg.showoutline   = 'yes'
        %
        % figure; ft_multiplotER(cfg, base_avg.alpha, exp_avg.alpha);
        % title('within-participant preASMR-Baseline vs ASMR');
        % cfg.freqlow = 1:2
        % cfg.freqhigh = 9:13
        % cfg.method = coh
        % %
        % % ft_crossfrequencyanalysis
    else
    end
end
end
%%
%%
%%
%%
%%

%% Run group comparisons for each state (e.g., EO_Start or EO_End)



for BaselineState = 1:2
    
    for j = 1:5

        freqband = [1 4 8 13 30 80]
        freqbandlabel = {'delta' 'theta' 'alpha' 'beta' 'gamma'}
        FOIlower = freqband(j)
        FOIhigher = freqband(j+1)

        cfg = [];
        cfg.latency          = 'all'; % time interval over which the experimental
                                         % conditions must be compared (in seconds)
        cfg.frequency        = [FOIlower FOIhigher];
        cfg.parameter = 'powspctrm_b'
        cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
        cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
                                       % evaluate the effect at the sample level
        cfg.correctm         = 'cluster'; % uses cluster method eg compared to Holms/Bonferonni etc
        cfg.clusteralpha     = 0.05; % alpha level of the sample-specific test statistic that
                                        % will be used for thresholding
        cfg.computeprob = 'yes'
        cfg.computecritval = 'yes'

        cfg.clusterstatistic = 'maxsum' %'maxsum'; % test statistic that will be evaluated under the
        %                                % permutation distribution.  how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
        %                           the option 'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity; 
        %                           see Hayasaka & Nichols (2004) NeuroImage for details
        cfg.clusterthreshold = 'nonparametric_common' %method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
        cfg.minnbchan        = 2; % minimum number of neighborhood channels that is
                                       % required for a selected sample to be included
                                       % in the clustering algorithm (default=0).
        cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
        cfg.clustertail      = cfg.tail;
        cfg.correcttail = 'alpha'
        cfg.alpha            = 0.05;  % alpha level of the permutation test
        cfg.numrandomization = 5000; %  In general you should set this to a higher value (1000 or up).
        % If it turns out that estimated p-value is very close to the the critical alpha-level (0.05 or 0.01), you should increase this number.

        cfg_neighb = [];
        cfg_neighb.method    = 'distance' %  'distance'; 'triangulation'%
        cfg.neighb.layout = 'biosemi64.lay'
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, ft_tf_EEG.(Group1).(statelist{BaselineState}).(ppt)); % the neighbours specify for each sensor with
                                         % which other sensors it can form clusters

        design = zeros(1,size(ba_vs_exp_tf_ga.(grouplist{1}).(statelist{BaselineState}).powspctrm_b,1) + size(ba_vs_exp_tf_ga.(grouplist{2}).(statelist{BaselineState}).powspctrm_b,1));
        design(1,1:size(ba_vs_exp_tf_ga.(grouplist{1}).(statelist{BaselineState}).powspctrm_b,1)) = 1;
        design(1,(size(ba_vs_exp_tf_ga.(grouplist{1}).(statelist{BaselineState}).powspctrm_b,1)+1):(size(ba_vs_exp_tf_ga.(grouplist{1}).(statelist{BaselineState}).powspctrm_b,1)+...
            size(ba_vs_exp_tf_ga.(grouplist{2}).(statelist{BaselineState}).powspctrm_b,1))) = 2;
        cfg.design           = design;
        cfg.ivar             = 1;  % number or list with indices indicating the independent variable(s)
        % We use cfg.ivar to indicate the row of the design matrix that contains the independent variable.
        % For a between-trials statistical analysis, the minimum design matrix contains only a single row,
        % and cfg.ivar is superfluous. However, in other statistical analyses (e.g., those for a within-UO design)
        % the design matrix must contain more than one row, and specifying cfg.ivar is essential.
        
        [statbetween.(statelist{BaselineState}).(freqbandlabel{j})] = ft_freqstatistics(cfg, ba_vs_exp_tf_ga.(grouplist{1}).(statelist{BaselineState}), ba_vs_exp_tf_ga.(grouplist{2}).(statelist{BaselineState})); % experimental goes first here
    end
    
 ClusterSummaryBetween.(statelist{BaselineState}) = {'Stat' 'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma'}
    ClusterSummaryBetween.(statelist{BaselineState})(2,1) = {'Pos Prob'} 
    ClusterSummaryBetween.(statelist{BaselineState})(3,1) = {'Neg Prob'}
    try 
    ClusterSummaryBetween.(statelist{BaselineState})(2,2) = {statbetween.(statelist{BaselineState}).delta.posclusters(1).prob(1)}
    end
    try
    ClusterSummaryBetween.(statelist{BaselineState})(3,2) = {statbetween.(statelist{BaselineState}).delta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummaryBetween.(statelist{BaselineState})(2,3) = {statbetween.(statelist{BaselineState}).theta.posclusters(1).prob(1)}
    end
    try
    ClusterSummaryBetween.(statelist{BaselineState})(3,3) = {statbetween.(statelist{BaselineState}).theta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummaryBetween.(statelist{BaselineState})(2,4) = {statbetween.(statelist{BaselineState}).alpha.posclusters(1).prob(1)}
    end
    try
    ClusterSummaryBetween.(statelist{BaselineState})(3,4) = {statbetween.(statelist{BaselineState}).alpha.negclusters(1).prob(1)}
    end
    try 
    ClusterSummaryBetween.(statelist{BaselineState})(2,5) = {statbetween.(statelist{BaselineState}).beta.posclusters(1).prob(1)}
    end
    try
    ClusterSummaryBetween.(statelist{BaselineState})(3,5) = {statbetween.(statelist{BaselineState}).beta.negclusters(1).prob(1)}
    end
    try 
    ClusterSummaryBetween.(statelist{BaselineState})(2,6) = {statbetween.(statelist{BaselineState}).gamma.posclusters(1).prob(1)}
    end
    try
    ClusterSummaryBetween.(statelist{BaselineState})(3,6) = {statbetween.(statelist{BaselineState}).gamma.negclusters(1).prob(1)}
    end
    
    ClusterSummaryBetween.(statelist{BaselineState})
end


