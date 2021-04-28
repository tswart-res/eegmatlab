%% EEGLAB preprocessing 


clear all; close all;
%% Paths

state1 = 'Controls' 
% state1 = 'ASMR-R' 
% condition1 = 'EO_start'
condition1 = 'EO_end'


cd C:\MATLAB;
pardir = pwd;

% EEG data path
EEGdatadir = sprintf('C:\\MATLAB\\EEGData\\%s\\', state1);
addpath(genpath(EEGdatadir));

% Scripts
scriptspath = [pardir + "\Thomas_Scripts\"];

addpath(genpath(scriptspath));

cd (EEGdatadir);

% EEGLAB path
% eeglabpath = ('C:\\Program Files\\MATLAB\\R2018b\\toolbox\\eeglab14_1_2b');
eeglabpath = ('C:\\MATLAB\\eeglab14_1_2b');
addpath(eeglabpath);
%% 
% Input these biosemi specific montage things
chEOG = [67, 68];
chECG = [69, 70];
chREF = [65, 66];

%% For Loop to import all BDF data and filter / triggger    

% Run 1

% Subjects to be analysed (From run 1)
% subj= [13, 16, 17, 18, 23, 38, 39, 40, 42, 44, 48, 51, 53, 55, 56, 58, 59, 61, 62, 64, 67];
% subj= [13, 16, 17, 18, 23, 38];
% subj= [39, 40, 42, 44, 48, 51, 53, 55, 56, 58, 59, 61, 62, 64, 67];

file_list = dir([EEGdatadir filesep() sprintf('*_%s.bdf', condition1)]); % generate a filename list based on 0/1/2/3 subfolder filenames
numsubjects = length(file_list);
file_name_list = {file_list.name};

for x=1:numsubjects                                 %  extract participant numbers from filename list
    file_name_list{x} = file_name_list{x}(3:4);
end

pnumber = unique(file_name_list, 'rows');

% Define output path



% Load EEGLAB
%% Start
for s = 1:numsubjects
% for s = 1
    eeglab;close all;
    p = pnumber{s}
    eegOUTpath = (sprintf('\\%s_%s\\', condition1, state1));
    filename = sprintf('P0%s_%s.bdf',p, condition1);
    EEG = pop_biosig(filename,'ref',chREF) %
    %EEG =  pop_chanedit(EEG1,[]);
    EEG = pop_chanedit(EEG, 'lookup','Standard-10-5-Cap385.sfp');
   % EEG=pop_chanedit(EEG, 'lookup','C:\\Program Files\\MATLAB\\R2018b\\toolbox\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_reref(EEG, [65, 66])
    % Highpass filter at 0.5Hz (filter order 3380)
     EEG = pop_eegfiltnew(EEG,0.5, [], 3380); 
    % Notch filter to remove powerline noise
     EEG = pop_eegfiltnew(EEG,48,52,[],1);    
    % EEG triggers converted to original intended number: e.g. 50, 93, 94...
   %  EEG = pop_editeventvals(EEG,'append',{5 [] [] []},'changefield',{6 'type' 3490},'changefield',{6 'type' 63490},'changefield',{6 'latency' 726});

    [EventType , Time, EEG] =   convert_trg(EEG) 
    % copy each event every 2 seconds till different event
    try
    EEG = EventReplicator(EEG) 
    end
    % Cut out 2 mins of fixation cross at 10-12 mins and 20-22 mins (+ 2
    % seconds for initial delay)
    EventLatency0 = cell2mat({EEG.event.latency})
    EventType0 = cell2mat({EEG.event.type})
    VideoEnd105IDX = find(EventType0 == 105)
    try
    VideoEnd107IDX = find(EventType0 == 107)
    end
    if exist('VideoEnd107IDX', 'var') == 1
        if length(VideoEnd105IDX) == 1
            if VideoEnd107IDX > VideoEnd105IDX
                VideoEnd105IDX(2) = VideoEnd107IDX
            elseif VideoEnd105IDX > VideoEnd107IDX
                VideoEnd105IDX(2) = VideoEnd105IDX
                VideoEnd105IDX(1) = VideoEnd107IDX
            end
        elseif length(VideoEnd105IDX) == 0
            VideoEnd105IDX = VideoEnd107IDX
        end
    end
    Video2StartIDX = find(EventType0 == 102)
    Video1StartIDX = find(EventType0 == 101)
    if Video2StartIDX > Video1StartIDX
        MiddleVideoStartIDX = Video2StartIDX
    elseif Video1StartIDX > Video2StartIDX
        MiddleVideoStartIDX = Video1StartIDX
    end
    try
    LatencyFirstEnd = EventLatency0(VideoEnd105IDX(1))/512
    end
    try
    LatencyFinalEnd = EventLatency0(VideoEnd105IDX(2))/512
    end
    try
    LatencyMiddleVideoStart = EventLatency0(MiddleVideoStartIDX)/512
    end
    EndofRecording = (length(EEG.times)/512)
    
   
    EEG.data(65,:) = diff(EEG.data([66, 65],:)); % calculate difference waveforms for HEOG
    EEG.data(66,:) = diff(EEG.data([68, 67],:)); % calculate difference waveforms for VEOG
    EEG.data(67,:) = diff(EEG.data([69, 70],:)); % calculate difference waveforms for HR
    
%     EEG = pop_select( EEG, 'nochannel',[65, 68:70] ); %%%% TEST
    EEG = pop_select( EEG, 'nochannel',[68:70] ); % remove reference channels now that difference waveforms have been calculated

%     TimePoints = length(EEG.times)./512;
% 
%     for m = 1:2:TimePoints
%         EEG = pop_editeventvals(EEG,'insert',{1 [] [] []},'changefield',{1 'type' 0},'changefield',{1 'latency' m});
%     end
%     EEG = eeg_checkset( EEG );

        %% EventLengthDelineator!
        try
        allEventTypes = {EEG.event}';
        allEventTypesMat = cell2mat(allEventTypes);
        N = length(allEventTypesMat); % number of events
        if N > 1
            eventsummary = cell2mat(permute(struct2cell(EEG.event), [3,1,2]));
            whilecounter = 0;
            eventsummary(:,4) = 0;
            eventsummary(:,2) = eventsummary(:,2)./512;
%             indices1 = find(eventsummary(:,1) == 240)
%             eventsummary(indices1,:) = []
            indices2 = find(eventsummary(:,1) == 10);
            eventsummary(indices2,:) = [];
            indices3 = find(eventsummary(:,1) == 13);
            eventsummary(indices3,:) = [];
            NewN = length(eventsummary);
            o=0;
            for o = 2:NewN
                if o+1 > NewN
                    eventsummary(o+1,:) = [4];
                end
                if eventsummary(o,1) ~= eventsummary(o-1,1)
                    if whilecounter > 0
                        eventsummary(o-1-whilecounter,4) = (eventsummary(o,2) - eventsummary(o-1-whilecounter,2));
                        whilecounter = 0;
                    else
                        eventsummary(o-1,4) = (eventsummary(o,2) - eventsummary(o-1,2))
                    end
                elseif eventsummary(o,1) == eventsummary(o-1,1)
                    whilecounter = whilecounter + 1;
                    continue
                end
            end
            eventsummary(NewN+1,:) = [];
            %% remove zero rows
            indices = find(eventsummary(:,4) == 0);
            eventsummary(indices,:) = [];
            
            %% Take a .5 second after/before
            eventsummary(:,5) = (eventsummary(:,2)+0.5);
            eventsummary(:,6) = (eventsummary(:,5)+eventsummary(:,4)-1);
            %% Updated duration scores
            eventsummary(:,7) = (eventsummary(:,6)-eventsummary(:,5));
            %% Put back into EEG struct
            EEG.DurationSummary = eventsummary;
        end
        end
        
        try
        [EEG.event2] = ERPLABTriggerConverter(EEG.event)
        end
        try
        EEG = pop_select( EEG,'notime',[LatencyFinalEnd EndofRecording] )
        end
        try
        EEG = pop_select( EEG,'notime',[LatencyFirstEnd (LatencyMiddleVideoStart - 0.5)] )
        end
        try
        EEG.event2((VideoEnd105IDX(2)+1):end) = []
        end
        try
        EEG.event2((VideoEnd105IDX(1)+1):(MiddleVideoStartIDX-1)) = []
        end

    % copy each event every 2 seconds till different event
    %EEG = EventReplicator(EEG) 
    % Cut out 2 mins of fixation cross at 10-12 mins and 20-22 mins (+ 2
    % seconds for initial delay)
    % Save
    % Creates folder for SPM output

    % Names
    eegoutlabel = fullfile(EEGdatadir, condition1);
    eegoutlabel = [eegoutlabel '_PSD']
    fileOUTname = filename(1:end-4);
    fileOUTname = [fileOUTname,'_filt'];
    % Go to place to save 
    if ~exist ([eegoutlabel],'dir')
        mkdir ([eegoutlabel]);
    end
    cd(eegoutlabel);
    % Save
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',fileOUTname,'saveold',filename(1:end-4),'gui','off'); 
%     eeglab redraw 
%%%% ASR remove bad channels
    originalEEG = EEG;
    EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','on','Distance','Euclidian','availableRAM_GB', 16);

    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',sprintf('%s_badchantemp',fileOUTname),'savenew',sprintf('%s_badchantemp.set',fileOUTname),'gui','off');
%     eeglab redraw
    %%%% Save channel mask outside of new EEG file created to be importable
    EEGCleanchannelmask = []
    try EEGCleanchannelmask = EEG.etc.clean_channel_mask
    end
    %%%% Ensure externals are preserved
    try EEGCleanchannelmask(65:end) = 1
    end
    if nnz(EEGCleanchannelmask > 0) 
        %%%% load previous filtered file with all channels
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',2,'study',0);
%         eeglab redraw 

        %%%% Omit noisey channels based on ASR
        EEG = pop_select(EEG, 'channel', find(EEGCleanchannelmask))
        EEG.etc.clean_channel_mask = EEGCleanchannelmask 

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',sprintf('%s_badchan_R',fileOUTname),'savenew',sprintf('%s_badchan_R.set',fileOUTname),'gui','off');
%         eeglab redraw 

        %%%% Interp from raw, filtered dataset (indicated by the ALLEEG(1) -
        %%%% indicating the indices of dataset used for interpolation)
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

%         EEG = pop_interp(EEG, ALLEEG(2).chanlocs(find(~EEGCleanchannelmask)), 'spherical');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',sprintf('%s_badchan_R_Interp',fileOUTname),'savenew',sprintf('%s_badchan_R_Interp.set',fileOUTname),'gui','off');
%         eeglab redraw
    else
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',2,'study',0);
%         eeglab redraw 
    end
    %%%%%% Epoch before ICA (remember ICA doesnt care if contigous -
    %%%%%% assumes data is all related
%     EEG = pop_epoch( EEG, {  '10'  '13'  }, [-0.1006         0.5], 'newname', sprintf('P0%s_Run2_epoch_trainingICA', p), 'epochinfo', 'yes');
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',sprintf('P0%s_Run2_epoch_trainingICA', p),'gui','off'); 
%     EEG = eeg_checkset( EEG );
%     EEG = pop_rmbase( EEG, [-100 0] ,[]);
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',sprintf('P0%s_Run2_epoch_trainingICA', p),'gui','off'); 

    fileOUTnameICA = fileOUTname(1:end-5);
    fileOUTnameICA = [fileOUTnameICA,'_ICA']
    
  %  ------------------------------------------------
%     %[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%     %EEG = pop_loadset('filename',fileOUTname,'filepath','C:\\MATLAB\\EEGData\\Run1_EEGLAB\\');
%     %[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
if isempty(EEGCleanchannelmask) == 1
    EEG = eeg_checkset( EEG );
    EEG = pop_runica(EEG, 'extended',1,'interupt','on');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
elseif isempty(EEGCleanchannelmask) == 0
    EEG = eeg_checkset( EEG );
    EEG = pop_runica(EEG, 'extended',1,'interupt','on', 'pca', nnz(EEGCleanchannelmask));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end
    EEG = eeg_checkset( EEG );
%     eeglab redraw
    EEG = eeg_checkset( EEG );
    EEG = pop_editset(EEG, 'setname', sprintf('P0%s_%s_ICA.set', p, condition1));
%     EEG = pop_editset(EEG, 'setname', sprintf('P0%s_%s_ICA.set', p, condition1), 'run', []);

    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',sprintf('P0%s_%s_ICA.set', p, condition1),'filepath', sprintf('C:\\MATLAB\\EEGData\\%s\\%s\\',state1, [condition1 '_PSD']));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%     eeglab redraw


end
%%%%%############%%%%%%%%%%%%%        AFTER ICA HERE



