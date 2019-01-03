function [ classificated_stimuli, nRounds, nStimuli, error ] = classificate_fixations(datapath, filename)
%CLASIFICATE_FIXATIONS Summary of this function goes here
%   
%   classificated_stimuli struct: One element per stimulus
%        - NHC: Patient ID
%        - Algorithm: Name of the algorithm used
%        - Round: Number of the test round
%        - Stimulus: ID of the stimulus
%        - StimPos: screen position of the stimulus (mm) and the coord 
%               start at bottom left corner.
%        - Timestamp: time in which the stamp was register
%        - RawVelocities: velocities at each timestamp (mm/s)
%        - Velocities: velocities for each point mm/s
%        - RawGazePos: gaze position on the screen of the stamp (mm)
%        - GazePos: gaze position on the screen of the stamp (mm)
%        - Classification: 1 if the stamp was classificated as Fixation, 
%              2 if saccade and 3 if glissade. 0 as non defined
% 
%
error = 0;

fDiscardOutliers = true;
f_rawData = true;
f_EyeSel = 0;

[test, avatar_test,f_emptyFile] = ff_load_and_setup(filename, datapath,...
    fDiscardOutliers, f_rawData, f_EyeSel);
nRounds = test.NROUNDS;
nStimuli = test.NSTIMULI;

if(f_emptyFile)
    error = 1;
    return;
end

%% BIT
ID = zeros(test.NPoints,1);
ID(:) = test.NHC;
distance = zeros(test.NPoints,1);
distance(:) = 55;
rawVelocity = computeLinVelocity(test.RawLeftPos(:,1), test.RawLeftPos(:,2), test.T, test.DPI);
timestamp = test.T*1000; % seconds to ms
stamps = [ ID, timestamp, test.RawLeftPos(:,1), test.RawRightPos(:,1), test.RawLeftPos(:,2), test.RawRightPos(:,2), distance];
[output_id,output_total] = BIT(stamps);

%%Divede the results into stimuli (ouptut skip the first stamp?)
%borrar NAN en validity dividir con find(onSets)
stimuliIdx = find(test.onsets);
filtDer = [1; -1]; % x - x_1
% Add the idx of the end of the last stimulus.
grad = conv(filtDer, test.sampleValidity);
aux = find(grad);
if(aux(end)>size(test.sampleValidity,1))
    aux(end) = size(test.sampleValidity,1);
end
stimuliIdx(size(stimuliIdx, 1)+1) = aux(end);

% stamps_position save the id of the fixation which the stamps belong
% 0 means is not a fixation, so classifications will mark as 1 the stamps
% that belong to a fixation.
classifications = output_total.stamps_positions(:,6) ~= 0;
idx_CS = 1;
idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'BIT';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
        classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI);
        classificated_stimuli(idx_CS).GazePos = [pix2mm(output_total.stamps_fixations(infIdx:supIdx,3), test.DPI), pix2mm(output_total.stamps_fixations(infIdx:supIdx,5), test.DPI)]; %left eye
        classificated_stimuli(idx_CS).Classification = classifications(infIdx:supIdx) %If not 0 is a fixation
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end    
end


%% NH 
%
% Execute NH with all gaze data
deleteInputNHFiles()
%Save the input file
tmp = test.PPos;
tmp(tmp(:,1)==0,:) = 0;

stamps = [timestamp, tmp];
csvwrite([pwd '\Classification_algorithm\NH\InputData\gazeData.csv' ], stamps);
%Execute NH algorithm
mainRun;
%read solution
NH_output = csvread([pwd '\Classification_algorithm\NH\DetectionResults\gazeData.csv' ], 1);
%
idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        %{
        % Execute NH algorithm with one stimulus       
        
        deleteInputNHFiles()
        tmp = test.PPos(infIdx:supIdx,:);
        tmp(tmp(:,1)==0,:) = 0;

        stamps = [timestamp(infIdx:supIdx), tmp];
        csvwrite([pwd '\Classification_algorithm\NH\InputData\gazeData.csv' ], stamps);
        %Execute NH algorithm
        mainRun;
        %read solution
        NH_output = csvread([pwd '\Classification_algorithm\NH\DetectionResults\gazeData.csv' ], 1);

        %}
        
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'NH';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
        classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).Velocities = NH_output(infIdx:supIdx,5);         
        %classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx); 
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI);       
        %classificated_stimuli(idx_CS).GazePos = [NH_output(:,2),NH_output(:,3)]; %left eye
        %classificated_stimuli(idx_CS).Classification = NH_output(:,4); % 1 = fixation; 2 = saccade; 3 = glissade; 4 = nan
    
        classificated_stimuli(idx_CS).GazePos = [NH_output(infIdx:supIdx,2),NH_output(infIdx:supIdx,3)]; %left eye
        classificated_stimuli(idx_CS).Classification = NH_output(infIdx:supIdx,4); % 1 = fixation; 2 = saccade; 3 = glissade; 4 = nan
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end    
end

%% NH 
%{
% Execute NH with all gaze data
deleteInputNHFiles()
%Save the input file
tmp = test.PPos;
tmp(tmp(:,1)==0,:) = 0;

stamps = [timestamp, tmp];
csvwrite([pwd '\Classification_algorithm\NH\InputData\gazeData.csv' ], stamps);
%Execute NH algorithm
mainRun;
%read solution
NH_output = csvread([pwd '\Classification_algorithm\NH\DetectionResults\gazeData.csv' ], 1);
%}

idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        %
        % Execute NH algorithm with one stimulus       
        
        deleteInputNHFiles()
        tmp = test.PPos(infIdx:supIdx,:);
        tmp(tmp(:,1)==0,:) = 0;

        stamps = [timestamp(infIdx:supIdx), tmp];
        csvwrite([pwd '\Classification_algorithm\NH\InputData\gazeData.csv' ], stamps);
        %Execute NH algorithm
        mainRun;
        %read solution
        NH_output = csvread([pwd '\Classification_algorithm\NH\DetectionResults\gazeData.csv' ], 1);

        %}
       
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'NH';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
        classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI);       
        
        classificated_stimuli(idx_CS).Velocities = NH_output(:,5);  
        classificated_stimuli(idx_CS).GazePos = [NH_output(:,2),NH_output(:,3)]; %left eye
        classificated_stimuli(idx_CS).Classification = NH_output(:,4); % 1 = fixation; 2 = saccade; 3 = glissade; 4 = nan
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end    
end

%% ID

stim = compile_stim_info(test);

idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI        
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'ID';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = stim{idx}.sId;
        classificated_stimuli(idx_CS).StimPos = [stim{idx}.sPosX, stim{idx}.sPosY];
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
        if(stim{idx}.valid == 1)
            classificated_stimuli(idx_CS).Timestamp = stim{idx}.T*1000;
            classificated_stimuli(idx_CS).Velocities = stim{idx}.V;
            classificated_stimuli(idx_CS).GazePos = [stim{idx}.X, stim{idx}.Y]; %left eye        
            classificated_stimuli(idx_CS).Classification = stim{idx}.fixatLabels; % 1 = fixation
        else
            classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = zeros(supIdx-infIdx,1);
            classificated_stimuli(idx_CS).GazePos = zeros(supIdx-infIdx,1);
            classificated_stimuli(idx_CS).Classification = zeros(supIdx-infIdx,1);
        end
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end
end


%% Classificator EMD 
stamps = [ (1:test.NPoints)', timestamp, mm2deg(test.PPos,500), zeros(test.NPoints,2), test.RPsize, mm2deg(test.PPos, 500), zeros(test.NPoints,1), 4*(test.sampleValidityIgnoreRows==0), test.LPsize, mm2deg(test.SPos,500)];

%csvwrite([pwd '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt' ], stamps);
dlmwrite([pwd '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt' ], stamps, 'delimiter', '\t', 'precision', '%10.6f'); 


classificator = {8};
method_str={8};
used = zeros(8,1);
data = cell(8,12);
method_name = {' IVT';' IKF';' IDT';' IMST';' IHMM'; 'IVVT'; 'IVDT'; 'IVMT'};
%IVT
    classificator{1} = classificator_IVT_class;
    used(1) = 1;
    method_str{1} = '_ivt';
    classificator{1}.saccade_detection_threshold = 70;
    
%IKF
    classificator{2} = classificator_IKF_class;
    used(2) = 1;
    method_str{2} = '_ikf';
    classificator{2}.chi_threshold =    15;
    classificator{2}.window_size =      5;
    classificator{2}.deviation =        1000;

%IDT
    classificator{3} = classificator_IDT_class;
    used(3) = 1;
    method_str{3} = '_idt';
    classificator{3}.dispersion_duration_sec_threshold =    0.1; %s
    classificator{3}.dispersion_threshold =                 1.35; %deg

%IMST
    classificator{4} = classificator_IMST_class;
    used(4) = 1;
    method_str{4} = '_imst';
    classificator{4}.saccade_detection_threshold =  0.6; %deg
    classificator{4}.window_size =                  200;

%IHMM
    classificator{5} = classificator_IHMM_class;
    used(5) = 1;
    method_str{5} = '_ihmm';
    classificator{5}.saccade_detection_threshold =      70;%deg/s
    classificator{5}.viterbi_sample_size =              200;
    classificator{5}.baum_welch_reiterations_count =    5;
%{
%IVVT
    classificator{6} = classificator_IVVT_class;
    used(6) = 1;
    method_str{6} = '_ivvt';
    classificator{6}.saccade_detection_threshold = 70; %deg/s
    classificator{6}.fixation_detection_threshold = 20; %deg/s

%IVDT
    classificator{7} = classificator_IVDT_class;
    used(7) = 1;
    method_str{7} = '_ivdt';
    classificator{7}.saccade_detection_threshold = 70; %deg/s
    classificator{7}.dispersion_duration_sec_threshold = 0.1;
    classificator{7}.idt_dispersion_threshold = 1.35; %deg

%IVMT
    classificator{8} = classificator_IVMT_class;
    used(8) = 1;
    method_str{8} = '_ivmt';
    classificator{8}.saccade_detection_threshold = 70; %deg/s
    classificator{8}.window_duration_sec_threshold = 0.5; %s
    classificator{8}.san_agustin_threshold = 0.1;
%}
    
for i=1:5
    data{ i,1 } = 'N/A';
    data{ i,2 } = 'N/A';
    data{ i,3 } = 'N/A';
    data{ i,4 } = 'N/A';
    data{ i,5 } = 'N/A';
    data{ i,6 } = 'N/A';
    data{ i,7 } = 'N/A';
    data{ i,8 } = 'N/A';
    data{ i,9 } = 'N/A';
    data{ i,10} = 'N/A';
    data{ i,11} = 'N/A';
    data{ i,12} = 'N/A';
    if( used(i) ~= 0 )
        classificator{i}.debug_mode = 0;
        classificator{i}.input_data_name = '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt'
        classificator{i}.x_field = 8;
        classificator{i}.y_field = 9;
        classificator{i}.v_field = 11;
        classificator{i}.header_count = 0;
        classificator{i}.delta_t_sec =  1/60; %sample rate = 60Hz
        classificator{i}.sample_rate = 60;
        classificator{i}.fields_count = 14;

        classificator{i}.use_degree_data_filtering_X = 0;
        classificator{i}.use_degree_data_filtering_Y = 0;
        classificator{i}.minimal_allowed_X_degree = 0;
        classificator{i}.maximal_allowed_X_degree = 0;
        classificator{i}.minimal_allowed_Y_degree = 0;
        classificator{i}.maximal_allowed_Y_degree = 0;

        classificator{i}.read_data();
        if( classificator{i}.error_code == 0 )

            classificator{i}.eye_tracker_data_filter_degree_range();
            classificator{i}.classify();
            classificator{i}.eye_tracker_data_filter_degree_range();
            classificator{i}.merge_fixation_time_interval = 75;
            classificator{i}.merge_fixation_distance = 0.5;
            classificator{i}.merge_records();

            classification = EMD_classification(classificator{i}, timestamp);
        else
            errordlg(classificator{i}.error_message,'File reader error.');
        end
    end
    
    idx = 1;
    for k =0:test.NROUNDS-1
        for j = 1:test.NSTIMULI
            infIdx = stimuliIdx(idx);
            supIdx = stimuliIdx(idx+1)-1;
            classificated_stimuli(idx_CS).NHC = test.NHC;
            classificated_stimuli(idx_CS).Algorithm = method_name{i};
            classificated_stimuli(idx_CS).Round = k + 1;
            classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
            classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
            classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye
            classificated_stimuli(idx_CS).Classification = classification(infIdx:supIdx); % 1 = fixation
            idx = idx + 1; 
            idx_CS = idx_CS + 1;
        end
    end
end

%% Manual classification
%{
% Victoria

idx = 1;
name = strsplit(filename, '.');
file = strcat('ManualClassifications\V\', name(1), '_classificated.mat')
try
vClass = load(char(file));

onset = getOnset(vClass.Round) + getOnset(vClass.Stimulus);
ranges = find(onset);
sID = vClass.Stimulus(ranges(1:end-1));
sRound = vClass.Round(ranges(1:end-1));
for k =1:test.NROUNDS
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        
        if(sum(sID == test.SId(infIdx)) && sum(sRound == k))
            range = ranges(find(sID == test.SId(infIdx))): ranges(find(sID == test.SId(infIdx))+1)-2;
            classificated_stimuli(idx_CS).NHC = test.NHC;
            classificated_stimuli(idx_CS).Algorithm = 'Victoria';
            classificated_stimuli(idx_CS).Round = k;
            classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
            classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
            classificated_stimuli(idx_CS).Timestamp = vClass.Timestamp(range);
            classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye            
            classificated_stimuli(idx_CS).Classification = vClass.Classification(range);
            
        else            
            infIdx = stimuliIdx(idx);
            supIdx = stimuliIdx(idx+1)-1;
            classificated_stimuli(idx_CS).NHC = test.NHC;
            classificated_stimuli(idx_CS).Algorithm = 'Victoria';
            classificated_stimuli(idx_CS).Round = k;
            classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
            classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
            classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye            
            classificated_stimuli(idx_CS).Classification = zeros(supIdx-infIdx, 1);

        end
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end
end 

catch     
      errorMessage = sprintf('Error: Victoria classification does not exist:\n');
      uiwait(warndlg(errorMessage));
end

% Esther

idx = 1;
name = strsplit(filename, '.');
file = strcat('ManualClassifications\E\', name(1), '_classificated.mat')
try
vClass = load(char(file));

onset = getOnset(vClass.Round) + getOnset(vClass.Stimulus);
ranges = find(onset);
sID = vClass.Stimulus(ranges(1:end-1));
sRound = vClass.Round(ranges(1:end-1));
for k =1:test.NROUNDS
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        
        if(sum(sID == test.SId(infIdx)) && sum(sRound == k))
            range = ranges(find(sID == test.SId(infIdx))): ranges(find(sID == test.SId(infIdx))+1)-2;
            classificated_stimuli(idx_CS).NHC = test.NHC;
            classificated_stimuli(idx_CS).Algorithm = 'Esther';
            classificated_stimuli(idx_CS).Round = k;
            classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
            classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
            classificated_stimuli(idx_CS).Timestamp = vClass.Timestamp(range);
            classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye            
            classificated_stimuli(idx_CS).Classification = vClass.Classification(range);
            
        else            
            infIdx = stimuliIdx(idx);
            supIdx = stimuliIdx(idx+1)-1;
            classificated_stimuli(idx_CS).NHC = test.NHC;
            classificated_stimuli(idx_CS).Algorithm = 'Esther';
            classificated_stimuli(idx_CS).Round = k;
            classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
            classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
            classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
            classificated_stimuli(idx_CS).Velocities = test.V(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye            
            classificated_stimuli(idx_CS).Classification = zeros(supIdx-infIdx, 1);

        end
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end
end
catch     
      errorMessage = sprintf('Error: Esther classification does not exist:\n');
      uiwait(warndlg(errorMessage));
end
%}

return;

end

function onset = getOnset(array)
filtDer = [1; -1]; % x - x_1
avatarLostSamples = zeros(size(array,1),1);
onset = conv(filtDer, array) ~= 0;

end

function deleteInputNHFiles()
    % Specify the folder where the files live.
    myFolder = [pwd '\Classification_algorithm\NH\InputData'];
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isdir(myFolder)
      errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
      uiwait(warndlg(errorMessage));
      return;
    end
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder, '*.csv'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(myFolder, baseFileName);
      fprintf(1, 'Now deleting %s\n', fullFileName);
      delete(fullFileName);
    end
end

function V = computeLinVelocity(X,Y,T,dpi)
% X and Y in pixels, T in s; V is returned in mm/s
Npoints = size(X,1);

% I compute d_x^{i-1} (avg_x_1) first to avoid a loop over the vector
avg_x_1 = zeros(Npoints,1);
avg_x_1(1,1) = X(1,1);
avg_x_1(2:end,1) = X(1:end-1,1);
Dx = X - avg_x_1;

% I compute d_y^{i-1} (avg_y_1) first to avoid a loop over the vector
avg_y_1 = zeros(Npoints,1);
avg_y_1(1,1) = Y(1,1);
avg_y_1(2:end,1) = Y(1:end-1,1);
Dy = Y - avg_y_1;

% Convert Dx and Dy from pix to mm
Dx = pix2mm(Dx,dpi);
Dy = pix2mm(Dy,dpi);

% Once I have displacements in x and y (Dx, Dy), I compute distances
D = sqrt(Dx.*Dx + Dy.*Dy);

% And with D one can compute velocities dividing by time [mm/sec]

V = D./T;
V(1) = 0; % the first one has zero velocity (otherwise is NaN)

return;
end


%%
% i = 1 mean that we work with fixations
% i = 2 mean that we work with saccades
% i = 3 mean that we work with pursuits
% i = 4 mean that we work with noise (not used)
function c = EMD_classification(classificator, timestamp)
    c = zeros(size(timestamp));
    tmp_lists = {classificator.fixation_records, classificator.saccade_records, classificator.pursuit_records, classificator.noise_records};
% For every list from our set
%Not taking into account the noise record
    for i=1:3
        idx = size(tmp_lists{i},2);        
        for j=1:idx
            c(tmp_lists{i}{j}( :,3)) = i;
        end
    end
end

