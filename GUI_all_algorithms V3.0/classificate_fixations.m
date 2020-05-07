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
if(f_emptyFile)
    classificated_stimuli = [];
    nRounds = 0;
    nStimuli = 0
    error = 1;
    return;
end

nRounds = test.NROUNDS;
nStimuli = test.NSTIMULI;
framerate = test.FRAMERATE;


velocities = mm2deg(test.V, test.params.VD);


%% BIT
bit_param = alg_parameters('BIT', framerate);
ID = zeros(test.NPoints,1);
ID(:) = test.NHC;
distance = zeros(test.NPoints,1);
distance(:) = bit_param.distance_to_screen;
rawVelocity = computeLinVelocity(test.RawLeftPos(:,1), test.RawLeftPos(:,2), test.T, test.DPI);
timestamp = test.T*1000; % seconds to ms
stamps = [ ID, timestamp, test.RawLeftPos(:,1), test.RawRightPos(:,1), test.RawLeftPos(:,2), test.RawRightPos(:,2), distance];
[output_id,output_total] = BIT(stamps, bit_param);

%% Split the results into stimuli (ouptut skip the first stamp?)
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
        classificated_stimuli(idx_CS).Velocities = velocities(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI);
        classificated_stimuli(idx_CS).GazePos = [pix2mm(output_total.stamps_fixations(infIdx:supIdx,3), test.DPI), pix2mm(output_total.stamps_fixations(infIdx:supIdx,5), test.DPI)]; %left eye
        classificated_stimuli(idx_CS).Classification = classifications(infIdx:supIdx) %If not 0 is a fixation
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end    
end


%% NH
NH_param = alg_parameters('NH', framerate);

idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        %        
        
        deleteInputNHFiles()
        tmp = test.PPos(infIdx:supIdx,:);
        tmp(tmp(:,1)==0,:) = 0;

        stamps = [timestamp(infIdx:supIdx), mm2pix(tmp, test.DPI)];
        csvwrite([pwd '\Classification_algorithm\NH\InputData\gazeData.csv' ], stamps, 1,0);
        %Execute NH algorithm
        mainRun(NH_param);
        %read solution
        NH_output = csvread([pwd '\Classification_algorithm\NH\DetectionResults\gazeData.csv' ], 1);
        
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'NH';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = test.SId(stimuliIdx(idx));
        classificated_stimuli(idx_CS).StimPos = test.SPos(stimuliIdx(idx),:);
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).Velocities = velocities(infIdx:supIdx); 
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI);       
        classificated_stimuli(idx_CS).GazePos = [tmp(:,1),tmp(:,2)]; %left eye
        classificated_stimuli(idx_CS).Classification = NH_output(:,4); % 1 = fixation; 2 = saccade; 3 = glissade; 4 = nan
    
        %classificated_stimuli(idx_CS).GazePos = [NH_output(infIdx:supIdx,2),NH_output(infIdx:supIdx,3)]; %left eye
        %classificated_stimuli(idx_CS).Classification = NH_output(infIdx:supIdx,4); % 1 = fixation; 2 = saccade; 3 = glissade; 4 = nan
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end    
end


%% DIVE
dive_param = alg_parameters('our_ID_IVT', framerate);
stim = compile_stim_info(test, dive_param);

idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI        
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'our_ID_IVT';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = stim{idx}.sId;
        classificated_stimuli(idx_CS).StimPos = [stim{idx}.sPosX, stim{idx}.sPosY];
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
        
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).Velocities = velocities(infIdx:supIdx);
        classificated_stimuli(idx_CS).GazePos = test.PPos(infIdx:supIdx, :);  

        if(stim{idx}.valid == 1)     
            classificated_stimuli(idx_CS).Classification = stim{idx}.fixatLabels; % 1 = fixation
            classificated_stimuli(idx_CS).Classification(stim{idx}.initialSaccade) = 2; %Saccade
            if j == test.NSTIMULI && i ==test.NROUNDS-1
                classificated_stimuli(idx_CS).Classification = classificated_stimuli(idx_CS).Classification(1:end-1);
            end
        else
            classificated_stimuli(idx_CS).Classification = zeros(supIdx-infIdx+1,1);
        end
    
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end
end

%% Classificator EMD 
stamps = [ (1:test.NPoints)', timestamp, mm2deg(test.PPos,650), zeros(test.NPoints,2), test.RPsize, mm2deg(test.PPos, 650), zeros(test.NPoints,1), 4*(test.sampleValidityIgnoreRows==0), test.LPsize, mm2deg(test.SPos,650)];

%csvwrite([pwd '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt' ], stamps);
dlmwrite([pwd '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt' ], stamps, 'delimiter', '\t', 'precision', '%10.6f'); 


classificator = {8};
method_str={8};
used = zeros(8,1);
data = cell(8,12);
method_name = {' IVT';' IKF';' IDT';' IMST';' IHMM'; 'IVVT'; 'IVDT'; 'IVMT'};
%IVT
    IVT_param = alg_parameters('IVT', framerate);
    classificator{1} = classificator_IVT_class;
    used(1) = 1;
    method_str{1} = '_ivt';
    classificator{1}.saccade_detection_threshold = IVT_param.velocityThreshold;
    
%IKF
    IKF_param = alg_parameters('IKF', framerate);
    classificator{2} = classificator_IKF_class;
    used(2) = 1;
    method_str{2} = '_ikf';
    classificator{2}.chi_threshold =    IKF_param.chi_threshold;
    classificator{2}.window_size =      IKF_param.window_size;
    classificator{2}.deviation =        IKF_param.deviation;

%IDT
    IDT_param = alg_parameters('IDT', framerate);
    classificator{3} = classificator_IDT_class;
    used(3) = 1;
    method_str{3} = '_idt';
    classificator{3}.dispersion_duration_sec_threshold =    IDT_param.dispersion_duration_sec_threshold; %s
    classificator{3}.dispersion_threshold =                 IDT_param.dispersion_threshold; %deg

%IMST
    IMST_param = alg_parameters('IMST', framerate);
    classificator{4} = classificator_IMST_class;
    used(4) = 1;
    method_str{4} = '_imst';
    classificator{4}.saccade_detection_threshold =  IMST_param.saccade_detection_threshold; %deg
    classificator{4}.window_size =                  IMST_param.window_size;

%IHMM
    IHMM_param = alg_parameters('IHMM', framerate);
    classificator{5} = classificator_IHMM_class;
    used(5) = 1;
    method_str{5} = '_ihmm';
    classificator{5}.saccade_detection_threshold =      IHMM_param.saccade_detection_threshold;%deg/s
    classificator{5}.viterbi_sample_size =              IHMM_param.viterbi_sample_size;
    classificator{5}.baum_welch_reiterations_count =    IHMM_param.baum_welch_reiterations_count;
   
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
        classificator{i}.input_data_name = '\Classification_algorithm\EMDv1.5\InputData\gazeData.txt';
        classificator{i}.x_field = 8;
        classificator{i}.y_field = 9;
        classificator{i}.v_field = 11;
        classificator{i}.header_count = 0;
        classificator{i}.delta_t_sec =  1/framerate; %sample rate = 60Hz
        classificator{i}.sample_rate = framerate;
        classificator{i}.fields_count = 14;

        classificator{i}.use_degree_data_filtering_X = 0;
        classificator{i}.use_degree_data_filtering_Y = 0;
        classificator{i}.minimal_allowed_X_degree = 0;
        classificator{i}.maximal_allowed_X_degree = 0;
        classificator{i}.minimal_allowed_Y_degree = 0;
        classificator{i}.maximal_allowed_Y_degree = 0;

        classificator{i}.read_data();
        if( classificator{i}.error_code == 0 )
            EMD_general_param = alg_parameters('EMD', framerate);
            classificator{i}.eye_tracker_data_filter_degree_range();
            classificator{i}.classify();
            classificator{i}.eye_tracker_data_filter_degree_range();
            classificator{i}.merge_fixation_time_interval = EMD_general_param.merge_fixation_time_interval;
            classificator{i}.merge_fixation_distance = EMD_general_param.merge_fixation_distance;
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
            classificated_stimuli(idx_CS).Velocities = velocities(infIdx:supIdx);
            classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
            classificated_stimuli(idx_CS).GazePos = pix2mm(test.RawLeftPos(infIdx:supIdx,:), test.DPI); %left eye
            classificated_stimuli(idx_CS).Classification = classification(infIdx:supIdx); % 1 = fixation
            idx = idx + 1; 
            idx_CS = idx_CS + 1;
        end
    end
end

%I2MC
I2MC_param = alg_parameters('I2MC', framerate);
fix = I2MC([datapath '\' filename], I2MC_param);

classification = zeros(length(timestamp),1);
for i = 1:length(fix.start)
    classification(fix.start(i):fix.end(i)) = 1;        
end

idx = 1;
for i =0:test.NROUNDS-1
    for j = 1:test.NSTIMULI        
        infIdx = stimuliIdx(idx);
        supIdx = stimuliIdx(idx+1)-1;
        classificated_stimuli(idx_CS).NHC = test.NHC;
        classificated_stimuli(idx_CS).Algorithm = 'I2MC';
        classificated_stimuli(idx_CS).Round = i + 1;
        classificated_stimuli(idx_CS).Stimulus = stim{idx}.sId;
        classificated_stimuli(idx_CS).StimPos = [stim{idx}.sPosX, stim{idx}.sPosY];
        classificated_stimuli(idx_CS).RawVelocities = rawVelocity(infIdx:supIdx);
        classificated_stimuli(idx_CS).RawGazePos = pix2mm(test.RawLeftPos(infIdx:supIdx, :), test.DPI); 
        
        classificated_stimuli(idx_CS).Timestamp = timestamp(infIdx:supIdx);
        classificated_stimuli(idx_CS).Velocities = velocities(infIdx:supIdx);
        classificated_stimuli(idx_CS).GazePos = test.PPos(infIdx:supIdx, :);  
        
        classificated_stimuli(idx_CS).Classification = classification(infIdx:supIdx); % 1 = fixation
        
        idx = idx + 1; 
        idx_CS = idx_CS + 1;
    end
end
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

