function [row] = load_and_compile(filename, datapath, fDiscardOutliers, ...
    f_rawData, f_EyeSel, plot_eyes)
%%load_and_compile 
% Load a file and compile its information to a cell array. This cell array
% must comply with the cell matrix used in store_stim_info.m . As such,
% this function is highly dependant on it. This function is meant to make
% running each file in parallel easier using a parfor loop.

Names = {{'SRT'}; {'SP'}; {'NF'}; {'DF'}; 
    {'fBCEA68'}; {'fBCEA95'}; {'fBCEA99'}; {'fKDE68'}; {'fKDE95'}; {'fKDE99'};
    {'fCGV1'}; {'fCGV2'};
    {'sBCEA68'}; {'sBCEA95'}; {'sBCEA99'}; {'sKDE68'}; {'sKDE95'}; {'sKDE99'};
    {'sCGV1'}; {'sCGV2'}};
NMEAS = length(Names);

%% Load the file, cleaning and filtering values
[test, avatar_test, f_emptyFile] = ff_load_and_setup(filename, datapath,...
    fDiscardOutliers, f_rawData, f_EyeSel);

%% Create the row for the test log analysis
row = repmat(-1, 1, (1+48)*NMEAS+3);
if (~f_emptyFile)
    bee_stim = compile_stim_info(test);
    if(isempty(avatar_test))
        stim = bee_stim;
    else
        avatar_stim = compile_avatar_info(avatar_test);
        stim{1} = avatar_stim;
        stim(2:size(bee_stim,1)+1,:) = bee_stim;
    end
else
    % log file was empty before or after outlier rejection
    row(1) = test.NHC;
    row(3) = 0; % col.2: validPatient: 0 if not valid, 1 if valid
    return
end

if plot_eyes
    ploteyes(test);
end

% stimSaccades(stim);

%% Fill in the row, putting each value of stim in the right place
row(1) = test.NHC;
row(2) = 0;
row(3) = 1; % col.2: validPatient: 0 if not valid, 1 if valid

for j=1:size(stim,1)
    if stim{j}.valid && (stim{j}.cyc<=6) 
        cyc = stim{j}.cyc;
        sId = stim{j}.sId;
        SRT = stim{j}.SRT;
        if isnan(SRT)
            SRT = -1;
        end
        SP = stim{j}.SP;
        if isnan(SP)
            SP = -1;
        end
        Nf = stim{j}.Nf;
        fixDurProm = nanmean(stim{j}.fixDur);
        BCEAProm = 0; 
        KDEProm = 0;
        CGv1Prom = nanmean([stim{j}.CGV1.meanDistanceToCG]);
        CGv2Prom = nanmean([stim{j}.CGV2.meanDistanceToCG]);
        actualFixations = 0;
        for k=1:stim{j}.Nf
            if stim{j}.BCEA(k,1) ~= -1
                actualFixations = actualFixations + 1;
                BCEAProm = BCEAProm + stim{j}.BCEA(k,:);
                % If any fixation was invalid (NaN) -> KDEProm will be invalid
                KDEProm = KDEProm + stim{j}.KDEG(k).area;
            end
        end
        if actualFixations > 0
            BCEAProm = BCEAProm / actualFixations; % media de la magnitud en cuestion
            KDEProm = KDEProm / actualFixations;
        else
            BCEAProm = [-1 -1 -1];
            KDEProm = [-1 -1 -1];
        end
        BCEAStim = stim{j}.BCEAStim;
        KDEGStim = stim{j}.KDEGStim.area;
        CGv1Stim = stim{j}.CGStimV1.meanDistanceToCG;
        CGv2Stim = stim{j}.CGStimV2.meanDistanceToCG;
        % 3: first columns, 1: Matlab indexes from 1, 6: maximum numCycles
        if(sId == -1000)
            colIdx = 3 + 1;
        else
            colIdx = 3 + 1 + NMEAS +(NMEAS * sId * 6) + (NMEAS * (cyc-1)); 
        end
        if colIdx>((1 + 48)*NMEAS)
            disp('Error - check');
            pause;
        end
        % Place each value in its corresponding column
        
        %% Stimulus-wide statistics
        row(colIdx) = SRT;
        row(colIdx + 1) = SP;
        row(colIdx + 2) = Nf;
        row(colIdx + 3) = fixDurProm;
        
        %% Average fixation statistics
        for k=1:size(BCEAProm,2)
            row(colIdx + 4 + k -1) = BCEAProm(1,k);
        end

        for k=1:size(KDEProm,2)
            row(colIdx + 7 + k -1) = KDEProm(1,k);
        end
        
        row(colIdx + 10) = CGv1Prom;
        row(colIdx + 11) = CGv2Prom;
        
        %% Stimulus-wide statistics
        for k=1:size(BCEAStim,2)
            row(colIdx + 12 + k -1) = BCEAStim(1,k);
        end

        for k=1:size(KDEGStim,2)
            row(colIdx + 15 + k -1) = KDEGStim(1,k);
        end
        
        row(colIdx + 18) = CGv1Stim;
        row(colIdx + 19) = CGv2Stim; 
    end
end
