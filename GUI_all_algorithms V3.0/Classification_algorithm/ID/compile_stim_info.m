function stim = compile_stim_info(test, extra_param)

f_verbose = false;
%% Create a cell with one item per stimuli. Each one will be a struct.
stim = cell(test.NSTIMULI*test.NROUNDS,1);

%% Start filling in the cell

%% 1. cycle number (cyc) [scalar]
k = 1;
for i=1:test.NROUNDS
    for j=1:test.NSTIMULI
        stim{k}.cyc = i;
        k = k + 1;
    end
end
clear i j k

%% 2. Order number (ord) [scalar]
k = 1;
for i=1:test.NROUNDS
    for j=1:test.NSTIMULI
        stim{k}.ord = j;
        k = k + 1;
    end
end
clear i j k

%% 3. Stimulus: number (sId), validity (valid), position (sPosX, sPosY) and size (sSizeDeg) [scalars]
auxSId = test.SId(test.onsets);
auxSSizeDeg = test.SSizeDeg(test.onsets);
auxSPosX = test.SPos(test.onsets,1);
auxSPosY = test.SPos(test.onsets,2);
onsetT = test.T(test.onsets);
dur = test.STIMDURATION/1000;

for i=1:size(stim,1)
    stim{i}.valid = 0;
end

for k=1:size(onsetT,1)
    idx = k;
    stim{idx}.sId = auxSId(k);
    stim{idx}.sSizeDeg = auxSSizeDeg(k);
    stim{idx}.sPosX = auxSPosX(k);
    stim{idx}.sPosY = auxSPosY(k);
    if mod(round(onsetT(k)),dur) == 0 % A stimulus starts in this sample
        stim{idx}.valid = 1;
    else
        % I really don't know what this is discarding... Gaze recovery
        % attempts?
        if f_verbose
            fprintf('Something is wrong with element %d of onsetT\n', k);
        end
    end
end
clear k i idx auxSId auxSSizeDeg auxSPosX auxSPosY

%% 4. X, Y, V, Vang, VangFilt, T, LPsize, RPsize, sampleValidity of each stimuli [vectors]
onsetIds = find(test.onsets==1);
for k=1:size(onsetT,1)-1 % size(stim,1)
    if mod(round(onsetT(k)),dur) == 0 % si multiplo de 3
        idx = k;
        if stim{idx}.valid && (idivide(int32(round(onsetT(k+1))),dur)+1-idx == 1)
            range = onsetIds(k):(onsetIds(k+1)-1);
            stim{idx}.X = test.PPos(range,1);
            stim{idx}.Y = test.PPos(range,2);
            stim{idx}.T = test.T(range,1);
            stim{idx}.V = test.V(range,1);
            stim{idx}.Vang = test.Vang(range,1);
            stim{idx}.VangFilt = test.VangFilt(range,1);
            stim{idx}.LPsize = test.LPsize(range,1);
            stim{idx}.RPsize = test.RPsize(range,1);
            stim{idx}.sampleValid = test.sampleValidity(range,1);
        else
            if f_verbose
                fprintf('Something is wrong with element %d of onsetT\n', k);
            end
        end
    end
end
% the last item
k = size(onsetT,1);
idx = k;
if stim{idx}.valid
    if(test.SId(end) == -1000)
        filtDer = [1; -1]; % x - x_1
        grad = conv(filtDer, test.sampleValidity);
        aux = find(grad);
        range = onsetIds(k):aux(end)-1;
    else
        range = onsetIds(k):size(test.PPos, 1);
    end
    stim{idx}.X = test.PPos(range,1);
    stim{idx}.Y = test.PPos(range,2);
    stim{idx}.V = test.V(range,1);
    stim{idx}.Vang = test.Vang(range,1);
    stim{idx}.VangFilt = test.VangFilt(range,1);
    stim{idx}.T = test.T(range,1);
    stim{idx}.LPsize = test.LPsize(range,1);
    stim{idx}.RPsize = test.RPsize(range,1);
    stim{idx}.sampleValid = test.sampleValidity(range,1);
elseif f_verbose
    fprintf('Something is wrong with element %d of onsetT\n', k);
end
clear k idx onsetIds onsetT
%% Some outlier rejection (i.e. labeling some stim as non-valid)
% If for some reason X Y V Vang or T fields are empty, stim is non valid
for i=1:size(stim,1)
    if stim{i}.valid
        if (~isfield(stim{i},'X') || ~isfield(stim{i},'Y') || ~isfield(stim{i},'T') || ...
                ~isfield(stim{i},'V') || ~isfield(stim{i},'Vang'))
            stim{i}.valid = 0;
        end
    end
end
% If less than 25% of samples are valid, stim is invalid
for i=1:size(stim,1)
    if stim{i}.valid
        if (sum(stim{i}.sampleValid(:)) < (0.25 * 180))
            stim{i}.valid = 0;
        end
    end
end
% % Sanity check of X
% figure,hold on;
% plot(stim{1}.X, 'r');
% plot(stim{2}.X, 'g');
% plot(stim{3}.X, 'b');
% plot(stim{4}.X, 'k');

for i=1:size(stim,1)
    if stim{i}.valid
        [stim{i}.initialSaccade, stim{i}.saccadicPeakV] = ff_findInitialSaccade(stim{i}, test.params.VD, extra_param.velocityThreshold);
    end
end

%% 5. Saccadic Reaction Time (SRT) per stimulus [scalar]
for i=1:size(stim,1)
    if stim{i}.valid
        [SRT] = ff_compute_SRT_oneStim(stim{i}.initialSaccade, test.FRAMERATE);
        stim{i}.SRT = SRT;
        stim{i}.TSRT = SRT + stim{i}.T(1);
    end
end
clear i SRT f_success
% Sanity check for SRT
if test.params.f_plotVerbose
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim{i}.valid
            figure, plot(stim{i}.T,stim{i}.X, 'b'); hold on;
            plot(stim{i}.T,stim{i}.Y, 'r');
            plot(stim{i}.T,stim{i}.sampleValid,'*m');
            line([stim{i}.TSRT stim{i}.TSRT],[min([stim{i}.X; stim{i}.Y]) max([stim{i}.X; stim{i}.Y])],'Color',[.75 .2 .75],'LineStyle','-','LineWidth',2);
        end
    end
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim{i}.valid
            figure, plot(stim{i}.T,stim{i}.VangFilt, 'b'); hold on;
            plot(stim{i}.T,stim{i}.Vang, 'r');
            plot(stim{i}.T,stim{i}.sampleValid,'*m');
            line([stim{i}.T(1) stim{i}.T(end)],[stim{i}.limitV stim{i}.limitV],'Color',[.75 .2 .75],'LineStyle','-','LineWidth',2);
        end
    end
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim{i}.valid
            figure, plot(stim{i}.T,stim{i}.V, 'r'); hold on;
            plot(stim{i}.T,stim{i}.sampleValid,'*m');
        end
    end
end


%% 6. Number of fixations (Nf) per stimuli [scalar] (as a side product, fixatLabels, fixatExtremes and fixatCentroids)
test.params.maxDispersion = extra_param.maxDispersion; % [deg]
test.params.minDuration = extra_param.minFixDuration; % [s]
test.params.maxTMergeFix = extra_param.maxSecondsToMergeFix; % [s], 75 ms.
test.params.maxSMergeFix = extra_param.maxDistanceToMergeFix; % [deg]
for i=1:size(stim,1)
    if stim{i}.valid
        [~,fixatCentroids,fixatLabels,fixatExtremes,fixatIni,fixatEnd] = ff_computeFixations(test.params,stim{i}.X,stim{i}.Y,stim{i}.T,test.params.minDuration,test.params.maxDispersion);
        NfAnt = size(fixatCentroids,1);
        [fixatCentroids,fixatLabels,fixatExtremes] = ff_mergeFixations(test.params,stim{i}.X,stim{i}.Y,stim{i}.T,fixatCentroids,fixatLabels,fixatExtremes,test.params.maxTMergeFix,test.params.maxSMergeFix);
        NfAct = size(fixatCentroids,1);
        while (NfAct < NfAnt)
            %fprintf('Next pass:\n');
            [fixatCentroids,fixatLabels,fixatExtremes] = ff_mergeFixations(test.params,stim{i}.X,stim{i}.Y,stim{i}.T,fixatCentroids,fixatLabels,fixatExtremes,test.params.maxTMergeFix,test.params.maxSMergeFix);
            NfAnt = NfAct;
            NfAct = size(fixatCentroids,1);
            %fprintf('===========\n');
        end
        % Before assigning it to stim{i} we should check that there is at
        % least one valid sample in each detected fixation (sampleValidity)
        initialNf = size(fixatCentroids,1);
        rr2del = []; % will contain rows to be deleted
        for j=1:initialNf
            sampleValidFix = stim{i}.sampleValid(fixatExtremes(j,1):fixatExtremes(j,2));
            nValidSamp = sum(sampleValidFix(:));
            if (nValidSamp <= 1) % I remove it from fixatCentroids, fixatLabels, fixatExtremes
                rr2del = [rr2del; j];
            end
        end
        fixatCentroids(rr2del,:) = [];
        fixatExtremes(rr2del,:) = [];
        fixatLabels(stim{i}.sampleValid==0)=0; % if non-valid it is neither sacc nor fixation, => =0
        stim{i}.Nf = size(fixatCentroids,1);
        stim{i}.fixatCentroids = fixatCentroids;
        stim{i}.fixatLabels = fixatLabels;
        stim{i}.fixatExtremes = fixatExtremes;
        stim{i}.fixatIni = fixatIni; % indices de inicio y fin eliminando primeros y últimos 0.5 segundos (see ff_computeFixations_v2)
        stim{i}.fixatEnd = fixatEnd;
    end
end
clear i fixatCentroids fixatLabels fixatExtremes rr2del sampleValidFix nValidSamp initialNf fixatIni fixatEnd

%% 6.5 Saccadic precision
stim{1}.SP = NaN;   % The gaze should already be on top of the first stimulus
for i=2:size(stim,1)
    if stim{i}.valid
        stim{i}.SP = ff_compute_SP([stim{i}.X, stim{i}.Y], [stim{i}.sPosX, stim{i}.sPosY],...
            stim{i}.sSizeDeg, [stim{i-1}.sPosX, stim{i-1}.sPosY], stim{i-1}.sSizeDeg, stim{i}.initialSaccade, stim{i}.fixatLabels);
    else
        stim{i}.SP = NaN;
    end
end

%% 7. Fixation duration (per fixation per stimuli) [vector per stimuli of size Nf]
for i=1:size(stim,1)
    if stim{i}.valid
        stim{i}.fixDur = zeros(stim{i}.Nf,1);
        for j=1:stim{i}.Nf
            stim{i}.fixDur(j) = stim{i}.T(stim{i}.fixatExtremes(j,2)) - stim{i}.T(stim{i}.fixatExtremes(j,1));
        end
    end
end
clear i j

%% 8. BCEA and KDEG (per fixation per stimuli) [vector per stimuli of size Nf]
test.params.P_FixStability = [0.68 0.95 0.99];
for i=1:size(stim,1)
    if stim{i}.valid
        stim{i}.BCEA = zeros(stim{i}.Nf,1);
        % Initialize variables in case stim{i}.Nf == 0
        stim{i}.CGV1 = ff_compute_CG_oneFixv1([],[],[],[]);
        stim{i}.CGV2 = ff_compute_CG_oneFixv2([],[],[],[]);

        StimX = stim{i}.sPosX;
        StimY = stim{i}.sPosY;
        
        for j=1:stim{i}.Nf
            XFix = stim{i}.X(stim{i}.fixatExtremes(j,1):stim{i}.fixatExtremes(j,2));
            YFix = stim{i}.Y(stim{i}.fixatExtremes(j,1):stim{i}.fixatExtremes(j,2));
            sampleValidFix = stim{i}.sampleValid(stim{i}.fixatExtremes(j,1):stim{i}.fixatExtremes(j,2));
            test.params.P_FixStability_type = 1; % 0.68
            stim{i}.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
            test.params.P_FixStability_type = 2; % 0.95
            stim{i}.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
            test.params.P_FixStability_type = 3; % 0.99
            stim{i}.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
            stim{i}.KDEG(j) = ff_compute_KDEG_oneFix(test.params,XFix,YFix,sampleValidFix, StimX, StimY);
            stim{i}.CGV1(j) = ff_compute_CG_oneFixv1(test.params,XFix,YFix,sampleValidFix);
            % Calcula el centro de gravedad con respecto de la posición del
            % estímulo solo si ésta es conocida.
            if i == 1
                stim{i}.CGV2(j) = ff_compute_CG_oneFixv2(test.params,XFix,YFix,sampleValidFix, StimX,StimY);
            end
        end
        XSampl = stim{i}.X(stim{i}.fixatIni:stim{i}.fixatEnd);
        YSampl = stim{i}.Y(stim{i}.fixatIni:stim{i}.fixatEnd);
        sampleValidSampl = stim{i}.sampleValid(stim{i}.fixatIni:stim{i}.fixatEnd);
        test.params.P_FixStability_type = 1; % 0.68
        stim{i}.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
        test.params.P_FixStability_type = 2; % 0.95
        stim{i}.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
        test.params.P_FixStability_type = 3; % 0.99
        stim{i}.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
        stim{i}.KDEGStim = ff_compute_KDEG_oneFix(test.params,XSampl,YSampl,sampleValidSampl, StimX,StimY);
        stim{i}.CGStimV1(1) = ff_compute_CG_oneFixv1(test.params,XSampl,YSampl,sampleValidSampl);
        stim{i}.CGStimV2(1) = ff_compute_CG_oneFixv2(test.params,XSampl,YSampl,sampleValidSampl,StimX,StimY);
    end
end
clear i j XFix YFix sampleValidFix dur StimX StimY XSampl YSampl sampleValidSampl
