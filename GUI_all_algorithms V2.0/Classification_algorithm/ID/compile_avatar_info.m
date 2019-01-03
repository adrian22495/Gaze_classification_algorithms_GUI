function stim = compile_stim_info(test, time_constriction)

f_verbose = false;
%% Create a cell with one item per stimuli. Each one will be a struct.
stim = cell(test.NSTIMULI*test.NROUNDS,1);

%% Start filling in the cell

%% 1. cycle and order number (cyc) [scalar]
stim.cyc = 1;
stim.ord = 1;

%% 3. Stimulus: number (sId), validity (valid), position (sPosX, sPosY) and size (sSizeDeg) [scalars]
auxSId = test.SId(test.onsets);
auxSSizeDeg = test.SSizeDeg(test.onsets);
auxSPosX = test.SPos(test.onsets,1);
auxSPosY = test.SPos(test.onsets,2);
onsetT = test.T(test.onsets);
dur = test.STIMDURATION/1000;

stim.valid = 0;

stim.sId = auxSId(1);
stim.sSizeDeg = auxSSizeDeg(1);
stim.sPosX = auxSPosX(1);
stim.sPosY = auxSPosY(1);
if mod(round(onsetT(1)),dur) == 0 % A stimulus starts in this sample
    stim.valid = 1;
else
    % I really don't know what this is discarding... Gaze recovery
    % attempts?
    if f_verbose
        fprintf('Something is wrong with element %d of onsetT\n', 1);
    end
end

clear k i idx auxSId auxSSizeDeg auxSPosX auxSPosY

%% 4. X, Y, V, Vang, VangFilt, T, LPsize, RPsize, sampleValidity of each stimuli [vectors]
onsetIds = find(test.onsets==1);
if stim.valid && (idivide(int32(round(test.T(end))),dur) == 1)
    stim.X = test.PPos(:,1);
    stim.Y = test.PPos(:,2);
    stim.T = test.T(:,1);
    stim.V = test.V(:,1);
    stim.Vang = test.Vang(:,1);
    stim.VangFilt = test.VangFilt(:,1);
    stim.LPsize = test.LPsize(:,1);
    stim.RPsize = test.RPsize(:,1);
    stim.sampleValid = test.sampleValidity(:,1);
else
    if f_verbose
        fprintf('Something is wrong with element %d of onsetT\n', k);
    end
end

clear k idx onsetIds onset

%% Some outlier rejection (i.e. labeling some stim as non-valid)
% If for some reason X Y V Vang or T fields are empty, stim is non valid
if stim.valid
    if (~isfield(stim,'X') || ~isfield(stim,'Y') || ~isfield(stim,'T') || ...
            ~isfield(stim,'V') || ~isfield(stim,'Vang'))
        stim.valid = 0;
    end
end
% If less than 25% of samples are valid, stim is invalid
if stim.valid
    if (sum(stim.sampleValid(:)) < (0.25 * 180))
        stim.valid = 0;
    end
end
% % Sanity check of X
% figure,hold on;
% plot(stim{1}.X, 'r');
% plot(stim{2}.X, 'g');
% plot(stim{3}.X, 'b');
% plot(stim{4}.X, 'k');

    if stim.valid
        [stim.initialSaccade, stim.saccadicPeakV] = ff_findInitialSaccade(stim);
    end

%% 5. Saccadic Reaction Time (SRT) per stimulus [scalar]
if stim.valid
    [SRT] = ff_compute_SRT_oneStim(stim.initialSaccade);
    stim.SRT = SRT;
    stim.TSRT = SRT + stim.T(1);
end
clear i SRT f_success
% Sanity check for SRT
if test.params.f_plotVerbose
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim.valid
            figure, plot(stim.T,stim.X, 'b'); hold on;
            plot(stim.T,stim.Y, 'r');
            plot(stim.T,stim.sampleValid,'*m');
            line([stim.TSRT stim.TSRT],[min([stim.X; stim.Y]) max([stim.X; stim.Y])],'Color',[.75 .2 .75],'LineStyle','-','LineWidth',2);
        end
    end
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim.valid
            figure, plot(stim.T,stim.VangFilt, 'b'); hold on;
            plot(stim.T,stim.Vang, 'r');
            plot(stim.T,stim.sampleValid,'*m');
            line([stim.T(1) stim.T(end)],[stim.limitV stim.limitV],'Color',[.75 .2 .75],'LineStyle','-','LineWidth',2);
        end
    end
    for i=1:test.NSTIMULI%9:test.NSTIMULI*2%1:test.NSTIMULI%
        if stim.valid
            figure, plot(stim.T,stim.V, 'r'); hold on;
            plot(stim.T,stim.sampleValid,'*m');
        end
    end
end


%% 6. Number of fixations (Nf) per stimuli [scalar] (as a side product, fixatLabels, fixatExtremes and fixatCentroids)
test.params.maxDispersion = 3; % [deg]
test.params.minDuration = 0.160; % [s]
test.params.maxTMergeFix = 0.075; % [s], 75 ms.
test.params.maxSMergeFix = 0.5; % [deg]
if stim.valid
    [~,fixatCentroids,fixatLabels,fixatExtremes,fixatIni,fixatEnd] = ff_computeFixations(test.params,stim.X,stim.Y,stim.T,test.params.minDuration,test.params.maxDispersion);
    NfAnt = size(fixatCentroids,1);
    [fixatCentroids,fixatLabels,fixatExtremes] = ff_mergeFixations(test.params,stim.X,stim.Y,stim.T,fixatCentroids,fixatLabels,fixatExtremes,test.params.maxTMergeFix,test.params.maxSMergeFix);
    NfAct = size(fixatCentroids,1);
    while (NfAct < NfAnt)
        %fprintf('Next pass:\n');
        [fixatCentroids,fixatLabels,fixatExtremes] = ff_mergeFixations(test.params,stim.X,stim.Y,stim.T,fixatCentroids,fixatLabels,fixatExtremes,test.params.maxTMergeFix,test.params.maxSMergeFix);
        NfAnt = NfAct;
        NfAct = size(fixatCentroids,1);
        %fprintf('===========\n');
    end
    % Before assigning it to stim we should check that there is at
    % least one valid sample in each detected fixation (sampleValidity)
    initialNf = size(fixatCentroids,1);
    rr2del = []; % will contain rows to be deleted
    for j=1:initialNf
        sampleValidFix = stim.sampleValid(fixatExtremes(j,1):fixatExtremes(j,2));
        nValidSamp = sum(sampleValidFix(:));
        if (nValidSamp <= 1) % I remove it from fixatCentroids, fixatLabels, fixatExtremes
            rr2del = [rr2del; j];
        end
    end
    fixatCentroids(rr2del,:) = [];
    fixatExtremes(rr2del,:) = [];
    fixatLabels(stim.sampleValid==0)=0; % if non-valid it is neither sacc nor fixation, => =0
    stim.Nf = size(fixatCentroids,1);
    stim.fixatCentroids = fixatCentroids;
    stim.fixatLabels = fixatLabels;
    stim.fixatExtremes = fixatExtremes;
    stim.fixatIni = fixatIni; % indices de inicio y fin eliminando primeros y últimos 0.5 segundos (see ff_computeFixations_v2)
    stim.fixatEnd = fixatEnd;
end
clear i fixatCentroids fixatLabels fixatExtremes rr2del sampleValidFix nValidSamp initialNf fixatIni fixatEnd

%% 6.5 Saccadic precision
stim.SP = NaN;   % The gaze should already be on top of the first stimulus


%% 7. Fixation duration (per fixation per stimuli) [vector per stimuli of size Nf]
if stim.valid
    stim.fixDur = zeros(stim.Nf,1);
    for j=1:stim.Nf
        stim.fixDur(j) = stim.T(stim.fixatExtremes(j,2)) - stim.T(stim.fixatExtremes(j,1));
    end
end

clear i j

%% 8. BCEA and KDEG (per fixation per stimuli) [vector per stimuli of size Nf]
test.params.P_FixStability = [0.68 0.95 0.99];
if stim.valid
    stim.BCEA = zeros(stim.Nf,1);
    % Initialize variables in case stim.Nf == 0
    stim.CGV1 = ff_compute_CG_oneFixv1([],[],[],[]);
    stim.CGV2 = ff_compute_CG_oneFixv2([],[],[],[]);

    if test.params.f_appv1 && i == 1
        % Even in the old version, the first stimulus is in the middle
        StimX = pix2mm(test.RESX / 2, test.DPI);
        StimY = pix2mm(test.RESY / 2, test.DPI);
    else
        StimX = stim.sPosX;
        StimY = stim.sPosY;
    end
    for j=1:stim.Nf
        XFix = stim.X(stim.fixatExtremes(j,1):stim.fixatExtremes(j,2));
        YFix = stim.Y(stim.fixatExtremes(j,1):stim.fixatExtremes(j,2));
        sampleValidFix = stim.sampleValid(stim.fixatExtremes(j,1):stim.fixatExtremes(j,2));
        test.params.P_FixStability_type = 1; % 0.68
        stim.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
        test.params.P_FixStability_type = 2; % 0.95
        stim.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
        test.params.P_FixStability_type = 3; % 0.99
        stim.BCEA(j,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XFix,YFix,sampleValidFix);
        stim.KDEG(j) = ff_compute_KDEG_oneFix(test.params,XFix,YFix,sampleValidFix, StimX, StimY);
        stim.CGV1(j) = ff_compute_CG_oneFixv1(test.params,XFix,YFix,sampleValidFix);
        % Calcula el centro de gravedad con respecto de la posición del
        % estímulo solo si ésta es conocida.
        if ~test.params.f_appv1 || i == 1
            stim.CGV2(j) = ff_compute_CG_oneFixv2(test.params,XFix,YFix,sampleValidFix, StimX,StimY);
        end
    end
    XSampl = stim.X(stim.fixatIni:stim.fixatEnd);
    YSampl = stim.Y(stim.fixatIni:stim.fixatEnd);
    sampleValidSampl = stim.sampleValid(stim.fixatIni:stim.fixatEnd);
    test.params.P_FixStability_type = 1; % 0.68
    stim.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
    test.params.P_FixStability_type = 2; % 0.95
    stim.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
    test.params.P_FixStability_type = 3; % 0.99
    stim.BCEAStim(1,test.params.P_FixStability_type) = ff_compute_BCEA_oneFix(test.params,XSampl,YSampl,sampleValidSampl);
    stim.KDEGStim = ff_compute_KDEG_oneFix(test.params,XSampl,YSampl,sampleValidSampl, StimX,StimY);
    stim.CGStimV1(1) = ff_compute_CG_oneFixv1(test.params,XSampl,YSampl,sampleValidSampl);
    stim.CGStimV2(1) = ff_compute_CG_oneFixv2(test.params,XSampl,YSampl,sampleValidSampl,StimX,StimY);
end
clear i j XFix YFix sampleValidFix dur StimX StimY XSampl YSampl sampleValidSampl
