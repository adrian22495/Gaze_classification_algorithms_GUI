function [stimuli_test, avatar_test,f_emptyFile] = ff_load_and_setup(filename, datapath, fDiscardOutliers, f_rawData, f_EyeSel)
% =================================
% THIS NEEDS UPDATING -- TODO
% ff_load_and_setup_v2.m
%
% Loads the data from txt file, computes velocities for each point, and
% from them fixations and saccades.
% -----------
% In:  txt file (path, name)
% Out: T [s], V [mm/s], Vang [deg/s]
% Out: PPos info: PPos = [X Y] in [mm] of gaze position per timestamp
% from bottom left corner of screen
% Out: SPos info: SPos = [X Y] in [mm] of stimulus position per timestamp
% Out: SSize info: SSize in [mm], diam. of a circle centered at SPos, and
% SSizeDeg contains the same in [deg]
% Out: SId and onsets info: SId stores the stimulus id per timestamp; in
% onsets, timestamps in which stimulus appears are ==1, rest are ==0 (also
% gradSId, gradPosX and gradPosY are computed (related to onsets)
% Out: NSTIMULI, DPI, NROUNDS, STIMDURATION
% Out: sampleValidity: 1 if sample is valid, 0 otherwise
% -----------
% Created: 16/03/2016
% ??/??/2016: Added outlier rejection
% 22/08/2016: Converted it into a function (returns a struct, test, with
% all the info)
% 22/08/2016: Improved outlier rejection with tags
% =================================

%% Initial params
f_verbose = false;
ETSamplingFreq = 120; % Hz, eye tracker sampling frequency
test.filename = filename;
test.datapath = datapath;

params.f_Tablet = 1; % we will remove all of Cheves, so all the rest should be tablet

if nargin < 5
    f_EyeSel = 2; % meaning we should use both eyes by default
    if nargin < 4
        f_rawData = 1; % use raw data by default
        if nargin < 3
            fDiscardOutliers = 1; % 1 if outliers (see below) are to be deleted from the data, 0 otherwise
        end
    end
end
params.fDiscardOutliers = fDiscardOutliers;
params.f_rawData = f_rawData;
params.f_EyeSel = f_EyeSel;

% Compute f_appv1 from filename (by date)
% To see if appv1 or not (if first version of logs, or second one)
% TODO - Check what f_appv_first_header == 1 or 2 means, and in general
% check all the versions of the log that we have been dealing with
Ch = textscan(filename,'%s %s %s', 'Delimiter','_');
forDate = Ch{1,3}{1,1};
mm = str2double(forDate(5:6));
dd = str2double(forDate(7:8));
yy = str2double(forDate(1:4));
params.f_appv1 = 0;
params.f_appv_first_header = 0;
if (yy == 2016)
    if (mm <= 2) || ((mm == 3) && (dd < 17))
        params.f_appv1 = 1;
        if f_verbose
            disp(['ff_load_and_setup_v2: * This test was done with the old version of the log.']);
        end
    end
    if ((mm == 8) && (dd >= 2) && (dd <= 5)) || (mm == 10 || mm == 11) 
        params.f_appv_first_header = 1;
    end
end

if ((yy == 2016) && (mm == 12)) || (yy >= 2017)
    params.f_appv_first_header = 2;    
end
%Added to the header a new parameter to know the eye tracker used
if (((yy >= 2018) && (mm >= 7) && (dd >= 24)) | ((yy >= 2018) && (mm >= 8)))
    params.f_appv_first_header = 3; 
end
%Added validity info to the saamples
if (((yy >= 2018) && (mm >= 12) && (dd >= 4)) |  yy > 2018)
    params.f_appv_first_header = 4; 
end

NHC = str2double(Ch{1,2}{1,1});
forTestNr = Ch{1,1}{1,1};
TestNr = str2double(forTestNr(5:end));

test.date.dd = dd;
test.date.mm = mm;
test.date.yy = yy;
test.NHC = NHC;
test.TestNr = TestNr;

% Screen resolution
if params.f_Tablet
    test.RESX = 2160;
    test.RESY = 1440; % -TODO: eventually this will be in the header of the txt file
else
    test.RESX = 1600; %1200;
    test.RESY = 900; %672; % -TODO: eventually this will be in the header of the txt file
end

params.MAX_GAP_LENGTH = 0.075; % [seconds]

params.VD = 500; % [mm] distance from screen

% Filtering parameters
params.f_plotVerbose = 0; % 1 if we want all the plots from the filters and so on
params.filterSGolay = 0; % 1 if we want to filter using Salvitzky-Golay. See its specific parameters below.
params.SGolay_k = 3; % k parameter for S-Golay filter (default = 3)
params.SGolay_f = 5; % f parameter for S-Golay filter (default = 5)
params.filterMedian_supp = 3; % support for Median filter (default = 3)

test.params = params;


%% Load data from log file (cells: Ch1 & Cd)
% cannot use csvread cause there are numeric and non-numeric values
% M = csvread([datapath '\' filename], 4, 0);
if f_verbose
    fprintf('\t%s\\%s\n', datapath, filename);
end
fileID = fopen([datapath '\' filename]);
if fileID == -1
    fprintf('Error: can not read the file!\n');
    if f_verbose
        fprintf('ff_load_and_setup: File %s does not exist.',...
            ' Adding it to rejected_logs.txt (Code 3) and ',...
            'exiting function...\n', filename);
    end
    fileIDRej = fopen('./rejected_logs.txt','a');
    fprintf(fileIDRej,'%s - Code 3\n',filename);
    fclose(fileIDRej);
    f_emptyFile = 1;
    return;
end
[header1Data, header2, columnData] = readFromLog(fileID, params.f_appv1,...
    params.f_appv_first_header);
fclose(fileID);

%% Hardcode all the columns for the different data to extract them later (struct: cols)
cols.Tdata = 8;
cols.Xdata = 14;
cols.Ydata = 15;
cols.LXdata = 18;
cols.LYdata = 19;
cols.RXdata = 25;
cols.RYdata = 26;
cols.LPsize = 20;
cols.RPsize = 27;
cols.label = 11;
if params.f_rawData
    cols.Xdata = 12;
    cols.Ydata = 13;
    cols.LXdata = 16;
    cols.LYdata = 17;
    cols.RXdata = 23;
    cols.RYdata = 24;
end
if params.f_appv1
    cols.Xdata = cols.Xdata - 2;
    cols.Ydata = cols.Ydata - 2;
    cols.Tdata = cols.Tdata - 2;
    cols.LXdata = cols.LXdata - 2;
    cols.LYdata = cols.LYdata - 2;
    cols.RXdata = cols.RXdata - 2;
    cols.RYdata = cols.RYdata - 2;
    cols.LPsize = cols.LPsize - 2;
    cols.RPsize = cols.RPsize - 2;
    cols.label = cols.label - 2;
end

if (params.f_appv_first_header == 3 && strcmp(header1Data{size(header1Data,2)}, 'Tobii') || ...
    params.f_appv_first_header == 4 && strcmp(header1Data{size(header1Data,2)}, 'Tobii'))
    cols.Tdata = 7;
    cols.Xdata = 8;
    cols.Ydata = 9;
    cols.LXdata = 8;
    cols.LYdata = 9;
    cols.RXdata = 20;
    cols.RYdata = 21;
    cols.LPsize = 13;
    cols.RPsize = 25;
    cols.label = 11;    
end

if (params.f_appv_first_header == 4 && strcmp(header1Data{size(header1Data,2)}, 'Tobii'))
    eye_string = header1Data{17};
    if(strcmp(eye_string, 'right'))
        params.f_EyeSel = 1;
        f_EyeSel = 1;
    else
        params.f_EyeSel = 0;
        f_EyeSel = 0;
    end
    
end

scale_factor = 96/144;
if yy >=2018
    scale_factor = header1Data{16};
end

clear forDate forTestNr mm dd yy NHC Ch

%% Get the precise date from the first header data
dateCell = textscan(header1Data{3}{:}, '%d/%d/%d %d:%d:%d');
[test.date.dd, test.date.mm, test.date.yy, test.date.hh, test.date.min,...
    test.date.ss] = deal(dateCell{:});

%% Correct order of timestamp of data in Cd
[columnData] = correctTimeStamp(columnData,cols);

%% Store general data in test struct
NSTIMULI = double(header1Data{4});
NROUNDS = double(header1Data{5});
STIMDURATION = double(header1Data{6}); % ms

DPI = 216;
if params.f_appv_first_header == 2 || params.f_appv_first_header == 3 || params.f_appv_first_header == 4
    ISDURATIONVAR = header1Data{9};
    if strcmp(ISDURATIONVAR,'yes')
        if f_verbose
            disp(('ff_load_and_setup_v2: * Duration variable.'));
        end
    end
    test.ISDURATIONVAR = ISDURATIONVAR;
end
test.NSTIMULI = NSTIMULI;
test.NROUNDS = NROUNDS;
test.STIMDURATION = STIMDURATION; % ms
test.DPI = DPI;

%% Reject outliers (based on zero and based on labels)
% Rejection will imply a 0 in the test.sampleValidity vector (1 is ok sample)
% Also, if there are no good points left in the log file, a message is
% output and the function returns.

if (isempty(columnData{1,1}))
    if f_verbose
        fprintf('ff_load_and_setup_v2: File %s is empty (before outlier',...
            'rejection). Adding it to rejected_logs.txt (Code 1) and ',...
            'exiting function...\n', filename);
    end
    fileIDRej = fopen('./rejected_logs.txt','a');
    fprintf(fileIDRej,'%s - Code 1\n',filename);
    fclose(fileIDRej);
    f_emptyFile = 1;
    return;
else
    Npoints = size(columnData{1},1);
    test.NOutliersRejected = 0;
    if (params.fDiscardOutliers)
        % test.sampleValidity will be: 0 if not valid, 1 if valid
        % test.sampleValidityIgnoreRows will be: 0 if not valid, 1 if
        % valid. Contiene marcas correspondientes solo al avatar y rondas
        % con id 555 (usado en el calculo de onsets).
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points prior to outlier rejection.\n', Npoints);
        end
        params = discardLostGaze(columnData,params);
        test.sampleValidity = params.sampleValidity;
        test.sampleValidityIgnoreRows = test.sampleValidity;
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze.\n', sum(test.sampleValidity(:)));
        end
        sampleValidityAvatarRows = discardLastAvatarRows(columnData);
        test.sampleValidity = test.sampleValidity .* sampleValidityAvatarRows;
        test.sampleValidityIgnoreRows = test.sampleValidityIgnoreRows .* sampleValidityAvatarRows;
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze & Avatar.\n', sum(test.sampleValidity(:)));
        end
        sampleValidityOutliersZero = discardOutliersZero(columnData,cols, f_EyeSel);
        test.sampleValidity = test.sampleValidity .* sampleValidityOutliersZero;
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze & Avatar & Zero.\n', sum(test.sampleValidity(:)));
        end
        sampleValidityLTZero = discardOutliersLTZero(columnData,cols, f_EyeSel);
        test.sampleValidity = test.sampleValidity .* sampleValidityLTZero;
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze & Avatar & Zero & LT Zero.\n', sum(test.sampleValidity(:)));
        end
        sampleValidityGTRes = discardOutliersGTRes(columnData,cols, test.RESX, test.RESY, f_EyeSel);
        test.sampleValidity = test.sampleValidity .* sampleValidityGTRes;
        if f_verbose
            fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze & Avatar & Zero & LT Zero & GT Res.\n', sum(test.sampleValidity(:)));
        end
        %if the test wasn't realized with tobii
        if (~((params.f_appv_first_header == 3 || params.f_appv_first_header == 4) && strcmp(header1Data{size(header1Data,2)}, 'Tobii')))
            sampleValidityLabels = discardOutliersLabels(columnData{1,cols.label});
            test.sampleValidity = test.sampleValidity .* sampleValidityLabels;
            if f_verbose
                fprintf('ff_load_and_setup_v2: %d points after outlier rejection - LostGaze & Avatar & Zero & LT Zero & GT Res & Labels.\n', sum(test.sampleValidity(:)));
            end
        end
    end
    test.NOutliersRejected = Npoints - sum(test.sampleValidity(:));
end

test.NPoints = Npoints;
NvalidPoints = sum(test.sampleValidity(:));
nSamplesInFirstStimulus = test.STIMDURATION / 1000 * ETSamplingFreq; % 3*60=180
if ((Npoints < nSamplesInFirstStimulus) || (NvalidPoints < nSamplesInFirstStimulus)) % at least one stimulus complete
    if f_verbose
        disp(['ff_load_and_setup_v2: File ' filename ' is empty after outlier rejection (less than 1 sitmulus). Adding it to rejected_logs.txt (Code 2) and exiting function...']);
    end
    fileIDRej = fopen('./rejected_logs.txt','a');
    fprintf(fileIDRej,'%s - Code 2\n',filename);
    fclose(fileIDRej);
    stimuli_test = 0;
    avatar_test =0;
    f_emptyFile = 1;
    return;
end


%% Process timestamp data to create T vector (& store it in test struct)
T = columnData{cols.Tdata} - columnData{cols.Tdata}(1,1);
if params.f_appv_first_header
    % If the avatar is shown, the first instant will be the first stimulus
    % after it ends.
    [avatarRows,~] =  find(columnData{1,1} < 0);
    index = 1;
    if size(avatarRows,1)>0
        index = size(avatarRows,1);
    end
    for i=1:size(avatarRows,1)-1
        if (avatarRows(i)+1 ~= avatarRows(i+1) && columnData{1,1}(avatarRows(i)+1)~= -1000 ...
                && columnData{1,1}(avatarRows(i)+1)~= 555)
            index = avatarRows(i)+1;
        end
    end
    T = columnData{1,cols.Tdata} - columnData{1,cols.Tdata}(index,1);
    % To deal with 555 rounds, all samples after such rounds will see their
    % timestamp reduced by the duration of the rounds before it.
    if isfield(params,'invalidGaps')
        for i=2:size(params.invalidGaps,1)
            durT = T(params.invalidGaps(i,2)) - T(params.invalidGaps(i,1));
            for j=params.invalidGaps(i,2)+1:size(T,1)
                T(j) = T(j) - durT;
            end
        end
    end
    
end
% Convert time to seconds
if ((params.f_appv_first_header == 3 || params.f_appv_first_header == 4) &&strcmp(header1Data{size(header1Data,2)}, 'Tobii')) %Tobii
    T = T / 1000000; %microseconds to seconds
else
    T = T / 1000; %miliseconds to seconds
end
test.T = T; % sec.

%% Fill in gaps
% use params.MAX_GAP_LENGTH (in sec.) and T (in sec.)

% 1. Fill in gaps in which sampleValidity==0
% Fill in each eye separately [Tobii] and also X and Y separately [-CHECK]
[LX,SVLX] = fillInGaps(T, columnData{1,cols.LXdata}, test.sampleValidity, params.MAX_GAP_LENGTH);
[LY,SVLY] = fillInGaps(T, columnData{1,cols.LYdata}, test.sampleValidity, params.MAX_GAP_LENGTH);
[RX,SVRX] = fillInGaps(T, columnData{1,cols.RXdata}, test.sampleValidity, params.MAX_GAP_LENGTH);
[RY,SVRY] = fillInGaps(T, columnData{1,cols.RYdata}, test.sampleValidity, params.MAX_GAP_LENGTH);

% 2. Discard small groups of valid samples surrounded by invalid ones
[LX,SVLX] = discardSparseInvalid(T, LX, SVLX, params.MAX_GAP_LENGTH);
[LY,SVLY] = discardSparseInvalid(T, LY, SVLY, params.MAX_GAP_LENGTH);
[RX,SVRX] = discardSparseInvalid(T, RX, SVRX, params.MAX_GAP_LENGTH);
[RY,SVRY] = discardSparseInvalid(T, RY, SVRY, params.MAX_GAP_LENGTH);

test.sampleValidity = test.sampleValidity .* SVLX .* SVLY .* SVRX .* SVRY;

% When I fill in a sample, test.sampleValidity must be 1 for that sample.
if params.f_plotVerbose
    ss=1670;ee=1725;figure, plot(T(ss:ee),test.sampleValidity(ss:ee),'k*-'); title('sample validity before gap fill in');
    fprintf('Prior to filling in gaps we had %d points.\n', sum(test.sampleValidity(:)));
    aggr = SVLX .* SVLY .* SVRX .* SVRY;
    test.sampleValidity(aggr==1) = 1;
    fprintf('After filling in gaps we have %d points.\n', sum(test.sampleValidity(:)));
    ss=1670;ee=1725;figure, plot(T(ss:ee),test.sampleValidity(ss:ee),'k*-'); title('sample validity after gap fill in');
    % PRINT OUT NPOINTS RECOVERED AND SO ON....
    ss=1670; ee=1725; % filename = Test1_1118382_201604081118239446.txt
    figure, plot(T(ss:ee),columnData{1,cols.RYdata}(ss:ee),'b-o'); hold on; plot(T(ss:ee),RY(ss:ee),'r-o'); title('R Y data');
    grid on; legend('original','after fill in'); xlabel('Time [s]'); ylabel('Position [pix]');
end

%% Compute xavg and yavg
% At this point, samples may be outside the display's resolution.
% Those that persist after averaging will be discarded
X = zeros(size(RX)); Y = zeros(size(RY));
if f_EyeSel == 0 % left eye
    X = LX;
    Y = LY;
    test.sampleValidity = test.sampleValidity .* SVLX .* SVLY;
elseif f_EyeSel == 1 % right eye
    X = RX;
    Y = RY;
    test.sampleValidity = test.sampleValidity .* SVRX .* SVRY;
elseif f_EyeSel == 2 % both eyes
    for i=1:size(RX,1)
        % Discard sample when any value is invalid
        if SVRX(i) == 0 || SVRY(i) == 0 || SVLX(i) == 0 ||  SVLY(i) == 0
            X(i) = 0;
            Y(i) = 0;
        else
            X(i) = 1/2 .* (LX(i) + RX(i)); % [pix]
            Y(i) = 1/2 .* (LY(i) + RY(i)); % [pix]
        end
    end
else
    disp('ff_load_and_setup_v2: ERROR: You must choose one of the three modes for eye processing (PAUSING...).');
    pause;
end


% Set invalid values to zero in both X and Y
X(test.sampleValidity == 0) = 0;
Y(test.sampleValidity == 0) = 0;

% Fill new possible gaps
[X,SVX] = fillInGaps(T, X, test.sampleValidity, params.MAX_GAP_LENGTH);
[Y,SVY] = fillInGaps(T, Y, test.sampleValidity, params.MAX_GAP_LENGTH);
test.sampleValidity = SVX .* SVY;

% Since some invalid values were set to zero before (discarding values 
% outside the display's boundaries), consider zero values to be invalid too
test.sampleValidity(X == 0 | Y == 0) = 0;

% If there are less than 30 valid samples contained in this log, reject it
if (sum(test.sampleValidity) < 30)
    if f_verbose
        disp(['ff_load_and_setup_v2: File ' filename ' is empty after outlier rejection (less than 1 sitmulus). Adding it to rejected_logs.txt (Code 2) and exiting function...']);
    end
    fileIDRej = fopen('./rejected_logs.txt','a');
    fprintf(fileIDRej,'%s - Code 2\n',filename);
    fclose(fileIDRej);
    f_emptyFile = 1;
    return;
end

%% Filter X and Y data for noise
if params.filterSGolay
    % Salvitzky-Golay filter
    % Parameter Set chosen ((3,5), as opposed to (5,11))
    k = params.SGolay_k; % default = 3
    f = params.SGolay_f; % default = 5
    Xfs = filterNoiseSGolay(X,test.sampleValidity,k,f);
    Yfs = filterNoiseSGolay(Y,test.sampleValidity,k,f);
    if params.f_plotVerbose
        figure, plot(T(ss:ee),X(ss:ee),'b-o'); hold on; plot(T(ss:ee),Xfs(ss:ee),'r-*'); title(['X data - Salvitzky-Golay filtering - params.: k=' num2str(k) ' and f=' num2str(f)]);
        grid on; legend('original','after filtering'); xlabel('Time [s]'); ylabel('Position [pix]');
        figure, plot(T(ss:ee),Y(ss:ee),'b-o'); hold on; plot(T(ss:ee),Yfs(ss:ee),'r-*'); title(['Y data - Salvitzky-Golay filtering - params.: k=' num2str(k) ' and f=' num2str(f)]);
        grid on; legend('original','after filtering'); xlabel('Time [s]'); ylabel('Position [pix]');
    end
    X = Xfs;
    Y = Yfs;
    clear Xfs Yfs
end

% The filters can make some samples go out of range!
invalid = X < 0 | Y < 0 | X > test.RESX | Y > test.RESY;
X(invalid) = 0;
Y(invalid) = 0;
test.sampleValidity(invalid) = 0;

%% Compute velocities for each point (including filtered angular velocity Vang)
rw2zero = find(test.sampleValidity == 0);
V = computeLinVelocity(X,Y,T,test.DPI,rw2zero); % V is returned in [mm/s]
test.V = V;
Vang = computeAngVelocity(X,Y,T,params.VD,test.DPI,rw2zero); % Vang is returned in [deg/s]
test.Vang = Vang;
VangFilt = filterVelocityAvg3(Vang);
test.VangFilt = VangFilt;
clear V Vang VangFilt
if (size(columnData{1},1) ~= size(X,1)) || (size(test.sampleValidity,1) ~= size(X,1))
    if f_verbose
        disp('PROBLEMAS! in ff_load_and_setup_v2');
    end
end

%% Compute positions on screen (from bottom left corner) in mm => X and Y in PPos

% Convert X and Y from top left corner of screen to bottom left corner
% X = X; no need to change anything in X
Y = test.RESY - Y;

% %% Plot all samples with a rectangle showing the display for reference
% figure, hold on,
% title(sprintf('Zeros: %d', sum(X == 0)));
% plot([0 0 test.RESX test.RESX 0], [0 test.RESY test.RESY 0 0], 'k-');
% axis equal
% scatter(X(test.sampleValidity == 1), Y(test.sampleValidity == 1),'g.');
% scatter(X(test.sampleValidity == 0), Y(test.sampleValidity == 0),'r.');

% Convert X and Y from pix to mm
X = pix2mm(X,test.DPI); % [mm]
Y = pix2mm(Y,test.DPI); % [mm]

% Build PPos
PPos = [X Y]; % [mm] from bottom left corner
test.PPos = PPos;
test.RawLeftPos = [LX, test.RESY- LY];
test.RawRightPos = [RX, test.RESY - RY];
clear X Y PPos

%% Compute radius of stimuli in mm => SSize and SSizeDeg
% Stimulus diameter (or side of a square) is col. 2 (stimulus_size)
% Note: SSize (diameter) is already in mm, no need to convert it
Diam = columnData{1,2};
SSize = Diam;
SSizeDeg = mm2deg(Diam,params.VD);
test.SSize = SSize;
test.SSizeDeg = SSizeDeg;
clear SSize SSizeDeg Diam

%% Compute positions of stimulus on screen (from bottom left corner) in mm => X and Y in SPos
if params.f_appv1
    SPos  = NaN*ones(size(T,1),2);
else
    % X and Y are cols. 3 and 4 (stimulus_x and stimulus_y)
    X = columnData{1,3} / (scale_factor); % Offset correction
    Y = columnData{1,4} / (scale_factor); % Offset correction
    
    % Convert X and Y from top left corner of screen to bottom left corner
    % X = X; % no need to change anything in X
    Y = test.RESY - Y; % [pix]
    
    % Convert X and Y from pix to mm
    X = pix2mm(X, test.DPI) + test.SSize / 2;
    Y = pix2mm(Y, test.DPI) - test.SSize / 2;
    % Build SPos
    SPos = [X Y]; % [mm] from bottom left corner
end
test.SPos = SPos;
clear X Y SPos

%% Compute onsets of stimuli (not only id changes, also pos changes need to be taken into account)
SId = double(columnData{1,1});
filtDer = [1; -1]; % x - x_1

% in gradSId, timestamps in which stimulus appears are !=0, rest are =0
gradSId = conv(filtDer,SId);

if params.f_appv1
    onsets = (gradSId ~= 0);
    onsets(1) = 1; % otherwise the first one does not get counted. It does if f_appv1 = 0 because then it uses SPos in X and Y
else
    gradPosX = conv(filtDer,round(test.SPos(:,1)));
    gradPosY = conv(filtDer,round(test.SPos(:,2)));
    onsets = (gradSId ~= 0 | gradPosX ~= 0 | gradPosY ~= 0);
end

onsets = onsets(1:end-1); % from the derivative it gets 1 more
% Se ponen a cero los valores correspondientes al avatar y rondas 555 y
% previa que se ignoran
onsets(test.sampleValidityIgnoreRows == 0) = 0;
for i=3:size(test.sampleValidityIgnoreRows,1)
    if test.sampleValidityIgnoreRows(i)==0
        onsets(i-1) = 0;
        onsets(i-2) = 0;
        %onsets(i+1) = 0;
    end
end
% Se mira si hay dos 1's seguidos en el vector onsets, en ese caso se pone
% a cero el segundo.
idx = find(onsets~=0);
for i=1:size(idx,1)-1
    if idx(i)+1 == idx(i+1)
        onsets(idx(i+1)) = 0;
    end
end

test.SId = SId;
test.onsets = onsets;
clear SId onsets

%% Plot pupil size for both eyes as a function of time throughout the whole test (for blink)

test.LPsize = columnData{cols.LPsize};
test.RPsize = columnData{cols.RPsize};
if params.f_plotVerbose
    figure,hold on
    for i=1:(size(test.T,1)-1)
        if test.sampleValidity(i)
            plot(test.T(i:i+1),test.LPsize(i:i+1),'b.');
            plot(test.T(i:i+1),test.RPsize(i:i+1),'r.');
        else
            plot(test.T(i:i+1),test.LPsize(i:i+1),'c.');
            plot(test.T(i:i+1),test.RPsize(i:i+1),'m.');
        end
    end
    grid on; axis([min(test.T((1):(end-1))) max(test.T((1):(end-1))) 0 46]); title(['Pupil size (NHC = ' num2str(test.NHC) ') - whole test']); legend('left','right'); xlabel('Time [sec]'); ylabel('Pupil size [mm?]');
    %saveas(gcf,['./pSizeWhole_' filename(1:end-4) '.png'],'png');
end

%% End of function
clear T
f_emptyFile = 0;

[avatar_test, stimuli_test] = splitAvatarStimuli(test)

return;

function [Cd] = correctTimeStamp(Cd,cols)
% Funcion que mira si el orden de los valores del timestamp sigue un valor
% creciente. Sino permuta las filas de los datos para que lo siga.
%Recorrer timestamps
for i=2:size(Cd{1,cols.Tdata},1)
    % Si timestamp siguiente es menor que el anterior se intercambian
    % filas
    k = i;
    while (k > 1) && (Cd{1,cols.Tdata}(k) < Cd{1,cols.Tdata}(k-1))
        %Intercambio de filas de todos los datos
        for j=1:size(Cd,2)
            Cd{1,j}([k-1 k]) = Cd{1,j}([k k-1]);
        end
        k = k - 1;
    end
end
return;

function [columnHeader1Data, columnHeader2, columnData] = readFromLog(fileID, f_appv1, f_appv_first_header)

% Read and store header info
if f_appv1
    header1Data = fgetl(fileID);
    columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d', 'Delimiter',','); % -CHECK if in fact stimulus duration and DPI are integers
else
    header1  = fgetl(fileID);
    Ch1 = textscan(header1,'%s', 'Delimiter',',');
    header1Data = fgetl(fileID);
    if f_appv_first_header == 1
        columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d %s %d', 'Delimiter',','); % -CHECK if in fact stimulus duration and DPI are integers
    elseif f_appv_first_header == 2
        columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d %d %s %s %f %f %f %d %d %f', 'Delimiter',','); % -CHECK if in fact stimulus duration and DPI are integers
    elseif f_appv_first_header == 3
        columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d %d %s %s %f %f %f %d %d %f %s', 'Delimiter',','); % Eye tracker used at the end
    elseif f_appv_first_header == 4
        columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d %d %s %s %f %f %f %d %d %f %s %s %s', 'Delimiter',','); % Eye tracker used at the end
    else
        columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d', 'Delimiter',','); % -CHECK if in fact stimulus duration and DPI are integers
    end
    
end
if f_appv_first_header == 1
    fgetl(fileID); fgetl(fileID); fgetl(fileID); % Ignore duplicated header
end
header2 = fgetl(fileID);
columnHeader2 = textscan(header2, '%s', 'Delimiter',',');

% Read and store data
if f_appv1
    columnData = textscan(fileID, '%d %f %s %s %s %f %s %d %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter',','); %18 doubles al final
else
    if f_appv_first_header == 0
        columnData = textscan(fileID, '%d %f %f %f %s %s %s %f %s %d %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter',','); %18 doubles al final
    elseif (f_appv_first_header == 3 && strcmp(columnHeader1Data{size(columnHeader1Data,2)}, 'Tobii')) %Tobii
        columnData = textscan(fileID, '%d %f %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter',','); %26 doubles al final
    elseif (f_appv_first_header == 4 && strcmp(columnHeader1Data{size(columnHeader1Data,2)}, 'Tobii')) %Tobii
        columnData = textscan(fileID, '%d %f %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %f', 'Delimiter',','); %26 doubles al final + 6 booleanos
	else %f_appv_first_header == 1 || 2 || 3 and EyeTribe
        columnData = textscan(fileID, '%d %f %f %f %s %s %s %f %s %d %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d', 'Delimiter',','); %18 doubles y track status al final
    end
end
return;

function params = discardLostGaze(C,params)
% Marks as outliers samples where gaze is lost (id=555) and the previus
% round
rr2del = []; % will contain rows to be deleted
[rr2delLostGaze,~] = find(C{1,1}==555); % Find rows with id = 555
% Find the previous rows of rows with 555 id (discard rounds)
index = 1;
%Almacenamos índices iniciales y finales de rondas con id 555 y la
%previa
invalidGap = [];
while index < size(rr2delLostGaze,1)
    if(C{1,1}(rr2delLostGaze(index)-1,1) ~= 555) %Marca como no validos files de ronda previa
        id = C{1,1}(rr2delLostGaze(index)-1,1);
        indexRoundIgnore = rr2delLostGaze(index)-1;
        
        while C{1,1}(indexRoundIgnore,1) == id
            rr2del = [rr2del; indexRoundIgnore];
            indexRoundIgnore = indexRoundIgnore - 1;
            if indexRoundIgnore == 0
                break;
            end
        end
        index = index + 1; %Avanza hasta la siguiente zona distinta de mirada perdida
        while C{1,1}(rr2delLostGaze(index)-1,1) == 555 && index < size(rr2delLostGaze,1)
            index = index + 1;
        end
        invalidGap = [invalidGap; [indexRoundIgnore+1 rr2delLostGaze(index-1)]];
    end
end
params.invalidGaps = invalidGap;
rr2del = [rr2del; rr2delLostGaze];
validSamples = ones(size(C{1,1},1),1);
validSamples(rr2del) = 0;
params.sampleValidity = validSamples;
return;

function validSamples = discardAvatarRows(C)
% Marks as outliers samples where the avatar appears (id < 0)
[rr2del,~] =  find(C{1,1} < 0);
validSamples = ones(size(C{1,1},1),1);
validSamples(rr2del) = 0;
return;

function validSamples = discardLastAvatarRows(C)
% Detects the start and end index of the two avatar appears and marks as
% outliers the last one.
filtDer = [1; -1]; % x - x_1
avatarLostSamples = zeros(size(C{1,1},1),1);
avatarLostSamples(((C{1,1} == -1000) + (C{1,1} == 555)) > 0) =  1;
grad = conv(filtDer, avatarLostSamples);
aux = find(grad);

validSamples = ones(size(C{1,1},1),1);
% If there are two avatar appears
if(size(aux,1) >= 4)
    validSamples(aux(end-1):aux(end)-1) = 0;
end
return;

function validSamples = discardOutliersZero(C,col, f_EyeSel)
% Marks as outliers samples where some of the data is zero (see
% testForOutliers vector to see which data)
rr2del = []; % will contain rows to be deleted
% Cols. which will be tested in outlier rejection are given by testForOutliers
testForOutliers = [col.Xdata, col.Ydata, col.LXdata, col.LYdata, col.RXdata, col.RYdata];
testForOutliers = [col.LXdata, col.LYdata, col.RXdata, col.RYdata];
% Instead of marking as not valid those samples in which any of the
% columns of testForOutliers is zero (that is the commented code below)
% now we mark as not valid only samples in which both X data are zero
% or both Y data are zero (both means left and right eyes)

%UPDATE: we mark as not valid only samples from the eye that is going to be
%used
rr2del = [];
for i=1:size(C{1,1},1)
    if f_EyeSel == 0 %left eye
        rr = (C{1,col.LXdata}(i)==0) && (C{1,col.LYdata}(i)==0) ;
    elseif f_EyeSel == 1 %right eye            
        rr = (C{1,col.RXdata}(i)==0) && (C{1,col.RYdata}(i)==0) ;        
    else %both eye
        rr = ((C{1,col.LXdata}(i)==0) && (C{1,col.LYdata}(i)==0)) || ...
        ((C{1,col.RXdata}(i)==0) && (C{1,col.RYdata}(i)==0));
    end
    %rr = 0;
    %    rr = ((C{1,col.LXdata}(i)==0) && (C{1,col.RXdata}(i)==0)) || ...
    %   ((C{1,col.LYdata}(i)==0) && (C{1,col.RYdata}(i)==0));
    
    %         rr = ((C{1,col.RXdata}(i)==0) && (C{1,col.RYdata}(i)==0)) || ...
    %         ((C{1,col.LXdata}(i)==0) && (C{1,col.LYdata}(i)==0));
    if rr
        rr2del = [rr2del; i];
    end
end
%     for i = testForOutliers
%         [rr2delAdd,~] = find(C{1,i}==0);
%         rr2del = [rr2del; rr2delAdd];
%     end
validSamples = ones(size(C{1,1},1),1);
validSamples(rr2del) = 0;
return;

function validSamples = discardOutliersLTZero(C, col, f_EyeSel)
% Marks as outliers samples where some of the data is less than zero (see
% testForOutliers vector to see which data)
rr2del = []; % will contain rows to be deleted
% Cols. which will be tested in outlier rejection are given by testForOutliers
if f_EyeSel == 0 %left eye
    testForOutliers = [col.LXdata, col.LYdata];
elseif f_EyeSel == 1 %right eye        
    testForOutliers = [col.RXdata, col.RYdata];
else %both eye
    testForOutliers = [col.LXdata, col.LYdata, col.RXdata, col.RYdata];
end

for i = testForOutliers
    [rr2delAdd,~] = find(C{1,i}<0);
    rr2del = [rr2del; rr2delAdd];
end
validSamples = ones(size(C{1,1},1),1);
validSamples(rr2del) = 0;
return;

function validSamples = discardOutliersGTRes(C,col, ResX, ResY, f_EyeSel)
% Marks as outliers samples where some of the data is greater than the
% display's resolution (see testForOutliers vector to see which data)
rr2del = []; % will contain rows to be deleted
% Cols. which will be tested in outlier rejection are given by testForOutliers
if f_EyeSel == 0 %left eye
    testForOutliersX = [col.LXdata];
    testForOutliersY = [col.LYdata];
elseif f_EyeSel == 1 %right eye            
    testForOutliersX = [col.RXdata];
    testForOutliersY = [col.RYdata];
else %both eye
    testForOutliersX = [col.LXdata, col.RXdata];
    testForOutliersY = [col.LYdata, col.RYdata];
end

for i = testForOutliersX
    [rr2delAdd,~] = find(C{1,i}>ResX);
    rr2del = [rr2del; rr2delAdd];
end
for i = testForOutliersY
    [rr2delAdd,~] = find(C{1,i}>ResY);
    rr2del = [rr2del; rr2delAdd];
end

validSamples = ones(size(C{1,1},1),1);
validSamples(rr2del) = 0;
return;

function validSamples = discardOutliersLabels(LabelColumn)
% Marks as outliers samples in which the state_string is not 'STATE_TRACKING_GAZE | STATE_TRACKING_EYES | STATE_TRACKING_PRESENCE'
[rr2del,~] = find(~strcmp(LabelColumn,'STATE_TRACKING_GAZE | STATE_TRACKING_EYES | STATE_TRACKING_PRESENCE'));
validSamples = ones(size(LabelColumn,1),1);
validSamples(rr2del) = 0;
return;

function [Xf,newSV] = fillInGaps(T,X,SV,MAX_GAP_LENGTH)
% SV is sampleValidity (0 if not valid, 1 if valid)
T = T * 1000;
MAX_GAP_LENGTH = MAX_GAP_LENGTH * 1000;
if size(SV,1)~=size(X,1)
    disp('PROBLEMAS! en fillInGaps en ff_load_and_setup_v2');
    pause;
end
Xf = X;
newSV = SV;
if (~SV(1))
    Tgap = 0;
    idxStart = 1;
    tStart = T(1);
    xStart = X(1);
else
    Tgap = 0;
end
for i=2:size(X,1)
    if (~SV(i)) % in gap
        Tgap = Tgap + (T(i)-T(i-1));
        if (SV(i-1)) % start of gap
            idxStart = i;
            tStart = T(i-1); % timestamp and x value of last sample...
            xStart = X(i-1); % ...before the gap (Tobii I-VT filter doc.)
        end
    else        % not in gap
        if (~SV(i-1)) % end of gap
            idxEnd = i-1;
            tEnd = T(i); % timestamp of first sample after gap
            xEnd = X(i); % x value of first sample after gap
            if (Tgap <= MAX_GAP_LENGTH) % fill in
                % Compute a and b of linear function x = a * t + b
                a = (xEnd - xStart) / (tEnd - tStart);
                b = xEnd - a * tEnd;
                % Fill the gap with the interpolated values
                for j=idxStart:idxEnd
                    Xf(j) = a * T(j) + b;
                end
                newSV(idxStart:idxEnd) = 1;
            end
            Tgap = 0; % for the next gap
        end
    end
end
f_verbose = 0;
if f_verbose
    fprintf('We have filled in %d points.\n', -(sum(SV(:))-sum(newSV(:))));
end
return

function [X, SV] = discardSparseInvalid(T, X, SV, MAX_GAP_LENGTH)
%%DISCARDSPARSEINVALID Removes small groups of valid values that are
%%between invalid values.
%% Discard any valid samples that are surrounded by invalid samples
i = 2;
while i < length(SV)
    if SV(i) == 1 && SV(i-1) == 0
        start = i;
        while i < length(SV) &&  SV(i+1) == 1
            i = i + 1;
        end
        if T(i) - T(start) < MAX_GAP_LENGTH
            SV(start:i) = 0;
        end
    end
    i = i + 1;
end

return

function V = computeLinVelocity(X,Y,T,dpi,rw2zero)
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
% and set to 0 if one of the two elements in the compute of distances
% is a discard value
D = sqrt(Dx.*Dx + Dy.*Dy);
D(rw2zero(1:end-1,1) + 1) = NaN;
D(rw2zero) = NaN;

% And with D one can compute velocities dividing by time [mm/sec]

V = D./T;
V(1) = 0; % the first one has zero velocity (otherwise is NaN)

return;

function Vang = computeAngVelocity(X,Y,T,VD,dpi,rw2zero)
% X and Y in pixels, T in s, VD in mm; V is returned in deg/s
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
% and set to 0 if one of the two elements in the compute of distances
% is a discard value
D = sqrt(Dx.*Dx + Dy.*Dy);
D(rw2zero(1:end-1,1) + 1) = NaN;
D(rw2zero) = NaN;

% Compute also angular distances, and angular velocity (in [deg/sec])
Dang = mm2deg(D, VD); % -CHECK if using atan is ok in the general case

Vang = Dang./T;
Vang(1) = 0; % the first one has zero velocity (otherwise is NaN)

return;

function vF = filterVelocityAvg3(v)
% TO DO: the filter support (3) should be a parameter
vF = zeros(size(v));
vF(1) = 0.5*(v(1) + v(2));
vF(end) = 0.5*(v(end-1) + v(end));
for i=2:(length(v)-1)
    vF(i) = nanmean([v(i-1) v(i) v(i+1)]); %mean ignoring nan values
end

return


function [avatarTest, stimuliTest] = splitAvatarStimuli(test)
    % Detects the start and end index of the two avatar appears and marks as
% outliers the last one.
filtDer = [1; -1]; % x - x_1
avatarLostSamples = zeros(size(test.SId,1),1);
avatarLostSamples(((test.SId == -1000) + (test.SId == 555)) > 0) =  1;
grad = conv(filtDer, avatarLostSamples);
range = find(grad);
if(sum(avatarLostSamples)==0)
    avatarTest = {};
    stimuliTest = test;
    return;
end

avatarTest.filename = test.filename;
avatarTest.datapath = test.datapath;
avatarTest.date = test.date;
avatarTest.NHC = test.NHC;
avatarTest.TestNr = test.TestNr;
avatarTest.RESX = test.RESX;
avatarTest.RESY = test.RESY;
avatarTest.params = test.params;
avatarTest.ISDURATIONVAR = test.ISDURATIONVAR;
avatarTest.NSTIMULI = 1;
avatarTest.NROUNDS = 1;
avatarTest.STIMDURATION = 12000; %12 seconds
avatarTest.DPI = test.DPI;
avatarTest.NOutliersRejected = test.NOutliersRejected;
avatarTest.sampleValidity = test.sampleValidity(range(1):range(2)-1);
avatarTest.sampleValidityIgnoreRows = test.sampleValidityIgnoreRows(range(1):range(2)-1);
avatarTest.NPoints = range(2)-1;
avatarTest.T = test.T(range(1):range(2)-1) + (-1 * test.T(1));
avatarTest.V = test.V(range(1):range(2)-1);
avatarTest.Vang = test.Vang(range(1):range(2)-1);
avatarTest.VangFilt = test.VangFilt(range(1):range(2)-1);
avatarTest.PPos = test.PPos(range(1):range(2)-1,:);
avatarTest.RawLeftPos = test.RawLeftPos(range(1):range(2)-1,:);
avatarTest.RawRightPos = test.RawRightPos(range(1):range(2)-1,:);
avatarTest.SSize = test.SSize(range(1):range(2)-1);
avatarTest.SSizeDeg = test.SSizeDeg(range(1):range(2)-1);
avatarTest.SPos = test.SPos(range(1):range(2)-1,:);
avatarTest.SId = test.SId(range(1):range(2)-1);
avatarTest.onsets = test.onsets(range(1):range(2)-1);
avatarTest.onsets(1) = 1;
avatarTest.LPsize = test.LPsize(range(1):range(2)-1);
avatarTest.RPsize = test.RPsize(range(1):range(2)-1);

stimuliTest = test;
stimuliTest.T = stimuliTest.T;

stimuliTest.sampleValidity(range(1):range(2)-1) = 0;
stimuliTest.sampleValidityIgnoreRows(range(1):range(2)-1) = 0;
stimuliTest.onsets(range(1):range(2)-1) = 0;
return;
