%% FIXATION DETECTION USING THE IDENTIFICATION BY 2-MEANS CLUSTERING (I2MC) ALGORITHM
% Description:
% The I2MC algorithm was designed to accomplish fixation detection in data
% across a wide range of noise levels and when periods of data loss may
% occur.
% 
% Cite as:
% Hessels, R.S., Niehorster, D.C., Kemner, C., & Hooge, I.T.C., (2016).
% Noise-robust fixation detection in eye-movement data - Identification by 
% 2-means clustering (I2MC). Behavior Research Methods.
% 
% For more information, questions, or to check whether we have updated to a
% better version, e-mail: royhessels@gmail.com / dcnieho@gmail.com. I2MC is
% available from www.github.com/royhessels/I2MC
% 
% Most parts of the I2MC algorithm are licensed under the Creative Commons
% Attribution 4.0 (CC BY 4.0) license. Some functions are under MIT 
% license, and some may be under other licenses.
% 
% Quick start guide for adopting this script for your own data:
% 1) Build an import function specific for your data (see importTobiiTX300
% for an example). 
% 
% 2) Change line 106 to use your new import function. The format should be:
%
% data.time for the timestamp
% data.left.X & data.left.Y for left gaze coordinates
% data.right.X & data.right.Y for right gaze coordinates
% data.average.X & data.average.Y for average gaze coordinates
% 
% You may provide coordinates from both eyes, only the left, only the
% right, or only the average.
% Gaze coordinates should be in pixels, timestamps should be in milliseconds
% 
% 3) Adjust the variables in the "necessary variables" section to match your
%    data
% 4) Run the algorithm
% 
% Note: Signal Processing Toolbox is required for the default downsampling
% procedure. If not available, set opt.downsampFilter to 0. This will use a
% different downsampling procedure.
% 
% Tested on MATLAB R2012a, R2014b & R2016a

function fix = I2MC(datapath, param)

%% INITIALIZE 
commandwindow;

%% NECESSARY VARIABLES

% General variables for eye-tracking data
opt.xres                        = param.xres; % maximum value of horizontal resolution in pixels
opt.yres                        = param.yres; % maximum value of vertical resolution in pixels
opt.missingx                    = -opt.xres; % missing value for horizontal position in eye-tracking data (example data uses -xres). used throughout functions as signal for data loss
opt.missingy                    = -opt.yres; % missing value for vertical position in eye-tracking data (example data uses -yres). used throughout functions as signal for data loss
opt.freq                        = param.framerate; % sampling frequency of data (check that this value matches with values actually obtained from measurement!)

% Variables for the calculation of visual angle
% These values are used to calculate noise measures (RMS and BCEA) of
% fixations. The may be left as is, but don't use the noise measures then.
opt.scrSz                       = param.scrSz; % screen size in cm
opt.disttoscreen                = param.disttoscreen; % distance to screen in cm.

% Folders
% Data folder should be structured by one folder for each participant with
% the eye-tracking data in textfiles in each folder.
folders.data                    = 'example data'; % folder in which data is stored (each folder in folders.data is considered 1 subject)
folders.output                  = 'output'; % folder for output (will use structure in folders.data for saving output)

% Plot results
do.plots                        = 0; % if set to 1, plot of fixation detection for each trial will be saved as png-file in output folder.
% the figures works best for short trials (up to around 20 seconds)

%% OPTIONAL VARIABLES
% The settings below may be used to adopt the default settings of the
% algorithm. Do this only if you know what you're doing. Uncomment the
% settings below and run the algorithm.

% % STEFFEN INTERPOLATION
% opt.windowtimeInterp            = 0.1;  % max duration (s) of missing values for interpolation to occur
% opt.edgeSampInterp              = 2;    % amount of data (number of samples) at edges needed for interpolation
% opt.maxdisp                     = opt.xres*0.2*sqrt(2); % maximum displacement during missing for interpolation to be possible
% 
% % K-MEANS CLUSTERING
opt.windowtime                  = param.windowtime;        % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
opt.steptime                    = param.steptime;       % time window shift (s) for each iteration. Use zero for sample by sample processing
opt.maxerrors                   = param.maxerrors;        % maximum number of errors allowed in k-means clustering procedure before proceeding to next file
opt.downsamples                 = param.downsamples;
opt.downsampFilter              = param.downsampFilter;          % use chebychev filter when downsampling? 1: yes, 0: no. requires signal processing toolbox. is what matlab's downsampling functions do, but could cause trouble (ringing) with the hard edges in eye-movement data
opt.chebyOrder                  = param.chebyOrder;     % order of cheby1 Chebyshev downsampling filter, default is normally ok, as long as there are 25 or more samples in the window (you may have less if your data is of low sampling rate or your window is small
% 
% % FIXATION DETERMINATION
opt.cutoffstd                   = param.cutoffstd; % number of standard deviations above mean k-means weights will be used as fixation cutoff
opt.maxMergeDist                = param.maxMergeDist; % maximum Euclidean distance in pixels between fixations for merging
opt.maxMergeTime                = param.maxMergeTime; % maximum time in ms between fixations for merging
opt.minFixDur                   = param.minFixDur; % minimum fixation duration after merging, fixations with shorter duration are removed from output

%% SET-UP FOLDERS

folders.func                = 'functions'; % folder for functions, will add to matlab path
addpath(genpath(folders.func));
if ~isdir(folders.output)
    mkdir(folders.output);
end

%% START ALGORITHM

% create textfile and open for writing fixations
fid = fopen(fullfile('Classification_algorithm','I2MC-1.1',folders.output,'allfixations.txt'),'w');
fprintf(fid,'FixStart\tFixEnd\tFixDur\tXPos\tYPos\tFlankedByDataLoss\tFraction Interpolated\tWeightCutoff\tRMSxy\tBCEA\tFixRangeX\tFixRangeY\tParticipant\tTrial\n');
       

%% IMPORT DATA

%[data.time,data.left.X,data.left.Y,data.right.X,data.right.Y] = importTobiiTX300(fullfile(folders.data,fold(e).name,file(f).name),1,[opt.xres opt.yres],opt.missingx,opt.missingy);
[data.time,data.left.X,data.left.Y,data.right.X,data.right.Y] = importDIVEdata(datapath,opt.missingx,opt.missingy);

%% RUN FIXATION DETECTION
fix = I2MCfunc(data,opt);

%% PLOT RESULTS

if do.plots
    % pre-allocate name for saving file
    savefile = fullfile('Classification_algorithm','I2MC-1.1',folders.output, 'plot');
    plotResults(data,fix,savefile,[opt.xres opt.yres]);
end

for g=1:numel(fix.start)
    fprintf(fid,'%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n',[fix.startT(g) fix.endT(g) fix.dur(g) fix.xpos(g) fix.ypos(g) fix.flankdataloss(g) fix.fracinterped(g) fix.cutoff, fix.RMSxy(g), fix.BCEA(g), fix.fixRangeX(g), fix.fixRangeY(g)]);
end


%% CLEAN UP

fclose(fid);
end