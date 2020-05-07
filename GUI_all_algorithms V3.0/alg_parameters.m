function [ parameters ] = alg_parameters( alg_name, framerate)
%Get an object with the requiered parameters for the algorithm alg_name
x_res = 2160; %px
y_res = 1440; %px
screen_size_x = 0.254; %meters
screen_size_y =0.1693; %meters
distance_to_screen = 0.65; %meters
dpi = 216;

%Saccadic
velocity_threshold = 90; %deg/s  
%Fixation var
maxDispersion = 1.5; %deg
minFixDuration = 0.16; %s
maxSecondsToMergeTwoFix = 0.075; %s
maxDegToMergeTwoFix = 0.5; %deg

switch alg_name
    case 'BIT'
        parameters.min_samples = 5; %%%%%%%%%%%% original: 3              
        parameters.n_lost = 3;                   
        parameters.perc_control = 1-sqrt(.001);  
        parameters.filename = 'BIT_ouput';       
        parameters.max_x = x_res;                 
        parameters.max_y = y_res;  
        parameters.distance_to_screen = distance_to_screen*100;
        
    case 'NH'
        parameters.screenSz = [x_res y_res]; % pixels
        parameters.screenDim = [screen_size_x, screen_size_y]; % meters
        parameters.viewingDist = distance_to_screen; % meters
        parameters.framerate = framerate; % Hz
        parameters.blinkVelocityThreshold = 1000;             % if vel > 1000 degrees/s, it is noise or blinks
        parameters.blinkAccThreshold = 100000;               % if acc > 100000 degrees/s^2, it is noise or blinks
        parameters.peakDetectionThreshold = velocity_threshold;              % Initial value of the peak detection threshold. deg/s

        parameters.minFixDur = minFixDuration; % in seconds // before 0.055
        parameters.minSaccadeDur = 0.03; % in seconds 

    case 'our_ID_IVT'
        parameters.velocityThreshold = velocity_threshold;      
        parameters.maxDispersion = maxDispersion; %deg
        parameters.minFixDuration = minFixDuration; %s
        parameters.maxSecondsToMergeFix = maxSecondsToMergeTwoFix; %s
        parameters.maxDistanceToMergeFix = maxDegToMergeTwoFix; %deg

    case 'IVT'
        parameters.velocityThreshold = velocity_threshold;

    case 'IKF'
        parameters.chi_threshold = 15;
        parameters.window_size = 5;
        parameters.deviation = 1000;

    case 'IDT'
        parameters.dispersion_duration_sec_threshold = minFixDuration; %s
        parameters.dispersion_threshold = maxDispersion; %deg

    case 'IMST'
        parameters.saccade_detection_threshold = 1; %deg
        parameters.window_size = 20;

    case 'IHMM'
        parameters.saccade_detection_threshold = velocity_threshold;%deg/s
        parameters.viterbi_sample_size = 200;
        parameters.baum_welch_reiterations_count = 5;
        
    case 'EMD' %General parameters for IVT, IKF, IDT, IMST and IHMM
        parameters.merge_fixation_time_interval = maxSecondsToMergeTwoFix*1000; %ms
        parameters.merge_fixation_distance = maxDegToMergeTwoFix; %deg
        
    case 'I2MC'
        parameters.xres = x_res;                        % maximum value of horizontal resolution in pixels
        parameters.yres = y_res;                        % maximum value of vertical resolution in pixels
        parameters.framerate = framerate;               % sampling frequency of data (check that this value matches with values actually obtained from measurement!)

        % Variables for the calculation of visual angle
        % These values are used to calculate noise measures (RMS and BCEA) of
        % fixations. The may be left as is, but don't use the noise measures then.
        parameters.scrSz = [screen_size_x*100 screen_size_y*100];          % screen size in cm
        parameters.disttoscreen = distance_to_screen*100;                  % distance to screen in cm.
                
        % % K-MEANS CLUSTERING
        if framerate == 120
            parameters.windowtime = 0.2;                 % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
            parameters.steptime = 0.02;                  % time window shift (s) for each iteration. Use zero for sample by sample processing
            parameters.maxerrors = 100;    
            parameters.downsamples = [2 5];              % Default [2 5 10] for 300Hz
            parameters.downsampFilter = 1;               % use chebychev filter when downsampling? 1: yes, 0: no. requires signal processing toolbox. is what matlab's downsampling functions do, but could cause trouble (ringing) with the hard edges in eye-movement data
            parameters.chebyOrder = 6;                   % order of cheby1 Chebyshev downsampling filter, default is normally ok, as long as there are 25 or more samples in the window (you may have less if your data is of low sampling rate or your window is small

        elseif framerate == 60
            parameters.windowtime = 0.2;                 % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
            parameters.steptime = 0.02;                  % time window shift (s) for each iteration. Use zero for sample by sample processing
            parameters.maxerrors = 100;                  % maximum number of errors allowed in k-means clustering procedure before proceeding to next file
            parameters.downsamples = [2];                % Default [2 5 10] for 300Hz
            parameters.downsampFilter = 1;               % use chebychev filter when downsampling? 1: yes, 0: no. requires signal processing toolbox. is what matlab's downsampling functions do, but could cause trouble (ringing) with the hard edges in eye-movement data
            parameters.chebyOrder = 3;                   % order of cheby1 Chebyshev downsampling filter, default (6) is normally ok, as long as there are 25 or more samples in the window (you may have less if your data is of low sampling rate or your window is small
        end
        
        % % FIXATION DETERMINATION
        parameters.cutoffstd = 2;                    % number of standard deviations above mean k-means weights will be used as fixation cutoff
        parameters.maxMergeDist = mm2pix(deg2mm(maxDegToMergeTwoFix, distance_to_screen*1000),dpi); %Default 30px % maximum Euclidean distance in pixels between fixations for merging
        parameters.maxMergeTime = maxSecondsToMergeTwoFix*1000;             % Default 30ms % maximum time in ms between fixations for merging
        parameters.minFixDur = minFixDuration;                   % minimum fixation duration after merging, fixations with shorter duration are removed from output

        %% OPTIONAL VARIABLES
        % There are additional setting variables in I2MC.m for STEFFEN
        % INTERPOLATION 
        
end   

return; 
end

