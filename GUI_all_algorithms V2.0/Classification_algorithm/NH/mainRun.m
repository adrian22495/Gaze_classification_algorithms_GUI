%clear all, close all, clc

global ETparams

%% Set parameters
ETparams.screenSz = [2160 1440]; % pixels
ETparams.screenDim = [0.254, 0.196]; % meters
ETparams.viewingDist = 0.50; % meters
ETparams.samplingFreq = 60; % Hz
ETparams.blinkVelocityThreshold = 1000;             % if vel > 1000 degrees/s, it is noise or blinks
ETparams.blinkAccThreshold = 100000;               % if acc > 100000 degrees/s^2, it is noise or blinks
ETparams.peakDetectionThreshold = 100;              % Initial value of the peak detection threshold. deg/s

ETparams.minFixDur = 0.055; % in seconds // before 0.055
ETparams.minSaccadeDur = 0.03; % in seconds /////before 0.01


%% Set paths to input data and where to write output data
folder_in = [cd,'\Classification_algorithm\NH\InputData\']; % Must be csv-files organized as the example
                                % files (one per trial)
                                % Columns names are 
                                % time (optional column), x-coord (pix), y-coord (pix)
folder_out = [cd,'\Classification_algorithm\NH\DetectionResults\'];

if ~exist(folder_out, 'dir')
    mkdir(folder_out)
end


%% Detect events. 

files = dir([folder_in, '*.csv']);
for f = files'
    
    % Read the data from csv-file
    disp([folder_in, f.name])
    ETdata = txt2mat([folder_in, f.name]);
    ETparams.data = ETdata;
    
    % Detect events
    eventDetection()
    
    % Write data to csv-file
    write2csv(ETparams, [cd,'\Classification_algorithm\NH\DetectionResults\', f.name])

end


%% Save and plot data
%{
% Calculate basic parameters
mean(cat(1,ETparams.data.avgNoise))
mean(cat(1,ETparams.data.stdNoise))
mean(cat(1,ETparams.glissadeInfo.duration))
mean(cat(1,ETparams.saccadeInfo.duration))
mean(cat(1,ETparams.fixationInfo.duration))


plotResultsVel(ETparams,1,1)

%Plot histograms
% figure
% hist(cat(1,ETparams.glissadeInfo.duration),100)
figure
hist(cat(1,ETparams.saccadeInfo.duration),40)
xlabel('Saccade duration (s)'),ylabel('Number of saccades')
figure
hist(cat(1,ETparams.fixationInfo.duration),100)
xlabel('Fixation duration (s)'),ylabel('Number of fixations')
%}