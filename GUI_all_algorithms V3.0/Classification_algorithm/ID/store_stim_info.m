function [M, numColumns] = store_stim_info(fdatapath, ffilename)
%% Builds the whole results table from the logs

%% Save a row per log file, each containing 48*NMEAS+3 columns. [48 = NSTIM*NROUNDS]
% Col. 1 tiene NHC, col. 2 tiene validTest (0 if not valid, 1 if valid)
Names = {{'SRT'}; {'SP'}; {'NF'}; {'DF'};
    {'fBCEA68'}; {'fBCEA95'}; {'fBCEA99'}; {'fKDE68'}; {'fKDE95'}; {'fKDE99'};
    {'fCGV1'}; {'fCGV2'};
    {'sBCEA68'}; {'sBCEA95'}; {'sBCEA99'}; {'sKDE68'}; {'sKDE95'}; {'sKDE99'};
    {'sCGV1'}; {'sCGV2'}};
NMEAS = length(Names);

%% 0. Clean txt from repeated tests -- TO DO

%% 1. Count lines in file of filenames
fid = fopen([fdatapath '\' ffilename]);

tline = fgetl(fid);
cont = 1;
while ischar(tline)
    tline = fgetl(fid);
    cont = cont + 1;
end
fclose(fid);
cont = cont - 1; % hay un salto de linea vacio al final

%% 2. Create the matrix, full of -1s and with the header
Npac = cont; % one row more than patients, for header

M = repmat(-1, Npac, (1+48)*NMEAS+3);
numColumns = size(M,2);
%% Get results from the gaze that is averaged from both eyes
% This should vary with the patient, taking into account if they suffer
% from amblyopia
f_EyeSel = 0; % 0 for left eye, 1 for right eye, 2 for both eyes
f_rawData = true;
f_discard_outliers = true;
f_plot_eyes = false;

%% 3. Open files and fill in M
fid = fopen([fdatapath '\' ffilename]);
filenames = textscan(fid, '%s'); % Read all lines in one go
filenames = filenames{:}; % Unpack the cell of cells to a cell array
fclose(fid);

% Loop over filenames saving each result to its corresponding row
if (length(filenames) > 5)
    parfor i = 1:length(filenames)
        M(i, :) = load_and_compile(filenames{i}, fdatapath, f_discard_outliers,...
            f_rawData, f_EyeSel, f_plot_eyes);
    end
else
    for i = 1:length(filenames)
        M(i, :) = load_and_compile(filenames{i}, fdatapath, f_discard_outliers,...
            f_rawData, f_EyeSel, f_plot_eyes);
    end
end
