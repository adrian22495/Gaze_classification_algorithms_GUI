function ETdata = txt2mat(filename)
%--------------------------------------------------------------------------
% Converts text files to matlab structure on the format
% ETdata(1, 1).X
% ETdata(1, 1).Y
% Where (1, 1) means participant 1 and trial 3
% Input:  filename - name of file where data are stored
          % must be file.csv containing the columns TRIAL_INDEX, gaze_x, gaze_y
%--------------------------------------------------------------------------

M = csvread(filename, 1);

% Check number of columns
if size(M,2) == 3
    ETdata(1, 1).T = M(:, 1)';
    ETdata(1, 1).X = M(:, 2)';
    ETdata(1, 1).Y = M(:, 3)';
elseif size(M,2) == 2
    ETdata(1, 1).T = ones(size(M,1),1);
    ETdata(1, 1).X = M(:, 1)';
    ETdata(1, 1).Y = M(:, 2)';
else
    error('Wrong format of input file')
end