function [f_success,fixatCentroids,fixatLabels,fixatExtremes,fixatIni,fixatEnd] = ff_computeFixations(params,X,Y,T,MINDURA,MAXDISP)    

% We will have a vector of length = NSamples, and tag 1 for
% fixation, 2 for saccade (and afterwards, 0 for nonvalid)

% 0. Set parameters

nonClassPoints = size(X,1); % non classified points
fixatLabels = zeros(size(X)); % 1 for points corresponding to a fixation, 2 for saccade
fixatOnsets = [];
fixatCentroids = [];
fixatExtremes = [];

% Pre-paso: Marcar inicio y fin eliminando primeros y ultimos 0.5 segundos
idx = 1;
while idx <= size(T,1) && T(idx)-T(1)<0.5
    idx = idx + 1;
end
fixatIni = idx;
idx = size(T,1);
while idx > 0 && (T(size(T,1))-T(idx))<0.5
    idx = idx - 1;
end
fixatEnd = idx;
clear idx

% 1. Create a window of min duration
idxIni = 1;
idxFin = idxIni + 1; % 1 more than idxIni
while nonClassPoints > 0 && idxFin <= size(T,1) && idxIni < size(T,1)
    % initialize wdw over points to cover duration threshold
    while (idxFin < size(T,1)) && (T(idxFin + 1) - T(idxIni) < MINDURA)
        idxFin = idxFin + 1;
    end
    XWdw = X(idxIni:idxFin);
    YWdw = Y(idxIni:idxFin);
    DWdw = mm2deg(max(XWdw(:)) - min(XWdw(:)), params.VD) + mm2deg(max(YWdw(:)) - min(YWdw(:)), params.VD); % [Salvucci and Goldberg 2000]
    % add additional points until dispersion threshold is reached
    if (DWdw < MAXDISP)
        while (DWdw < MAXDISP) && (idxFin < size(T,1))
            idxFin = idxFin + 1;
            XWdw = X(idxIni:idxFin);
            YWdw = Y(idxIni:idxFin);
            DWdw = mm2deg(max(XWdw(:))-min(XWdw(:)),params.VD) + mm2deg(max(YWdw(:))-min(YWdw(:)),params.VD);
        end
        % a fixation has been found - store it
        idxFin = idxFin - 1;
        nonClassPoints = nonClassPoints - (idxFin-idxIni+1);
        fixatOnsets = [fixatOnsets; idxIni];
        fixatLabels(idxIni:idxFin,:) = 1;
        %center = centroid([XWdw YWdw]);
        centroidX = mean(XWdw(:));
        centroidY = mean(YWdw(:));
        fixatCentroids = [fixatCentroids; centroidX centroidY];
        fixatExtremes = [fixatExtremes; idxIni idxFin];
        %disp(['Fixation found: ' num2str(idxIni) '-' num2str(idxFin)]);        
        % values for start of next iteration
        idxIni = idxFin + 1;
        idxFin = idxIni + 1; % 1 more than idxIni
    elseif min(XWdw) == 0
        % There is an invalid sample in the window. Consider the sample at
        % idxIni as invalid too since the information for saccades can't be
        % relevant unless all samples around it are valid.
        fixatLabels(idxIni, :) = 2;
        nonClassPoints = nonClassPoints - 1;
        idxIni = idxIni + 1;
    else
        fixatLabels(idxIni, :) = 2;
        nonClassPoints = nonClassPoints - 1;
        idxIni = idxIni + 1; 
    end
end
f_success = 1;