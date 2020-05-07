function [saccadicPrecision] = ff_compute_SP(PPos, targetPos, targetSizeDeg, prevTargetPos, prevTargetSizeDeg, initialSaccade, fixatLabels)
%FF_COMPUTE_SP Computes saccadic precision as the ratio between the actual 
% amplitude of the saccade and the desired amplitude (straight line between the 
% first sample of the saccade and the stimulus.
% If 75% of the samples after the saccade are not within the stimulus area
% or the initialSaccade is empty, the result will be NaN.

saccadicPrecision = NaN;

%% Discard empty saccades. (Might need to discard saccades with very few samples)
if isempty(initialSaccade)
    % No valid saccade
    return
end

%% Check if the samples before the saccade are near the previous stimulus
hSide = deg2mm(prevTargetSizeDeg) / 2;
validradiusmm = 1.1 * sqrt(hSide^2 + hSide^2); 
validSamples = 0;
for i = 1:initialSaccade(1)-1
    % Consider fixations only
    if fixatLabels(i) == 1
        distance = norm(PPos(i, :) - prevTargetPos);
        if distance < validradiusmm
            validSamples = validSamples + 1;
        end
    end
end
% At least 3/4 (arbitrary) of the samples must be near the previous stimulus to 
% consider the error to be small enough to calculate saccadic precision
if validSamples / sum(fixatLabels(1:initialSaccade(1)-1) == 1) < 3/4
    % Discard this stimulus
    return
end

%% Check if the samples after the saccade are near the stimulus
% Allow valid samples to be 10% away from the radius
hSide = deg2mm(targetSizeDeg) / 2;
validradiusmm = 1.1 * sqrt(hSide^2 + hSide^2); 
validSamples = 0;
for i = initialSaccade(end):length(PPos)
    % Consider fixations only
    if fixatLabels(i) == 1
        distance = norm(PPos(i, :) - targetPos);
        if distance < validradiusmm
            validSamples = validSamples + 1;
        end
    end
end

% At least 3/4 (arbitrary) of the samples must be near the stimulus to consider 
% the error to be small enough to calculate saccadic precision
if validSamples / sum(fixatLabels(initialSaccade(end):length(PPos)) == 1) < 3/4
    % Discard this stimulus
    return
end

%% Compute saccadic precision 
% Saccadic precision is the ratio between the actual amplitude of the saccade 
% and the desired amplitude (straight line between the first sample of the 
% saccade and the stimulus.
distPrevStimulusStimulus = norm(targetPos - prevTargetPos);
distPrevStimulusFinal = norm(PPos(initialSaccade(end), :) - prevTargetPos);
saccadicPrecision = distPrevStimulusFinal / distPrevStimulusStimulus;

end

