function [SRT] = ff_compute_SRT_oneStim(initialSaccade)

%% Default values
if isempty(initialSaccade)
    % There is not an easily differentiable saccade in this stimulus
    SRT = NaN;
else
    SRT = initialSaccade(1) * 1 / 60;
end
end
