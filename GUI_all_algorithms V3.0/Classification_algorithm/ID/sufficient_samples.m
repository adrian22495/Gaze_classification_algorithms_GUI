function [ output ] = sufficient_samples( validityVector )
%SUFFICIENT_SAMPLES Returns true if validityVector has enough samples for
% calculating BCEA, KDE and CG.
% The purpose of this function is merely to have a single script to modify
% in case the definition of "enough samples" changes.
    if sum(validityVector) < 10
        output = false;
    else
        output = true;
    end
end

