function [velocities, distances] = compute_velocities(positions, target, sampleValid, T)
%COMPUTE_VELOCITIES Calculates the velocity of each sample in positions
%   Computes velocity and distance to the target of each given sample. If 
%   Any invalid samples will be ignored, computing the speed for the first
%   sample after a valid sample after them.

velocities = [];
distances = [];

X = positions(:,1);
Y = positions(:,2);
i = 2;
while i <= length(positions)
    while i <= length(positions) - 2 && sampleValid(i) == 0
        % Skip the next sample since speed can't be computed for it
        i = i + 2;
        continue;
    end
    % Calculate vector from the last sample, distance and velocity
    vector = [X(i) Y(i)] - [X(i-1) Y(i-1)];
    distanceTraveled = norm(vector);
    velocity = distanceTraveled / (T(i) - T(i-1));
    
    % direction = vector / distanceTraveled;
    
    % Calculate distance to the stimulus
    vector = target - [X(i) Y(i)];
    distance = norm(vector);
    % Store the new velocity and distance to the stimulus
    velocities = [velocities; velocity];
    distances = [distances; distance];
    i = i + 1;
end

end

