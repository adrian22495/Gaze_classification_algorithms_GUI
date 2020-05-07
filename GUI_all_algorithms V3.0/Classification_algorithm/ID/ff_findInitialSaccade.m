function [saccadeSampleIds, peakVelocity] = ff_findInitialSaccade(stim, view_distance, velocityThreshold)
%FF_FINDINITIALSACCADE Finds the first saccade that is directed towards the
% stimulus. Only the first samples will be considered when looking for the
% saccade. If these first samples contain a certain ammount of invalid
% samples, or if there is no saccade directed at the stimulus in the range,
% the result will be an empty list. When the saccade is found, the result
% will be a list of all the samples that are a part of said saccade.

% A saccade is expected to happen some time after the stimulus onset
% This value is likely to change since it is completely arbitrary for now.
expectedTime = 0.800; % Expect the saccade 800ms after the onset (at most)

[~, lastSampleId] = min(abs(stim.T - stim.T(1) - expectedTime));
sampleTime = mean(gradient(stim.T)); % Average time between samples (should be close to 0.016 or 0.0083)

% Initialize return value
saccadeSampleIds = [];
peakVelocity = NaN;

MIN_CORRECT_DISTANCE = stim.sSizeDeg + 1; % degrees. 1 degree extra tolerance.
MIN_SACCADE_VELOCITY = velocityThreshold; % deg/s

% Don't try to find a saccade when there are more than 10% invalid samples
if sum(~stim.sampleValid(1:lastSampleId)) > lastSampleId * 0.1
    return;
end
invalidSampleIds = find(stim.sampleValid(1:lastSampleId) == 0);

% Gaze positions come in mm
X = stim.X;
Y = stim.Y;

%% Calculate distances and velocities
SPos = [stim.sPosX stim.sPosY];

[velocities, distances] = compute_velocities([X(1:lastSampleId), Y(1:lastSampleId)], ...
    SPos, stim.sampleValid(1:lastSampleId), stim.T(1:lastSampleId));

% % Consider filtering velocities and distances
% % velocities = sgolayfilt(velocities,3,5);
% % distances = sgolayfilt(distances,3,5);

% Convert mm into degrees 
velocities = mm2deg(velocities, view_distance);
distances = mm2deg(distances, view_distance);

if distances(1) < MIN_CORRECT_DISTANCE
    % The gaze begins over the stimulus, there can be no saccade
    return
end
if sum(findpeaks(velocities) > MIN_SACCADE_VELOCITY) > 5
    % The velocity has too many peaks, discard just to be sure
    return
end

% For even cleaner saccades, discard when velocities' average is over 40 deg/s.

%% Find the saccade. It should begin with a sample of a higher velocity and
%% lower distance than those before it

% Using the gradient of the distance we can find where the greatest approach happened
approachV = -gradient(distances) / sampleTime; % Velocity directed at the stimulus

% Find local maximums. Ideally, saccades will have a single velocity peak,
% however, children may perform two saccadic movements (could be more, but 
% we will consider those as failures for now).
[pks, ids] = findpeaks(approachV);
greaterPeaks = pks > MIN_SACCADE_VELOCITY;
peakIds = ids(greaterPeaks);

if isempty(peakIds)
    return
elseif length(peakIds) == 1 || peakIds(2) - peakIds(1) > 0.2 / sampleTime
    % Find the points that form the maximum peak
    startId = peakIds(1) - 1;
    endId = peakIds(1) + 1;
elseif length(peakIds) == 2
    startId = peakIds(1) - 1;
    endId = peakIds(2) + 1;
else
    % Too many peaks? Return a failure
    return
end

% Find the bases of the peak(s)
while startId > 1 && velocities(startId) > MIN_SACCADE_VELOCITY && velocities(startId) > velocities(startId - 1)
    startId = startId - 1;
end
while endId < length(velocities) && (velocities(endId) > MIN_SACCADE_VELOCITY || velocities(endId) > velocities(endId + 1))
    endId = endId + 1;
end

if min(distances(startId:endId)) > MIN_CORRECT_DISTANCE
    % The closest the gaze has gotten to the stimulus is not close enough
    % to consider it correct.
    return;
end

% Return the saccade Ids adjusting them by adding invalid positions to
% the left.
for i = 1:length(invalidSampleIds)
    if invalidSampleIds(i) <= startId
        startId = startId + 1;
        endId = endId + 1;
    elseif invalidSampleIds(i) <= endId
        endId = endId + 1;
    end
    if invalidSampleIds(i) < length(velocities)
        % Fill invalid ids for plot
        velocities = [velocities(1:invalidSampleIds(i) - 1); 0 ; velocities(invalidSampleIds(i):end)];
        approachV = [approachV(1:invalidSampleIds(i) - 1); 0 ; approachV(invalidSampleIds(i):end)];
        distances = [distances(1:invalidSampleIds(i) - 1); 0 ; distances(invalidSampleIds(i):end)];
    elseif length(velocities) < lastSampleId-1
        len = length(velocities);
        velocities = [velocities(1:min([invalidSampleIds(i)-1, len])); zeros(invalidSampleIds(i) - len, 1)];
        approachV = [approachV(1:min([invalidSampleIds(i)-1, len])); zeros(invalidSampleIds(i) - len, 1)];
        distances = [distances(1:min([invalidSampleIds(i)-1, len])); zeros(invalidSampleIds(i) - len, 1)];
    end
end

saccadeSampleIds = startId:endId;
peakVelocity = max(velocities);

%% Plot velocity and distance against time
if false
    zeroT = stim.T - stim.T(1);
    figure;
    subplot(3,1,1), hold on
    axis([0, expectedTime, 0, 200]);
    title(['Velocity/Time ', sprintf('%d, %3f', stim.sId, mean(velocities))]);
    plot(stim.T(1:lastSampleId-1) - stim.T(1), velocities)
    scatter(zeroT(1:lastSampleId-1) - zeroT(1), velocities, 60, 'r.')
    scatter([zeroT(startId), zeroT(endId)], [velocities(startId), velocities(endId)], 'g*');
    subplot(3,1,2), hold on
    axis([0, expectedTime, -50, 250]);
    title(['DistanceDecrement/Time ', sprintf('%d', stim.sId)]);
    plot(zeroT(1:lastSampleId-1) - zeroT(1), approachV)
    scatter(zeroT(1:lastSampleId-1), approachV, 60, 'r.')
    scatter([zeroT(startId), zeroT(endId)], [approachV(startId), approachV(endId)], 'g*');
    subplot(3,1,3), hold on
    axis([0, expectedTime, 0, 12]);
    title('Distance/Time');
    plot(zeroT(1:lastSampleId-1), distances)
    scatter(zeroT(1:lastSampleId-1), distances, 60, 'r.')
    scatter([zeroT(startId), zeroT(endId)], [distances(startId), distances(endId)], 'g*');
end

end
