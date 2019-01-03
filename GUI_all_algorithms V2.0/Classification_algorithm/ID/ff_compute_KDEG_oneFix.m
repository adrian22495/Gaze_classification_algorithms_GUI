function kdeg = ff_compute_KDEG_oneFix(params, XFix, YFix, sampleValidFix, StimX, StimY)
if nargin < 6 || StimX ~= StimX
    StimX = 0;
    StimY = 0;
end

kdeg.bw = -1;
kdeg.d = -1;
kdeg.x = -1;
kdeg.y = -1;
kdeg.area = -1;

if ~sufficient_samples(sampleValidFix)
    return;
end

f_verbose = false;
plotContours = false;

D = [XFix YFix];

% remove from fixation (XFix and YFix) non-valid samples
[rr2del,~] = find(sampleValidFix == 0);
D(rr2del,:) = [];

% convert x and y from mm to deg
D = mm2deg(D, params.VD);

%% Discard fixations with less than minSamples
if false %all(var(D) > 0)
    %% Compute density estimation
    [bandwidth, density, X, Y] = kde2d(D, 128);
    
    % Store outputs
    kdeg.bw = bandwidth;
    kdeg.d = density;
    kdeg.x = X;
    kdeg.y = Y;
    
    %% Compute the area in the three isoline 'heights':
    kdeg.area = zeros(1, size(params.P_FixStability, 2));

    % Discretize the 2D space, since the density is stimated in rectangles
    % we can get the area of a single one of this rectangles and multiply
    % it by the number of them that have a density greater than Y.
    rectangleArea = (X(1,2) - X(1,1)) * (Y(2,1) - Y(1,1));
    
    % For each threshold: [0.68 - 0.95 - 0.99]
    DensityInDataPoints = interpn(X', Y', density', D(:,1)', D(:,2)');
    prc = prctile(DensityInDataPoints, 100-params.P_FixStability * 100);
    
    for i=1:length(prc)
        % All density values greater than prc are within the desired area
        % The area is equal to the number of rectangles times their area
        kdeg.area(i) = rectangleArea * sum(sum(density > prc(i)));
    end

    if plotContours
        % Plot data and contours
        figure, hold on,
        title(sprintf('Areas: %.2f, %.2f, %.2f', kdeg.area));
        limDeg = [mm2deg(254, params.VD) mm2deg(169.3, params.VD)];
        set(gca,'YDir','Reverse')
        plot([0 0 limDeg(1) limDeg(1) 0], [0 limDeg(2) limDeg(2) 0 0], 'k-');
        axis equal;
        scatter(D(:,1), D(:,2),'r.');
        c = contourc(X(1,:), Y(:,1), density, 80);
        lines = getcontourlines(c);
        colors = {'g-', 'b-', 'r-'}; 
        for i=1:length(prc)
            [~, index] = min(abs([lines.v] - prc(i)));
            for j = 1:length(lines(index).x)
               plot(lines(index).x{j}, lines(index).y{j}, colors{i});
            end
        end
        scatter(mm2deg(StimX, params.VD), mm2deg(StimY, params.VD), 'g*');
    end
elseif f_verbose
     fprintf('compute_KDEG_oneFix: Rejected fixation\n');
end

return