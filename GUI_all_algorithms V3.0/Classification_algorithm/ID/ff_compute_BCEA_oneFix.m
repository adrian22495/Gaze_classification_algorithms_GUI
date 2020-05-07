function bcea = ff_compute_BCEA_oneFix(params,XFix,YFix,sampleValidFix)
bcea = -1;
if ~sufficient_samples(sampleValidFix)
    return
end

D = [XFix YFix];

% remove from fixation (XFix and YFix) non-valid samples
[rr2del,~] = find(sampleValidFix == 0);
D(rr2del,:) = [];

% convert x and y from mm to deg
D = mm2deg(D, params.VD);

% Compute the area of the ellipse
sigmax = std(D(:,1));
sigmay = std(D(:,2));

% probability (0.63,0.95,0.68) - multiple values have been used in the literature
P = params.P_FixStability(params.P_FixStability_type); %[0.68 0.90 0.99]; 
k = -log(-P+1);
R = corrcoef(D);
if (size(R,1) == 2) && (size(R,2) == 2)
    rho = R(1,2);
    bcea = 2 * k * pi * sigmax * sigmay * (1-rho^2)^(1/2);
end

return