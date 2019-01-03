function cg = ff_compute_CG_oneFixv1(params,XFix,YFix,sampleValidFix)
if ~sufficient_samples(sampleValidFix)
    cg.centroid = [-1 -1];
    cg.distanceToCG = -1;
    cg.meanDistanceToCG = -1;
    cg.meanAngularLoc = -1;
    return;
end
% Cálculo según paper de Ygge 2005
% Cálculo del CG = media coordenadas X e Y posición mirada
% Obtención de la distancia de la mirada respecto al CG (distancia en deg y
% localización angular en deg)

D = [XFix YFix];

% remove from fixation (XFix and YFix) non-valid samples
[rr2del,~] = find(sampleValidFix == 0);
D(rr2del,:) = [];

% convert x and y from mm to deg
D = mm2deg(D, params.VD);

% Calculo CG media datos X - Y
XCentroid = mean(D(:,1));
YCentroid = mean(D(:,2));

%Calculo distancia de cada mirada al centroide
XDistanceToCG = abs(D(:,1) - XCentroid*ones(size(D,1),1));
YDistanceToCG = abs(D(:,2) - YCentroid*ones(size(D,1),1));
DistanceToCG = sqrt(XDistanceToCG.^2 + YDistanceToCG.^2);
AngularLocalization = rad2deg(atan(deg2rad(YDistanceToCG./XDistanceToCG)));

% Almacenaje de datos
cg.centroid = [XCentroid YCentroid];
cg.distanceToCG = DistanceToCG;
cg.meanDistanceToCG = mean(DistanceToCG);
cg.meanAngularLoc = mean(AngularLocalization);

% %% Plot CG
% figure, hold on,
% title('CGV1');
% scatter(D(:,1), D(:,2), 'b.');
% scatter(XCentroid, YCentroid, 'r*');

end