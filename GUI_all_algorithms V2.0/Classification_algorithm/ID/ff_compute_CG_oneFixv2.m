function cg = ff_compute_CG_oneFixv2(params,XFix,YFix,sampleValidFix,StimX,StimY)
if ~sufficient_samples(sampleValidFix)
    cg.centroid = [-1 -1];
    cg.distanceToCG = -1;
    cg.meanDistanceToCG = -1;
    cg.meanAngularLoc = -1;
    return;
end

% Cálculo según paper de Kallmark 2005
% Cálculo del Centro de Gravedad = media coordenadas X e Y posición mirada relativa al
% estímulo
% Obtención de la distancia del CG respecto al estímulo (distancia en deg y
% localización angular en deg (hay que pasarlo a 30deg - bins)

D = [XFix YFix];

% remove from fixation (XFix and YFix) non-valid samples
[rr2del,~] = find(sampleValidFix == 0);
D(rr2del,:) = [];

% convert x and y from mm to deg
D = mm2deg(D, params.VD);

StimX = mm2deg(StimX, params.VD);
StimY = mm2deg(StimY, params.VD);

% I don't think this next part makes any sense whatsoever...
% %Posiciones mirada relativas al estimulo
% dataXStim = abs(D(:,1) - StimX*ones(size(D,1),1));
% dataYStim = abs(D(:,2) - StimY*ones(size(D,1),1));

% Calculo CG: media datos X - Y
XCentroid = mean(D(:,1));
YCentroid = mean(D(:,2));

% Calculo distancia entre el estimulo y el CG en grados
DistanceCGStimX = abs(XCentroid - StimX);
DistanceCGStimY = abs(YCentroid - StimY);
DistanceToCG = sqrt(DistanceCGStimX.^2 + DistanceCGStimY.^2);
% Calculo distancia angular en grados
AngularLocalization = rad2deg(atan(deg2rad(DistanceCGStimY./DistanceCGStimX)));

% Almacenaje de datos
cg.centroid = [XCentroid YCentroid];
cg.distanceToCG = DistanceToCG;
cg.meanDistanceToCG = DistanceToCG;
cg.meanAngularLoc = AngularLocalization;

% %% Plot CG
% figure, hold on,
% title('CGV2');
% axis equal
% axis([0 260 0 160]); % Make the plot as big as the original device's display
% scatter(deg2mm(D(:,1)), deg2mm(D(:,2)),'b.');
% scatter(deg2mm(cg.centroid(1)), deg2mm(cg.centroid(2)), 'r*');
% scatter(deg2mm(StimX), deg2mm(StimY), 'g*');

end