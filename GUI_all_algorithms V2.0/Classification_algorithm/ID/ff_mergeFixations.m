function [fixatCen,fixatLab,fixatExt] = ff_mergeFixations(params,X,Y,T,fixatCentroids,fixatLabels,fixatExtremes,MAXSEC,MAXDEG)    

% T should be the same size as fixatLabels
% fixatCentroids, X and Y are in mm (see ff_computeFixations and
% ff_load_and_setup)

NfIni = size(fixatCentroids,1);
fixatLab = fixatLabels;
if (NfIni==1) % si solo hay una la copio tal cual y punto
    fixatCen = fixatCentroids;
    fixatExt = fixatExtremes;
else
    fixatCen = [];
    fixatExt = [];
    i=1;
    while (i < NfIni)
        cA = fixatCentroids(i,:);
        cB = fixatCentroids(i+1,:);
        distSAB = sqrt((cA(1)-cB(1)).^2 + (cA(2)-cB(2)).^2); % 0.5 deg max entre los centroides
        distSAB = mm2deg(distSAB,params.VD);
        distTAB = T(fixatExtremes(i+1,1)) - T(fixatExtremes(i,2)); % 75 ms max entre el inicio de la segunda y el final de la primera
        if (distSAB < MAXDEG && distTAB < MAXSEC) % merge i and i+1
            % Label properly (mark as fixation all the gap in between)
            fixatLab(fixatExtremes(i,2):fixatExtremes(i+1,1)) = 1; 
            % Compute new centroid
            newX = X(fixatExtremes(i,1):fixatExtremes(i,2));
            newX = [newX; X(fixatExtremes(i+1,1):fixatExtremes(i+1,2))];
            newY = Y(fixatExtremes(i,1):fixatExtremes(i,2));
            newY = [newY; Y(fixatExtremes(i+1,1):fixatExtremes(i+1,2))];
            newCX = mean(newX);
            newCY = mean(newY);
            fixatCen = [fixatCen; newCX newCY];
            clear newX newY newCX newCY
            % Join extremes
            fixatExt = [fixatExt; fixatExtremes(i,1) fixatExtremes(i+1,2)];
            i = i+2;
        else % do not merge i and i+1 (we will try in next iter with i+1 and i+2
            % - fixatLab no hay que tocarlo
            % - fixatCen y fixatExt hay que copiar el i, que se queda igual (el
            % i+1 aún no sabemos si lo mergearemos)
            fixatCen = [fixatCen; fixatCentroids(i,:)];
            fixatExt = [fixatExt; fixatExtremes(i,:)];
            if (i == (NfIni-1)) % si i+1 es la ultima y ha salido que no se junta con i, la copio tal cual
                fixatCen = [fixatCen; fixatCentroids(i+1,:)];
                fixatExt = [fixatExt; fixatExtremes(i+1,:)];
            end
            i=i+1;
        end
    end
end
%fprintf('We initially had %d fix., after merging we have %d fix. \n', NfIni, size(fixatCen,1));

return