function write2csv(ETparams, outname)
% Writes output from ET params to csv-file
% output is x, y, type

nTrials = 1;
for i = 1:nTrials
    nRows = length(ETparams.data(1, i).X);
    if nRows < 1
        continue
    end
    D_temp = zeros(nRows, 5);
    D_temp(:, 1) = ETparams.data(1, i).T;
    D_temp(:, 2) = ETparams.data(1, i).X;
    D_temp(:, 3) = ETparams.data(1, i).Y;   
    D_temp(logical(ETparams.fixationIdx(1, i).Idx), 4) = 1;
    D_temp(logical(ETparams.saccadeIdx(1, i).Idx), 4) =  2;
    D_temp(logical(ETparams.glissadeIdx(1, i).Idx), 4) = 3;
    D_temp(logical(ETparams.nanIdx(1, i).Idx), 4) = 4;
    D_temp(:, 5) = ETparams.data(1, i).vel;   

end

% Write to csv-file
cHeader = {'t' 'x' 'y' 'event_type','vel'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
%write header to file
fid = fopen(outname,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
%write data to end of file
dlmwrite(outname,D_temp,'-append');

