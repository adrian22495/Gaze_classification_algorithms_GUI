clear all, close all, clc
load('1250Hz_3_Participants.mat')


sz = size(ETdata);

for part = 1:sz(1)
    for trial = 1:sz(2)
        disp([part, trial])
        x= ETdata(part,trial).X;
        y= ETdata(part,trial).Y;
        
        
        cHeader = {'x' 'y'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); %cHeader in text with commas
        %write header to file
        outname = strcat(string(part),'_', string(trial),'.csv');
        fid = fopen(outname,'w'); 
        fprintf(fid,'%s\n',textHeader);
        fclose(fid);
        %write data to end of file
        dlmwrite(outname,[x', y'],'-append');        
    end
end




