function [timestamp,lx,ly,rx,ry] = importDIVEdata( datapath, missingx,missingy)
%IMPORTDIVEDATA Summary of this function goes here
%   Detailed explanation goes here
    columnData = getColumnData(datapath);
    Tdata_idx = 8;
    LXdata_idx = 9;
    LYdata_idx = 10;
    RXdata_idx = 21;
    RYdata_idx = 22;
    
    timestamp = columnData{1,Tdata_idx};
    lx = columnData{1,LXdata_idx};
    ly = columnData{1,LYdata_idx};
    rx = columnData{1,RXdata_idx};
    ry = columnData{1,RYdata_idx};
    
    lx(lx==0) = missingx;
    ly(ly==0) = missingy;
    rx(rx==0) = missingx;
    ry(ry==0) = missingy;
end

function columnData = getColumnData(datapath)
    fileID = fopen([datapath]);
    
    header1  = fgetl(fileID);
    Ch1 = textscan(header1,'%s', 'Delimiter',',');
    header1Data = fgetl(fileID);
    columnHeader1Data = textscan(header1Data, '%s %d %s %d %d %d %d %d %s %s %f %f %f %d %d %f %s %s %d', 'Delimiter',','); % Eye tracker used at the end

    header2 = fgetl(fileID);
    columnHeader2 = textscan(header2, '%s', 'Delimiter',',');

    % Read and store data
    columnData = textscan(fileID, '%s %d %f %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %f', 'Delimiter',',');
end