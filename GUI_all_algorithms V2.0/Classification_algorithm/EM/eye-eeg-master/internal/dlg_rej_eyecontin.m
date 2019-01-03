% dlg_rej_eyecontin - pops dialogue, called by pop_rej_eyecontin
%                see >> help pop_rej_eyecontin
%
% Copyright (C) 2009-2018 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de 

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

function [chns minvals maxvals windowsize rejectionmode] = dlg_rej_eyecontin(callingFcn,EEG)

geometry(1:12) = {1 1 1 [0.5 0.4 0.2 0.4 0.4] [0.5 0.4 0.2 0.4 0.4] [0.5 0.4 0.2 0.4 0.4] [0.5 0.4 0.2 0.4 0.4] [0.5 0.4 0.2 0.4 0.4] 1 [0.5 0.4 1.0] 1 [0.5 1.0 0.4]};

rejectionmethods = {'1. remove bad intervals from data (cut data)'; sprintf('2. add \"bad_ET\" markers to EEG.event (keep data)')};

uilist(1:36) = {...
    {'Style', 'text', 'string', 'Detect and/or reject intervals with out-of-range eye track','fontweight', 'bold'},...
    {'Style', 'text', 'string', 'Note: eye blinks and missing data are usually characterized by zeros (values < 1)'},...
    {},...
    {'Style', 'text', 'string', ' '},...
    {'Style', 'text', 'string', 'ET channel index'},...
    {'Style', 'text', 'string', ''},...
    {'Style', 'text', 'string', 'Reject values <'},...
    {'Style', 'text', 'string', 'Reject values >'},...
    ...
    {'Style', 'text', 'string', 'Left eye, horiz. (X).:'},...
    {'Style', 'edit', 'string', '', 'tag', 'chans_LX' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(EEG.chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chans_LX''), ''string'',tmp); clear tmp tmpchanlocs tmpval'},...
    {'Style', 'edit', 'string', '1' },...
    {'Style', 'edit', 'string', '1024' },...
    ...
    {'Style', 'text', 'string', 'Left eye, verti. (Y):'},...
    {'Style', 'edit', 'string', '', 'tag', 'chans_LY' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(EEG.chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chans_LY''), ''string'',tmp); clear tmp tmpchanlocs tmpval' },...
    {'Style', 'edit', 'string', '1' },...
    {'Style', 'edit', 'string', '768' },...
    ...
    {'Style', 'text', 'string', 'Right eye, horiz. (X):'},...
    {'Style', 'edit', 'string', '', 'tag', 'chans_RX' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(EEG.chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chans_RX''), ''string'',tmp); clear tmp tmpchanlocs tmpval' },...
    {'Style', 'edit', 'string', '1' },...
    {'Style', 'edit', 'string', '1024' },...
    ...
    {'Style', 'text', 'string', 'Right eye, verti. (Y):'},...
    {'Style', 'edit', 'string', '', 'tag', 'chans_RY' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(EEG.chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chans_RY''), ''string'',tmp); clear tmp tmpchanlocs tmpval' },...
    {'Style', 'edit', 'string', '1' },...
    {'Style', 'edit', 'string', '768' },...
    {},...
    {'Style', 'text', 'string', 'Remove an extra '},...
    {'Style', 'edit', 'string', '100', 'tag', 'windowsize' }, ...
    {'Style', 'text', 'string', 'ms before and after out-of-range intervals'},...
    {},...
    {'Style', 'text', 'string', 'What to do with bad ET intervals? '},...
    {'Style', 'popupmenu', 'string',rejectionmethods,'tag','rejectionmode','value',2},...
    {'Style', 'text', 'string', ''}
    };
    % 
%%
[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Reject continuous intervals based on eye tracking data -- ', callingFcn]);

if isempty(results)
    return
end

%% prepare output
[chn1 chn1txt] = eeg_decodechan(EEG.chanlocs,results{1});
[chn2 chn2txt] = eeg_decodechan(EEG.chanlocs,results{4});
[chn3 chn3txt] = eeg_decodechan(EEG.chanlocs,results{7});
[chn4 chn4txt] = eeg_decodechan(EEG.chanlocs,results{10});

min1 = str2num(results{2});
min2 = str2num(results{5});
min3 = str2num(results{8});
min4 = str2num(results{11});

max1 = str2num(results{3});
max2 = str2num(results{6});
max3 = str2num(results{9});
max4 = str2num(results{12});

chns = []; minvals = []; maxvals = [];

% use inputs of fields were channel index was provided (max: 4 inputs)
if ~isempty(chn1) 
    if isempty(min1), min1 = -inf; end
    if isempty(max1), max1 =  inf; end
    chns = [chns chn1]; minvals = [minvals min1]; maxvals = [maxvals max1];
end
if ~isempty(chn2) 
    if isempty(min2), min2 = -inf; end
    if isempty(max2), max2 =  inf; end
    chns = [chns chn2]; minvals = [minvals min2]; maxvals = [maxvals max2];
end
if ~isempty(chn3) 
    if isempty(min3), min3 = -inf; end
    if isempty(max3), max3 =  inf; end
    chns = [chns chn3]; minvals = [minvals min3]; maxvals = [maxvals max3];
end
if ~isempty(chn4) 
    if isempty(min4), min4 = -inf; end
    if isempty(max4), max4 =  inf; end
    chns = [chns chn4]; minvals = [minvals min4]; maxvals = [maxvals max4];
end
windowsize = round(str2num(results{13})./(1000/EEG.srate));
rejectionmode  = outstruct.rejectionmode;

end