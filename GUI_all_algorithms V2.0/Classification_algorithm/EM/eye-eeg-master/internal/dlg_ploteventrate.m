% dlg_ploteventrate - pops user dialogue, called by pop_ploteventrate.m
%                see >> help pop_ploteventrate
%
% Copyright (C) 2009-2016 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [rate_event] = dlg_ploteventrate(callingFcn,EEG)

geometry = { 1 [2 1.3 0.4]};

%% callbacks
cbevent = ['if ~isfield(EEG.event, ''type'')' ...
    '   errordlg2(''No type field'');' ...
    'else' ...
    '   if isnumeric(EEG.event(1).type),' ...
    '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
    '   else,' ...
    '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
    '   end;' ...
    '   if ~isempty(tmps)' ...
    '       set(findobj(''parent'', gcbf, ''tag'', ''rate_event''), ''string'', tmpstr);' ...
    '   end;' ...
    'end;' ...
    'clear tmps tmpv tmpstr tmpfieldnames;' ];

%% main menu
uilist = {...
    {'style', 'text', 'string', 'Plot rate of which event?', 'fontweight','bold'},...
    ...
    {'style', 'text', 'string', 'Event type:'},...
    {'style', 'edit', 'string', 'saccade', 'tag', 'rate_event' },...
    {'style', 'pushbutton', 'string', '...', 'callback', cbevent}, ...
    };


%% make GUI
[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Plot event rate -- ', callingFcn]);

%% process user input (cancel)
if isempty(results)
    return
end

%% collect dialogue inputs
rate_event = outstruct.rate_event;

end