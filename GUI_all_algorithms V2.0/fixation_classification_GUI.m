
function varargout = fixation_classification_GUI(varargin)
% FIXATION_CLASSIFICATION_GUI MATLAB code for fixation_classification_GUI.fig
%      FIXATION_CLASSIFICATION_GUI, by itself, creates a new FIXATION_CLASSIFICATION_GUI or raises the existing
%      singleton*.
%
%      H = FIXATION_CLASSIFICATION_GUI returns the handle to a new FIXATION_CLASSIFICATION_GUI or the handle to
%      the existing singleton*.
%
%      FIXATION_CLASSIFICATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIXATION_CLASSIFICATION_GUI.M with the given input arguments.
%
%      FIXATION_CLASSIFICATION_GUI('Property','Value',...) creates a new FIXATION_CLASSIFICATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fixation_classification_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fixation_classification_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLESg
% Edit the above text to modify the response to help fixation_classification_GUI

% Last Modified by GUIDE v2.5 18-Oct-2018 11:26:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fixation_classification_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @fixation_classification_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before fixation_classification_GUI is made visible.
function fixation_classification_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fixation_classification_GUI (see VARARGIN)

% Choose default command line output for fixation_classification_GUI
handles.output = hObject;

handles.filename = '';
handles.path = '';
handles.stimuli = {};
handles.stimulusIdx = 0;
handles.validity = [];
handles.fixPoints = [];
handles.nRounds = 0;
handles.nStimuli = 0;
axes(handles.gazeAxes);
axis equal;

addpath(genpath('Classification_algorithm'));
addpath('ClassificationResults');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fixation_classification_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = fixation_classification_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.txt');
handles.filename = file;
handles.path = path;

set(handles.stateText, 'String', 'Loading...');
set(handles.stateText, 'ForegroundColor', [0.7 0.7 0.7]);
% Get the clasificated stimuli
[classificated_stimuli, nRounds, nStimuli, error] = classificate_fixations(path, file);

if(error ~= 0)
    set(handles.stateText, 'String', 'Error');
    set(handles.stateText, 'ForegroundColor', [1 0 0]);
    handles.stimulusIdx = 0;
else
    set(handles.stateText, 'String', 'Loaded');
    set(handles.stateText, 'ForegroundColor', [0 1 0]);
    handles.stimuli = classificated_stimuli;
    handles.stimulusIdx = 1;
    handles.nRounds = nRounds;
    handles.nStimuli = nStimuli;
    
    
    % Display the first stimuli
    displayInfo(handles)
    plotVelocities(handles);
    plotGazePosition(handles);
    plotFixationAxes(handles);
    plotAlgorithms(handles);
end

% Update handles structure
guidata(hObject, handles);
end


% --- Executes on button press in validButton.
function validButton_Callback(hObject, eventdata, handles)
% hObject    handle to validButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    handles.validity(handles.stimulusIdx) = 1;
    fprintf(handles.fileID, [num2str(handles.stimuli(handles.stimulusIdx).Algorithm) ', ' ...
                        num2str(handles.stimuli(handles.stimulusIdx).Round) ', '...
                        num2str(handles.stimuli(handles.stimulusIdx).Stimulus) ', 1' char(10)]);
    
    %Display next stimuli
    if(handles.stimulusIdx < size(handles.stimuli,2))
        handles.stimulusIdx = handles.stimulusIdx + 1;
        plotVelocities(handles);
        plotGazePosition(handles)
    else
        %Not more stimuli
        errordlg('Not more stimuli')
    end
    
    
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on button press in invalidButton.
function invalidButton_Callback(hObject, eventdata, handles)
% hObject    handle to invalidButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    handles.validity(handles.stimulusIdx) = 0;
    fprintf(handles.fileID, [num2str(handles.stimuli(handles.stimulusIdx).Algorithm) ', ' ...
                        num2str(handles.stimuli(handles.stimulusIdx).Round) ', '...
                        num2str(handles.stimuli(handles.stimulusIdx).Stimulus) ', 0' char(10)]);
                    
    %Display next stimuli
    if(handles.stimulusIdx < size(handles.stimuli,2))
        handles.stimulusIdx = handles.stimulusIdx + 1;
        plotVelocities(handles);
        plotGazePosition(handles);
        plotFixationAxes(handles);
        plotAlgorithms(handles);
    else
        %Not more stimuli
        errordlg('Not more stimuli')
    end

    % Update handles structure
    guidata(hObject, handles);
end

function plotDistanceFromTheStimuli(handles)

    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    gazePos = handles.stimuli(handles.stimulusIdx).GazePos;
    stimPos = repmat(handles.stimuli(handles.stimulusIdx).StimPos, [size(gazePos, 1), 1]);
    dif = stimPos - gazePos;
    distance = sqrt(sum((dif(:,1) + dif(:,2)) .^ 2, 2));
    
    % Divide into points classificated as fixation and others
    
    fixations = handles.stimuli(handles.stimulusIdx).Classification == 1;    
    fixIdx = [1; find(diff(fixations)); size(fixations,1)];
    
    axes(handles.distanceAxes);
    
    for(i=1:size(fixIdx)-1)
        fixation = fixations(fixIdx(i)+1);
        if(fixation)
            plot(handles.distanceAxes, handles.stimuli(handles.stimulusIdx).Timestamp(fixIdx(i):fixIdx(i+1)), distance(fixIdx(i):fixIdx(i+1)), 'r');
        else
            plot(handles.distanceAxes, handles.stimuli(handles.stimulusIdx).Timestamp(fixIdx(i):fixIdx(i+1)), distance(fixIdx(i):fixIdx(i+1)), 'b');
        end
        hold on;
    end
        
    xlabel(handles.distanceAxes, ['Time (ms)']);
    ylabel(handles.distanceAxes, 'Distance (mm)');
    handles.distanceAxes.XMinorGrid = 'on';
    handles.distanceAxes.XGrid = 'on';
    hold off;
end

function plotVelocities(handles)
    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    
    index = handles.stimulusIdx + handles.nRounds*handles.nStimuli*(handles.popupmenu1.Value-1);
    
    axes(handles.velocityAxes);
    %plot(handles.velocityAxes, handles.stimuli(index).Timestamp, handles.stimuli(index).Velocities, 'b');
        
    velocity = handles.stimuli(index).Velocities;
    % Divide into points classificated as fixation and others
    fixIdx = [1; find(diff(handles.stimuli(index).Classification))+1; size(velocity,1)];
    
    for(i=1:length(fixIdx)-1)
        if(handles.stimuli(index).Classification(fixIdx(i))==1)
            plot(handles.velocityAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), velocity(fixIdx(i):fixIdx(i+1)), 'r');
        elseif(handles.stimuli(index).Classification(fixIdx(i))==2)
            plot(handles.velocityAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), velocity(fixIdx(i):fixIdx(i+1)), 'g');
        else
            plot(handles.velocityAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), velocity(fixIdx(i):fixIdx(i+1)), 'b');
        end
        hold(handles.velocityAxes, 'on');
    end
    
    hold(handles.velocityAxes, 'on');
    xlabel(handles.velocityAxes, ['Time (ms)']);
    ylabel(handles.velocityAxes, 'Velocity (mm/s)');    
    
    title(handles.velocityAxes, 'Filtered velocity');
    xlim(handles.velocityAxes, [handles.stimuli(index).Timestamp(1), ...
        handles.stimuli(index).Timestamp(end)]);
    handles.velocityAxes.XMinorGrid = 'on';
    handles.velocityAxes.XGrid = 'on';
    hold(handles.velocityAxes, 'off');
end


function plotGazePosition(handles, markedPointIdx)
    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end    
    
    axes(handles.gazeAxes);    

    %plot stimuli
    plot(handles.gazeAxes, handles.stimuli(handles.stimulusIdx).StimPos(1),...
        handles.stimuli(handles.stimulusIdx).StimPos(2), 'o', ...
        'MarkerSize', 20, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    hold(handles.gazeAxes, 'on');
    %plot gaze position
    gazePos = handles.stimuli(handles.stimulusIdx).GazePos;
    for i = 1:(length(gazePos)-1)
        if(gazePos(i,1)~=0 && gazePos(i+1,1)~=0)
            color =stimulus2Color(i,length(gazePos));
            plot(handles.gazeAxes, gazePos(i:i+1,1), gazePos(i:i+1,2),'-', 'MarkerEdgeColor', 'none', 'Color', color);
        else
            %color = [0.7, 0.7, 0.7];
            %plot(handles.gazeAxes, gazePos(i:i+1,1), gazePos(i:i+1,2),'-', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'Color', color);
        end
    end
    %plot(handles.gazeAxes, handles.stimuli(handles.stimulusIdx).GazePos(:,1), handles.stimuli(handles.stimulusIdx).GazePos(:,2));
 
    if(nargin==2)
        
        frameGazePos = handles.stimuli(handles.stimulusIdx).GazePos(markedPointIdx,:);
        plot(handles.gazeAxes, frameGazePos(1), frameGazePos(2), 'k.', 'MarkerSize',10);
    end
    
    xlabel(handles.gazeAxes, ['mm']);
    ylabel(handles.gazeAxes, 'mm');
    axis(handles.gazeAxes, [0 254 0 169.33]);    
    
    hold(handles.gazeAxes, 'off');
end

function plotFixationAxes(handles)

    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    index = handles.stimulusIdx + handles.nRounds*handles.nStimuli*(handles.popupmenu1.Value-1);
    
    gazePos = handles.stimuli(index).GazePos;
    stimPos = repmat(handles.stimuli(index).StimPos, [size(gazePos, 1), 1]);
    dif = stimPos - gazePos;
    distance = sqrt(dif(:,1).^ 2 + dif(:,2) .^ 2);
    % Divide into points classificated as fixation and others
    fixIdx = [1; find(diff(handles.stimuli(index).Classification))+1; size(gazePos,1)];
    
    axes(handles.fixationAxes);
    
    for(i=1:length(fixIdx)-1)
        if(handles.stimuli(index).Classification(fixIdx(i))==1)
            plot(handles.fixationAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), distance(fixIdx(i):fixIdx(i+1)), 'r');
        elseif(handles.stimuli(index).Classification(fixIdx(i))==2)
            plot(handles.fixationAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), distance(fixIdx(i):fixIdx(i+1)), 'g');
        else
            plot(handles.fixationAxes, handles.stimuli(index).Timestamp(fixIdx(i):fixIdx(i+1)), distance(fixIdx(i):fixIdx(i+1)), 'b');
        end
        hold(handles.fixationAxes, 'on');
    end
    
    xlabel(handles.fixationAxes, ['Time (ms)']);
    ylabel(handles.fixationAxes, 'Distance (mm)');
    xlim(handles.fixationAxes, [handles.stimuli(index).Timestamp(1), ...
        handles.stimuli(handles.stimulusIdx).Timestamp(end)]);
    title(handles.fixationAxes, 'Filtered Distance');
    
    handles.fixationAxes.XGrid = 'on';
    handles.fixationAxes.XMinorGrid = 'on';
    hold(handles.fixationAxes, 'off');
end


function plotAlgorithms(handles)
    if(handles.stimulusIdx == 0)
        errordlg('Stimuli not loaded')
        return;
    end
    actualStimulusRound = handles.stimuli(handles.stimulusIdx).Round;
    actualStimulusId = handles.stimuli(handles.stimulusIdx).Stimulus;
    axes(handles.algorithmsAxes);
    
    graphicColumn = 1;
    for i=1:length(handles.stimuli)
        if(handles.stimuli(i).Round == actualStimulusRound && ...
                handles.stimuli(i).Stimulus == actualStimulusId)
          

            fixIdx = [1; find(diff(handles.stimuli(i).Classification))+1; size(handles.stimuli(i).Classification,1)];
            
            for j=1:size(fixIdx)-1
                row = ones(size(handles.stimuli(i).Timestamp(fixIdx(j):fixIdx(j+1))))*graphicColumn;
                if(handles.stimuli(i).Classification(fixIdx(j)) == 1)
                    plot(handles.algorithmsAxes, handles.stimuli(i).Timestamp(fixIdx(j):fixIdx(j+1)), row, 'r', 'LineWidth',5);
                elseif (handles.stimuli(i).Classification(fixIdx(j)) == 2)
                    plot(handles.algorithmsAxes, handles.stimuli(i).Timestamp(fixIdx(j):fixIdx(j+1)), row, 'g', 'LineWidth',5);
                elseif (handles.stimuli(i).Classification(fixIdx(j)) == 3)
                    plot(handles.algorithmsAxes, handles.stimuli(i).Timestamp(fixIdx(j):fixIdx(j+1)), row, 'b', 'LineWidth',5);
                else
                    plot(handles.algorithmsAxes, handles.stimuli(i).Timestamp(fixIdx(j):fixIdx(j+1)), row, 'b', 'LineWidth',5);
                end                
                hold(handles.algorithmsAxes, 'on');
            end
            graphicColumn = graphicColumn + 1;
        end
    end
    xlabel(handles.algorithmsAxes, ['Time (ms)']);
    ylim(handles.algorithmsAxes, [0, 11]);
    xlim(handles.algorithmsAxes, [handles.stimuli(handles.stimulusIdx).Timestamp(1), ...
    handles.stimuli(handles.stimulusIdx).Timestamp(end)]);
    handles.algorithmsAxes.XGrid = 'on';
    handles.algorithmsAxes.XMinorGrid = 'on';
    hold(handles.algorithmsAxes, 'off');
end

% --- Executes on button press in addPoint.
function addPoint_Callback(hObject, eventdata, handles)
% hObject    handle to addPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
fixTimestamp = get(0,'userdata');
fixIdx = find(handles.stimuli(handles.stimulusIdx).Timestamp == fixTimestamp);
handles.fixPoints = sort([handles.fixPoints; fixIdx]);


guidata(hObject, handles);
plotFixationAxes(handles);
end


% --- Executes on button press in deleteLastPoint.
function deleteLastPoint_Callback(hObject, eventdata, handles)
% hObject    handle to deleteLastPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if(size(handles.fixPoints, 2) == 0)
        errordlg('Fixation points list is empty')
        return;
    end
handles.fixPoints(end) = [];

guidata(hObject, handles);
plotFixationAxes(handles);
end


% --- Executes on button press in selectFixationsMode.
function selectFixationsMode_Callback(hObject, eventdata, handles)
% hObject    handle to selectFixationsMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.dcm_obj = datacursormode();
    datacursormode on;
    set(handles.dcm_obj,'UpdateFcn',{@myupdatefcn});

    guidata(hObject, handles);
end

function txt = myupdatefcn(~, event_obj)

     % get the handles structure for the figure/GUI
    handles = guidata( event_obj.Target );
    hAxesParent  = get(get(event_obj,'Target'),'Parent');
    pos = event_obj.Position;    

    drawnow nocallbacks;
    if(hAxesParent == handles.fixationAxes) 
        idx = find(handles.stimuli(handles.stimulusIdx).Timestamp ==  pos(1));
        plotGazePosition(handles, idx);
        plotVelocities(handles, idx); 
        %clear(handles.fixationAxes);
        %plotFixationAxes(handles);
        
        txt = {['Timestamp: ',num2str(pos(1)), ' ms'],...
            ['Disance: ',num2str(pos(2)), ' mm']};
    elseif(hAxesParent == handles.velocityAxes)        
        %{
        idx = find(handles.stimuli(handles.stimulusIdx).Timestamp ==  pos(1));
        plotGazePosition(handles, idx);
        plotFixationAxes(handles, idx);
        clear(handles.velocityAxes);
        plotVelocities(handles);
        %}
        
        
        txt = {['Timestamp: ',num2str(pos(1)), ' ms'],...
            ['Velocity: ',num2str(pos(2)), ' mm/s']};
    elseif(hAxesParent == handles.gazeAxes)        
        %{
        idx = find(ismember(handles.stimuli(handles.stimulusIdx).GazePos, pos, 'rows'));
        plotVelocities(handles, idx);
        plotFixationAxes(handles, idx);
        %clear(handles.gazeAxes);
        %plotGazePosition(handles);
        %}
        
        txt = {['X: ',num2str(pos(1)), ' mm'],...
            ['Y: ',num2str(pos(2)), ' mm']};
    end
    displayInfo(handles);
    set(0,'userdata',pos(1));
end


% --- Executes on button press in nextStimulusButton.
function nextStimulusButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextStimulusButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.stimuli(handles.stimulusIdx).definedFixations = handles.fixPoints;
    handles.fixPoints = [];
    
    if(handles.stimulusIdx < handles.nRounds*handles.nStimuli)
        handles.stimulusIdx = handles.stimulusIdx + 1;
        displayInfo(handles);
        plotVelocities(handles);
        plotGazePosition(handles);
        plotFixationAxes(handles);
        plotAlgorithms(handles);
    else
        %Not more stimuli
        errordlg('No more stimuli');
    end
    
guidata(hObject, handles);
end


% --- Executes on button press in preStimulusButton.
function preStimulusButton_Callback(hObject, eventdata, handles)
% hObject    handle to preStimulusButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if(handles.stimulusIdx > 1)        
        handles.stimulusIdx = handles.stimulusIdx - 1;        
        handles.provisionalPoint = [];
        handles.fixPoint = [];
        
        displayInfo(handles);
        plotVelocities(handles);
        plotGazePosition(handles);
        plotFixationAxes(handles);
        plotAlgorithms(handles);
    else
        %Not more stimuli
        errordlg('No more stimuli');
    end
    
guidata(hObject, handles);
end

% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    idx = 1;
    stimuli.Round = [];
    stimuli.Stimulus= [];
    stimuli.RawGazePosX = [];
    stimuli.RawGazePosY = [];
    stimuli.Timestamp = [];
    stimuli.Classification = [];
    algorithm = handles.stimuli(idx).Algorithm;
    while(idx<=size(handles.stimuli, 2))
        nSamples = length(handles.stimuli(idx).Timestamp);
        stimuli.Round = [stimuli.Round; ones(nSamples,1)*handles.stimuli(idx).Round];
        stimuli.Stimulus= [stimuli.Stimulus; ones(nSamples,1)*handles.stimuli(idx).Stimulus];
        stimuli.RawGazePosX = [stimuli.RawGazePosX; handles.stimuli(idx).RawGazePos(:,1)];
        stimuli.RawGazePosY = [stimuli.RawGazePosY; handles.stimuli(idx).RawGazePos(:,2)];
        stimuli.Timestamp = [stimuli.Timestamp; handles.stimuli(idx).Timestamp];
        stimuli.Classification = [stimuli.Classification; handles.stimuli(idx).Classification];
        idx = idx +1;
        
        % Stimuli of the current algorithm ended
        if(idx<=size(handles.stimuli, 2) && strcmp(algorithm, handles.stimuli(idx).Algorithm) == 0)             
            mkdir('ClassificationResults', algorithm);
            
            name = strsplit(handles.filename, '.');
            filename = char(strcat(cd, '\ClassificationResults\', algorithm, '\', name(1), '_classificated.mat'));
            save(filename, '-struct', 'stimuli');
                        
            algorithm = handles.stimuli(idx).Algorithm;
            stimuli.Round = [];
            stimuli.Stimulus= [];
            stimuli.RawGazePosX = [];
            stimuli.RawGazePosY = [];
            stimuli.Timestamp = [];
            stimuli.Classification = [];
        end
    end
end

function displayInfo(handles)
text = ['NHC: ' num2str(handles.stimuli(handles.stimulusIdx).NHC) char(10)...
    'Round: ' num2str(handles.stimuli(handles.stimulusIdx).Round) char(10)...
    'Stimulus: ' num2str(handles.stimuli(handles.stimulusIdx).Stimulus)];                   
                    
set(handles.stimuliText, 'String', text);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotVelocities(handles);
plotGazePosition(handles);
plotFixationAxes(handles);
plotAlgorithms(handles);
        
end
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
