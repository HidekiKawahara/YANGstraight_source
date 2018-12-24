function varargout = detailInspector(varargin)
% DETAILINSPECTOR MATLAB code for detailInspector.fig
%      DETAILINSPECTOR, by itself, creates a new DETAILINSPECTOR or raises the existing
%      singleton*.
%
%      H = DETAILINSPECTOR returns the handle to a new DETAILINSPECTOR or the handle to
%      the existing singleton*.
%
%      DETAILINSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETAILINSPECTOR.M with the given input arguments.
%
%      DETAILINSPECTOR('Property','Value',...) creates a new DETAILINSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detailInspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detailInspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detailInspector

% Last Modified by GUIDE v2.5 24-Dec-2018 19:05:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detailInspector_OpeningFcn, ...
                   'gui_OutputFcn',  @detailInspector_OutputFcn, ...
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

% --- Executes just before detailInspector is made visible.
function detailInspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to detailInspector (see VARARGIN)

% Choose default command line output for detailInspector
handles.output = hObject;
%handles
handles = setDefault(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes detailInspector wait for user response (see UIRESUME)
% uiwait(handles.detailInspectorGUI);
end


% --- Outputs from this function are returned to the command line.
function varargout = detailInspector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%--- private functions
function handles = setDefault(handles)
set(handles.freqAssignAxis, 'visible', 'off');
set(handles.snrAxis, 'visible', 'off');
set(handles.freqAxis, 'visible', 'off');
set(handles.waveAxis, 'visible', 'off');
set(handles.snrApplyButton, 'visible', 'off');
set(handles.freqApplyButton, 'visible', 'off');
set(handles.waveApplyButton, 'visible', 'off');
set(handles.playButton, 'enable', 'off');
set(handles.selectButton, 'enable', 'off');
set(handles.readButton, 'enable', 'on');
set(handles.reportGeneratorButton, 'enable', 'off');
set(handles.settingApplyButton, 'enable', 'off');
set(handles.quitButton, 'enable', 'on');
set(handles.stretchPopup, 'value', 2, 'enable', 'off');
set(handles.downsamplePopup, 'value', 1, 'enable', 'off');
set(handles.chpoctPopup, 'value', 3, 'enable', 'off');
handles.channels_per_octave = 12;
handles.downsampling = 1;
handles.stretching = 1.05;
end

%--- private functions

% --- Executes on button press in waveApplyButton.
function waveApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to waveApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlim = get(handles.waveAxis, 'xlim');
set(handles.snrAxis, 'xlim', xlim);
set(handles.freqAxis, 'xlim', xlim);
set(handles.playButton, 'enable', 'on');
set(handles.reportGeneratorButton, 'enable', 'on');
end

% --- Executes on button press in freqApplyButton.
function freqApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to freqApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlim = get(handles.freqAxis , 'xlim');
set(handles.snrAxis, 'xlim', xlim);
set(handles.waveAxis, 'xlim', xlim);
set(handles.playButton, 'enable', 'on');
set(handles.reportGeneratorButton, 'enable', 'on');
end

% --- Executes on button press in snrApplyButton.
function snrApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to snrApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlim = get(handles.snrAxis  , 'xlim');
set(handles.freqAxis, 'xlim', xlim);
set(handles.waveAxis, 'xlim', xlim);
set(handles.playButton, 'enable', 'on');
set(handles.reportGeneratorButton, 'enable', 'on');
end

% --- Executes on button press in readButton.
function readButton_Callback(hObject, eventdata, handles)
% hObject    handle to readButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[fileName, pathName, filterIndex] = uigetfile({'*.wav', '*.WAV'}, 'Select wave file:');
%disp(['fileName:' fileName]);
%disp(['pathName:' pathName]);
%disp(['filterIndex:' num2str(filterIndex)]);
set(handles.pathText, 'String', pathName);
set(handles.fileNameText, 'String', fileName);
[x, fs] = audioread([pathName fileName]);
if size(x, 2) > 1
    channelName = questdlg('Which channel to analyse?', ...
                         'Channel selection', ...
                         'Channel-1', 'Channel-2', 'Channel-2');
else
    channelName = 'Channel-1';
end
set(handles.channelText, 'String', channelName);
switch channelName
    case 'Channel-1'
        channelId = 1;
    case 'Channel-2'
        channelId = 2;
end
handles.sounddata = x(:, channelId);
x_sel = handles.sounddata + ...
    randn(length(handles.sounddata), 1) * std(handles.sounddata) * 2 ^ (-24);
handles.sounddata = x_sel;
handles.samplingFrequency = fs;
outputS = sourceInformationAnalysis(x_sel, fs, [1 length(x_sel)], 27.5, 1400, 6);
set(handles.freqAssignAxis, 'visible', 'on');
axes(handles.freqAssignAxis);
estSNR = 10 * log10(outputS.fixed_points_measure) - 18;
handles.fqSNRScatterHandle = semilogx(outputS.fixed_points_freq(:, 1:4), estSNR(:, 1:4), '.');
grid on;
xlabel('frequency (Hz)');
ylabel('estimated SNR (dB)');
set(handles.selectButton, 'enable', 'on');
set(handles.freqAxis, 'visible', 'off');
set(handles.snrAxis, 'visible', 'off');
set(handles.waveAxis, 'visible', 'off');
guidata(hObject, handles);
end

% --- Executes on button press in selectButton.
function selectButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
xlim = get(handles.freqAssignAxis, 'xlim');
handles.frequency_low = xlim(1);
handles.frequency_high = xlim(2);
set(handles.freqAssignAxis, 'visible', 'off');
for ii = 1:4
set(handles.fqSNRScatterHandle(ii), 'visible', 'off');
end
x_sel = handles.sounddata;
fs = handles.samplingFrequency;
%handles.channels_per_octave = 12;
%handles.downsampling = 1;
%handles.stretching = 1.05;
outputS = sourceInformationAnalysis(x_sel, fs, [1 length(x_sel)], ...
    handles.frequency_low, handles.frequency_high, handles.channels_per_octave, ...
    handles.downsampling, handles.stretching);
estSNR = 10 * log10(outputS.fixed_points_measure) - 18;
set(handles.freqAxis, 'visible', 'on');
set(handles.snrAxis, 'visible', 'on');
set(handles.waveAxis, 'visible', 'on');
axes(handles.freqAxis);
handles.currentAnalysisResult = outputS;
tx = outputS.time_axis_wavelet;
tx_audio = (1:length(x_sel)) / fs;
semilogy(tx, outputS.fixed_points_freq(:, 4:-1:1), '.');
grid on;
axis([tx([1 end]) handles.frequency_low handles.frequency_high]);
axes(handles.snrAxis);
plot(tx, estSNR(:, 4:-1:1), '.');
grid on;
axis([tx([1 end]) -15 50]) 
axes(handles.waveAxis);
plot(tx_audio, x_sel);
set(handles.selectButton, 'enable', 'off');
set(handles.snrApplyButton, 'visible', 'on');
set(handles.freqApplyButton, 'visible', 'on');
set(handles.waveApplyButton, 'visible', 'on');
set(handles.stretchPopup, 'enable', 'on');
set(handles.downsamplePopup, 'enable', 'on');
set(handles.chpoctPopup, 'enable', 'on');
set(gca, 'xlim', tx_audio([1 end]));
guidata(hObject, handles);
end

% --- Executes on button press in playButton.
function playButton_Callback(hObject, eventdata, handles)
% hObject    handle to playButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
xlim = get(handles.waveAxis, 'xlim');
x_sel = handles.sounddata;
fs = handles.samplingFrequency;
handles.playerObject = audioplayer(x_sel(round(xlim(1) * fs + 1):round(xlim(2) * fs)), fs);
playblocking(handles.playerObject);
end


% --- Executes on selection change in chpoctPopup.
function chpoctPopup_Callback(hObject, eventdata, handles)
% hObject    handle to chpoctPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chpoctPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chpoctPopup
handles = guidata(hObject);
contents = cellstr(get(hObject,'String'));
handles.channels_per_octave = str2double(contents{get(hObject,'Value')});
set(handles.settingApplyButton, 'enable', 'on');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function chpoctPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chpoctPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in downsamplePopup.
function downsamplePopup_Callback(hObject, eventdata, handles)
% hObject    handle to downsamplePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns downsamplePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from downsamplePopup
handles = guidata(hObject);
contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'on'
        handles.downsampling = 1;
    case 'off'
        handles.downsampling = 0;
end
set(handles.settingApplyButton, 'enable', 'on');
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function downsamplePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downsamplePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in stretchPopup.
function stretchPopup_Callback(hObject, eventdata, handles)
% hObject    handle to stretchPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stretchPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stretchPopup
handles = guidata(hObject);
contents = cellstr(get(hObject,'String'));
handles.stretching = str2double(contents{get(hObject,'Value')});
guidata(hObject, handles);
set(handles.settingApplyButton, 'enable', 'on');
end

% --- Executes during object creation, after setting all properties.
function stretchPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stretchPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in settingApplyButton.
function settingApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to settingApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
currentXlim = get(handles.waveAxis, 'xlim');
x_sel = handles.sounddata;
fs = handles.samplingFrequency;
%handles.channels_per_octave = 12;
%handles.downsampling = 1;
%handles.stretching = 1.05;
outputS = sourceInformationAnalysis(x_sel, fs, [1 length(x_sel)], ...
    handles.frequency_low, handles.frequency_high, handles.channels_per_octave, ...
    handles.downsampling, handles.stretching);
estSNR = 10 * log10(outputS.fixed_points_measure) - 18;
axes(handles.freqAxis);
tx = outputS.time_axis_wavelet;
%tx_audio = (1:length(x_sel)) / fs;
handles.currentAnalysisResult = outputS;
semilogy(tx, outputS.fixed_points_freq(:, 4:-1:1), '.');
grid on;
axis([currentXlim handles.frequency_low handles.frequency_high]);
axes(handles.snrAxis);
plot(tx, estSNR(:, 4:-1:1), '.');
grid on;
axis([currentXlim -15 50]) 
%axes(handles.waveAxis);
%plot(tx_audio, x_sel);
guidata(hObject, handles);
set(handles.settingApplyButton, 'enable', 'off');
set(handles.reportGeneratorButton, 'enable', 'on');
end

% --- Executes on button press in reportGeneratorButton.
function reportGeneratorButton_Callback(hObject, eventdata, handles)
% hObject    handle to reportGeneratorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
currentXlim = get(handles.waveAxis, 'xlim');
x_sel = handles.sounddata;
fs = handles.samplingFrequency;
pathName = get(handles.pathText, 'String');
fileName = get(handles.fileNameText, 'String');
channelName = get(handles.channelText, 'String');
figureHandle = figure;
%handles.currentAnalysisResult = outputS;
outputS = handles.currentAnalysisResult;
estSNR = 10 * log10(outputS.fixed_points_measure) - 18;
tx = outputS.time_axis_wavelet;
tx_audio = (1:length(x_sel)) / fs;
param_selector = tx > currentXlim(1) & tx < currentXlim(2);
wave_selector = tx_audio > currentXlim(1) & tx_audio < currentXlim(2);
subplot(312)
semilogy(tx(param_selector), outputS.fixed_points_freq(param_selector, 4:-1:1), '.');
ylabel('frequency (Hz)');
grid on;
axis([currentXlim handles.frequency_low handles.frequency_high]);
subplot(313)
plot(tx(param_selector), estSNR(param_selector, 4:-1:1), '.');
ylabel('est. SNR (dB)');
xlabel('time (s)');
grid on;
axis([currentXlim -15 50]) 
subplot(311);
plot(tx_audio(wave_selector), x_sel(wave_selector));
grid on;
axis([currentXlim max(abs(x_sel)) * [-1 1]]);
title(['file: ' fileName '  '  channelName], 'interpreter', 'none');
figure(figureHandle)
outfileNameRoot = ['rep' datestr(now, 30)];
outfileNameFig = [outfileNameRoot 'r.fig'];
%print('-depsc', outfileName);
savefig(outfileNameFig);
outfileNameTxt = [outfileNameRoot 'r.txt'];
fid = fopen(outfileNameTxt,'w');
fprintf(fid, 'Pathname: %s \n', pathName);
fprintf(fid, 'FileName: %s \n', fileName);
fclose(fid);
end


% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.detailInspectorGUI);
end


% --- Executes on button press in fullSoomOutButton.
function fullSoomOutButton_Callback(hObject, eventdata, handles)
% hObject    handle to fullSoomOutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
x_sel = handles.sounddata;
fs = handles.samplingFrequency;
xlim = [0 length(x_sel) / fs];
set(handles.freqAxis , 'xlim', xlim, 'ylim', [handles.frequency_low handles.frequency_high]);
set(handles.snrAxis, 'xlim', xlim, 'ylim', [-15 50]);
set(handles.waveAxis, 'xlim', xlim, 'ylim', max(abs(x_sel)) * [-1 1]);
set(handles.playButton, 'enable', 'off');
end
