function varargout = Cal_node_detect(varargin)
% CAL_NODE_DETECT MATLAB code for Cal_node_detect.fig
%      CAL_NODE_DETECT, by itself, creates a new CAL_NODE_DETECT or raises the existing
%      singleton*.
%
%      H = CAL_NODE_DETECT returns the handle to a new CAL_NODE_DETECT or the handle to
%      the existing singleton*.
%
%      CAL_NODE_DETECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAL_NODE_DETECT.M with the given input arguments.
%
%      CAL_NODE_DETECT('Property','Value',...) creates a new CAL_NODE_DETECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Cal_node_detect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Cal_node_detect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Cal_node_detect

% Last Modified by GUIDE v2.5 05-Sep-2025 15:01:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cal_node_detect_OpeningFcn, ...
                   'gui_OutputFcn',  @Cal_node_detect_OutputFcn, ...
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



% --- Executes just before Cal_node_detect is made visible.
function Cal_node_detect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cal_node_detect (see VARARGIN)

% Choose default command line output for Cal_node_detect
handles.output = hObject;
set(handles.axes1, 'XTick', [], 'YTick', [],'XColor', 'none', 'YColor', 'none'); 
handles.dis_therehold = 0;
handles.core_precent = 0.2;
set(handles.dis_there_slider,'Value',0);
set(handles.core_precent_slider,'Value',0.2);
set(handles.dis_there_edit,'String',0);
set(handles.core_precent_edit,'String',0.2);
%定义统计全局变量
handles.core_para = table( [], [], [] ,  'VariableNames', {'核心数量', '非核心数量', '核心占比'}  );
handles.core_contact = table([], [], [],'VariableNames', {'核心之间', '核心非核心之间', '非核心之间'});
handles.average_cor = table([], [], [],'VariableNames', {'核心之间', '核心非核心之间', '非核心之间'});
handles.average_dis = table([], [], [],'VariableNames', {'核心之间', '核心非核心之间', '非核心之间'});
handles.coredepart_threshold_all = [];
handles.cor_dis_record = {{};{}};
guidata(hObject, handles);

% UIWAIT makes Cal_node_detect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Cal_node_detect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    fprintf("正在读取...\n");
    [file, path] = uigetfile('*.*');
    set(handles.data_adress_text, 'String', strcat(path,'\',file));
    if sum(file ~= 0) ~= 0 || sum(path ~= 0) ~= 0
        data_path = [path, file];
        [~,~,ext] = fileparts(data_path);
        
         
        if strcmp(ext, '.avi')
            [data, dataInfo] = readAviSeq(data_path);
             contrast_value = get(handles.contrast_value, 'string');
             contrast_value = str2double(contrast_value);

            if contrast_value > 0
              contrastRatio = contrast_value;                       
            else
                contrastRatio = 0.15;%默认对比度0.15
                set(handles.contrast_value, 'String', 0.15);
            end
            handles.video.videoData = data;
            handles.video.videoInfo = dataInfo;
            handles.data_path = data_path;
            [path_dir, ~, ~] = fileparts(fileparts(data_path));
            handles.path_dir = path_dir;
            dataDouble = double(handles.video.videoData);
            meanData = uint8(mean(dataDouble, 3));
            handles.meanData = meanData;
            handles.rawData = handles.meanData;
            handles.contrastRatio = contrastRatio;
           
            handles.rawData_original = handles.rawData;%返回切割切割前细胞所用
            

            handles.contrastData = adapthisteq(handles.rawData, 'ClipLimit', contrastRatio);
            handles.contrastData_original = handles.contrastData;
            handles.video_original = handles.video;



            axes(handles.axes1)
            imshow(handles.contrastData)
             % 显示路径
            set(handles.data_adress_text, 'String', data_path);
            fprintf("读取完成\n")
        else
            fprintf("不支持的文件类型\n")
        end    
    end
catch me
    fprintf("读取失败\n")
    rethrow(me)
end
guidata(hObject, handles);



% --- Executes on button press in enlarge.
function enlarge_Callback(hObject, eventdata, handles)
% hObject    handle to enlarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%删除之前细胞参数
handles.contrastData = handles.contrastData_original;
handles.rawData = handles.rawData_original;
handles.video = handles.video_original;
set(handles.enlarge, 'backgroundcolor', [1.00, 0.00, 0.00])%设位置按键颜色
set([handles.main_sturrture,handles.main_sturrture,handles.mdbutton...
    ,handles.signal_detect,handles.node_detect],...
    'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'off') % 按键可用性设置

if isfield(handles, 'Cell')%清除该图片中上一个细胞数据
    handles = rmfield(handles, 'Cell');
end
if isfield(handles, 'mainsturcture')
    handles = rmfield(handles, 'mainsturcture');
end
if isfield(handles, 'microdomain')
    handles = rmfield(handles, 'microdomain');
end
% if isfield(handles, 'mask')
%     handles = rmfield(handles, 'mask');
% end
if isfield(handles, 'boundary')
    handles = rmfield(handles, 'boundary');
    handles = rmfield(handles, 'StructName');
end
if isfield(handles, 'rawSignal')
    handles = rmfield(handles, 'rawSignal');
    handles = rmfield(handles, 'brounary_micro');
    handles = rmfield(handles, 'rawSignal_micro');
    handles = rmfield(handles, 'dfSignal_micro');
    handles = rmfield(handles, 'calcium_info');
end


axes(handles.axes1)
cla
imshow(handles.contrastData)
contrastData = handles.contrastData;
rawData = handles.rawData;
video = handles.video;
roi = images.roi.Freehand;
roi = drawfreehand;
mask = createMask(roi);
[B,~,~,~] = bwboundaries(mask,8);
coordinates = B(1,1);
coordinates = cell2mat(coordinates);



x_1 = max(coordinates(:,1));
x_2 = min(coordinates(:,1));
y_1 = max(coordinates(:,2));
y_2 = min(coordinates(:,2));

% 保存截取图像位置，即xy初始截取长度
handles.cutlength.x = x_2-1;
handles.cutlength.y = y_2-1;
handles.cutlength.x_2 = x_2;
handles.cutlength.y_2 = y_2;
handles.cutlength.x_1 = x_1;
handles.cutlength.y_1 = y_1;

contrastData = contrastData(x_2:x_1,:);
contrastData = contrastData(:,y_2:y_1);
rawData = rawData(x_2:x_1,:);
rawData = rawData(:,y_2:y_1);
video.videoData = handles.video.videoData(x_2:x_1,y_2:y_1,:);


axes(handles.axes1)
clc
imshow(contrastData)
handles.contrastData = contrastData;
handles.rawData = rawData;
handles.video = video;
set(handles.enlarge, 'backgroundcolor', [0.94, 0.94, 0.94])
set([handles.cellbutton,handles.reset_sturct], 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on') % 按键可用性设置
guidata(hObject, handles);



% --- Executes on button press in cellbutton.
function cellbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cellbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.cellbutton, 'backgroundcolor', [1.00, 0.00, 0.00])%设位置按键颜色
axes(handles.axes1)
cla
imshow(handles.contrastData)
roi = images.roi.Freehand;
roi = drawfreehand;
mask = createMask(roi);
[B,~,~,~] = bwboundaries(mask,8);

% 单细胞，如果连通域过多，取最大的

if length(B) > 1
    k = 1;
    for i = 1 : length(B)
        if length(B{k}) < length(B{i})
            k = i;
        end
    end
    handles.curBoundary = B(k);
else
    handles.curBoundary = B;    % 保存当前的连通域
end

% handles.boundary = handles.curBoundary; % 所有连通域

handles.Cell.CellMask = mask;  % 单细胞掩膜

handles.Cell.CellBoundary = handles.curBoundary;   % 单细胞连通域
handles.Cell.StructName = {'Cell'}; % 保存名称
[boundary, StructName] = getAllBoudary(handles);
handles.boundary = boundary;
handles.StructName = StructName;
StructShow(hObject, eventdata, handles, boundary, StructName)
set(handles.cellbutton, 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'off') 
set(handles.main_sturrture, 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on') % 按键可用性设置

guidata(hObject, handles);



% --- Executes on button press in main_sturrture.
function main_sturrture_Callback(hObject, eventdata, handles)
% hObject    handle to main_sturrture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.main_sturrture, 'backgroundcolor', [1.00, 0.00, 0.00])%设位置按键颜色

roi = images.roi.Freehand;
roi = drawfreehand;
handmask = createMask(roi);
if ~isfield(handles, 'mask')
    mask = handmask.*handles.Cell.CellMask;
else
    mask = handmask.*handles.mask.cellDelMS;
end
[B,L,~,~] = bwboundaries(mask,8);

if length(B) > 1
    k = 1;
    for i = 1 : length(B)
        if length(B{k}) < length(B{i})
            k = i;
        end
    end
    handles.curBoundary = B(k);
else
    handles.curBoundary = B;    % 保存当前的连通域
end

if ~isfield(handles, 'mainsturcture')
    handles.mainsturcture.Mask = mask;  % 单细胞掩膜
    handles.mainsturcture.Boundary = handles.curBoundary;   % 单细胞连通域
    handles.mainsturcture.StructName = {'Mainsturcture'};
else
    handles.mainsturcture.Mask = cat(3, handles.mainsturcture.Mask, mask);
    handles.mainsturcture.Boundary = [handles.mainsturcture.Boundary; handles.curBoundary];
    handles.mainsturcture.StructName = [handles.mainsturcture.StructName; {'Mainsturcture'}]; % 保存名称
end

if ~isfield(handles, 'mask')
    handles.mask.cellDelMS = handles.Cell.CellMask -mask;
else
    handles.mask.cellDelMS = handles.mask.cellDelMS - mask;
end

[boundary, StructName] = getAllBoudary(handles);
handles.boundary = boundary;
handles.StructName = StructName;

axes(handles.axes1)
imshow(handles.contrastData), hold on

boundaryAndTextOnImage(boundary, 0,0, []), hold off
handles.structName = StructName;
guidata(hObject, handles);

set(handles.mdbutton, 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on') % 按键可用性设置


guidata(hObject, handles);





% --- Executes on button press in mdbutton.
function mdbutton_Callback(hObject, eventdata, handles)
% hObject    handle to mdbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.mdbutton, 'backgroundcolor', [1.00, 0.00, 0.00])%设位置按键颜色
handles.mainsturtboundary = handles.boundary;
zoom_vaule = 4;
zoom_mutiple = zoom_vaule/2;
zoom_mutiple =  zoom_mutiple*zoom_mutiple;

 %handles.mask.otherMask = handles.mask.cellDelSomaAndMainBranch;
%  if  strcmp(handles.StructName{end, 1} ,'End Feet')
%      handles.mask.otherMask = handles.mask.cellDelSoma;
%      for i = 1 : length(handles.mainBranch.mainBranchMask(1, 1, :))
%          handles.mask.otherMask = handles.mask.otherMask - handles.mainBranch.mainBranchMask(:, :, i);
%      end
%      for i = 1 : length(handles.endFeet.endFeetMask(1, 1, :))
%          handles.mask.otherMask = handles.mask.otherMask - handles.endFeet.endFeetMask(:, :, i);
%      end
% 
%  else
%      handles.mask.otherMask = handles.mask.cellDelSoma;
%      for i = 1 : length(handles.mainBranch.mainBranchMask(1, 1, :))
%          handles.mask.otherMask = handles.mask.otherMask - handles.mainBranch.mainBranchMask(:, :, i);
%      end
%  end
%         case 1  % 自动检测
if isfield(handles, 'microdomain')
    handles = rmfield(handles, 'microdomain');  % 清除
end
handles.otherData = double(handles.contrastData).*  handles.mask.cellDelMS;
tmp_max = double(max(handles.otherData, [], 'all'));    % 最大像素点
image_normalization = handles.otherData / tmp_max;   % 归一化
handles.image_normalization = image_normalization;

meanImage = mean(handles.image_normalization(handles.image_normalization ~= 0));  % 均值
stdImage = std(handles.image_normalization(handles.image_normalization ~= 0));    % 标准差


% 最小强度阈值
thresholdIntensityMin = meanImage + 0.5 * stdImage;
% 最小面积阈值
minAreaThreshold = 10*zoom_mutiple;
% 最大面积阈值
maxAreaThreshold  = 200*zoom_mutiple;

[B, mask] = autoDetectMicrodomain(handles.image_normalization, thresholdIntensityMin, minAreaThreshold, maxAreaThreshold);
StructName = fsyncSaveMicrodomainPara( B);
handles.microdomain.StructName = StructName;
handles.microdomain.microdomainBoundary = B;
handles.microdomain.para.thresholdIntensityMin = thresholdIntensityMin;
handles.microdomain.para.minAreaThreshold = minAreaThreshold;
handles.microdomain.para.maxAreaThreshold = maxAreaThreshold;
%         case 2  % 手动检测


[boundary, StructName] = getAllBoudary(handles);  % 保存所有boundary
handles.boundary = boundary;
handles.StructName = StructName;
axes(handles.axes1)
imshow(handles.contrastData), hold on
boundaryAndTextOnImage(boundary, 0,0, StructName), hold off

guidata(hObject, handles);

StructShow(hObject, eventdata, handles, boundary, StructName)
set([handles.mdbutton,handles.main_sturrture],'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'off') 
set(handles.signal_detect, 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on') % 按键可用性设置





% --- Executes on button press in signal_detect.
function signal_detect_Callback(hObject, eventdata, handles)
% hObject    handle to signal_detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% df/f
[rawSignal] = rawSignalOnROI(handles.video.videoData, handles.boundary);
% [rawSignal] = rawSignalOnROI(handles.video.videoData, edge_first_ROI);
calcium_signal = getDF(rawSignal, 5);
% 检测钙信号参数
for i = 1 : length(calcium_signal(1,:))
    dfSignal = calcium_signal(:,i);
    filter_low_start = 0.05;
    filter_low_finish = 0.05;
    [start_tmp_merge, peak_tmp_merge, finish_tmp_merge] = ...
        autoDetectCalcium(dfSignal,  filter_low_start, filter_low_finish);

    calcium_info{i}.signal.dfSignal = dfSignal;
    calcium_info{i}.calciumIndex.start = start_tmp_merge;
    calcium_info{i}.calciumIndex.peak = peak_tmp_merge;
    calcium_info{i}.calciumIndex.finish = finish_tmp_merge;
    %     calcium_info{i}.ROIInfo = classfied_microdomain_struct{i, 1}.ROIInfo;
    waveParams = WaveParams(dfSignal, start_tmp_merge, finish_tmp_merge, 5);
    calcium_info{i}.waveParams = waveParams;

end

handles.brounary_micro = handles.microdomain.microdomainBoundary; 
handles.rawSignal_micro = rawSignal(:,end-length(handles.microdomain.microdomainBoundary)+1:end);
handles.dfSignal_micro = calcium_signal(:,end-length(handles.microdomain.microdomainBoundary)+1:end);


handles.calcium_signal = calcium_signal;
handles.calcium_info = calcium_info;
handles.rawSignal = rawSignal;
% 展示
visibleCalcium(calcium_info);
fprintf("已绘制全部钙信号\n")


% handles.distance_soma_micro = distance_ROI(handles.brounary_soma_micro);

set(handles.node_detect, 'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on')
guidata(hObject, handles);



% --- Executes on button press in reset_sturct.
function reset_sturct_Callback(hObject, eventdata, handles)
% hObject    handle to reset_sturct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mainsturcture')%清除已画细胞数据
    handles = rmfield(handles, 'mainsturcture');
    handles = rmfield(handles, 'mainsturtboundary');
end
if isfield(handles, 'microdomain')
    handles = rmfield(handles, 'microdomain');
end
% if isfield(handles, 'mask')
%     handles = rmfield(handles, 'mask');
% end
if isfield(handles, 'boundary')
    handles = rmfield(handles, 'boundary');
    handles = rmfield(handles, 'StructName');
end


set([handles.main_sturrture,handles.mdbutton], ...
    'backgroundcolor', [0.94, 0.94, 0.94], 'enable', 'on')

axes(handles.axes1)
imshow(handles.contrastData);
% boundaryAndTextOnImage(handles.Cell.CellBoundary, 0,0, {0,1});hold off
set(handles.sturturetable, 'data', {})

[boundary, StructName] = getAllBoudary(handles);
handles.boundary = boundary;
handles.StructName = StructName;
StructShow(hObject, eventdata, handles, boundary, StructName)
guidata(hObject, handles);





% --- Executes on button press in node_detect.
function node_detect_Callback(hObject, eventdata, handles)
handles = astromicronet_gui_analyze(handles);
astromicronet_gui_show(handles,'core');
guidata(hObject,handles);

function core_sort_Callback(hObject, eventdata, handles)
assert(isfield(handles,'astromicronetResult'),'Run node detection first.');
astromicronet_gui_show(handles,'network');
guidata(hObject,handles);

function core_stats_Callback(hObject, eventdata, handles)
handles = astromicronet_gui_statistics(handles);
astromicronet_gui_show(handles,'statistics');
guidata(hObject,handles);

function para_active_Callback(hObject, eventdata, handles)
handles = astromicronet_gui_accumulate(handles);
guidata(hObject,handles);

function export_stats_Callback(hObject, eventdata, handles)
astromicronet_gui_export(handles,'statistics');

function plot_image_Callback(hObject, eventdata, handles)
astromicronet_gui_export_image(handles);

function export_signal_Callback(hObject, eventdata, handles)
astromicronet_gui_export(handles,'signals');

function export_node_Callback(hObject, eventdata, handles)
astromicronet_gui_export(handles,'nodes');

function contrast_value_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contrast_value as text
%        str2double(get(hObject,'String')) returns contents of contrast_value as a double


% --- Executes during object creation, after setting all properties.
function contrast_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function dis_there_slider_Callback(hObject, eventdata, handles)
% hObject    handle to dis_there_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
dis_therehold_Str = get(handles.dis_there_slider, 'Value');
dis_therehold = dis_therehold_Str;
set(handles.dis_there_edit, 'string', dis_therehold,'min', 0, 'max', 1);
handles.dis_therehold = dis_therehold;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function dis_there_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dis_there_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function dis_there_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dis_there_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dis_there_edit as text
%        str2double(get(hObject,'String')) returns contents of dis_there_edit as a double
dis_therehold_Str = get(handles.dis_there_edit, 'string');
dis_therehold = str2double(dis_therehold_Str);
set(handles.dis_there_slider, 'Value', dis_therehold,'min', 0, 'max', 1);
handles.dis_therehold = dis_therehold;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dis_there_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dis_there_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function core_precent_slider_Callback(hObject, eventdata, handles)
% hObject    handle to core_precent_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
core_precent_Str = get(handles.core_precent_slider, 'Value');
core_precent = core_precent_Str;
set(handles.core_precent_edit, 'string', core_precent,'min', 0, 'max', 1);
handles.core_precent = core_precent;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function core_precent_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to core_precent_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function core_precent_edit_Callback(hObject, eventdata, handles)
% hObject    handle to core_precent_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of core_precent_edit as text
%        str2double(get(hObject,'String')) returns contents of core_precent_edit as a double
core_precent_Str = get(handles.core_precent_edit, 'string');
core_precent = str2double(core_precent_Str);
set(handles.core_precent_slider, 'Value', core_precent,'min', 0, 'max', 1);
handles.core_precent = core_precent;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function core_precent_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to core_precent_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes when entered data in editable cell(s) in sturturetable.
function sturturetable_CellEditCallback(hObject, eventdata, handles)

StructName = get(handles.sturturetable, 'data');
StructShow(hObject, eventdata, handles,  handles.microdomain.microdomainBoundary, StructName)

guidata(hObject, handles);





% 使用函数
function StructShow(hObject, eventdata, handles, boundary, StructName)

md_StructName = StructName(cell2mat(StructName(:,2))==1,:);
set(handles.sturturetable, 'data', md_StructName)
axes(handles.axes1)
imshow(handles.contrastData), hold on

boundaryAndTextOnImage(boundary, 0,0, StructName), hold off
handles.structName = StructName;
guidata(hObject, handles);




function [boundary, StructName] = getAllBoudary(handles)

if isfield(handles, 'boundary')
    handles = rmfield(handles, 'boundary');
end
if isfield(handles, 'StructName')
    handles = rmfield(handles, 'StructName');
end

boundary = cell(0);  % 边界
StructName = cell(0);    % 名称 + 是否显示 

if isfield(handles, 'Cell')
    
    boundary = [boundary; handles.Cell.CellBoundary];
    for i = 1 : length(handles.Cell.StructName)
        StructName = [StructName; {handles.Cell.StructName{i}, true}];
    end
end

if isfield(handles, 'mainsturcture')
    
    boundary = [boundary; handles.mainsturcture.Boundary];
    for i = 1 : length(handles.mainsturcture.StructName)
        StructName = [StructName; {handles.mainsturcture.StructName{i}, false}];
    end
end

if isfield(handles, 'microdomain')
    
    boundary = [boundary; handles.microdomain.microdomainBoundary];
    for i = 1 : length(handles.microdomain.StructName)
        StructName = [StructName; {handles.microdomain.StructName{i},true}];
    end
end



function StructName = fsyncSaveMicrodomainPara( B)
% 保存微域的参数
if ~isempty(B)
    for i = 1 : length(B)
        if i == 1

            StructName = {'microdomain'}; % 保存名称
        else
            StructName = [StructName; {'microdomain'}]; % 保存名称
        end
    end
else
    StructName = [];
    fprintf("未检测到microdomain\n");
end
