function filePath=astromicronet_gui_export(handles,kind,filePath)
% Preserve the historical MAT export, including the final labels and provenance.
if nargin<3
    [file,folder]=uiputfile('*.mat','Save AstroMicroNet data',[kind '.mat']);
    if isequal(file,0),filePath='';return;end
    filePath=fullfile(folder,file);
end
if strcmp(kind,'signals')
    assert(isfield(handles,'calcium_signal'),'Analyze signals before exporting.');
    save_data=selectFields(handles,{'calcium_info','rawSignal','calcium_signal','boundary','StructName'}); %#ok<NASGU>
    save(filePath,'save_data','-v7.3');
    fprintf('Saved %s\n',filePath);return;
end
assert(isfield(handles,'astromicronetResult'),'Run node detection before exporting.');
if strcmp(kind,'statistics')
    keys={'core_para','core_contact','average_cor','average_dis','cor_dis_record', ...
        'coredepart_threshold_all','connection_density_all','analysis_records','network_metrics'};
    outputpara=selectFields(handles,keys); %#ok<NASGU>
    save(filePath,'outputpara','-v7.3');
else
    keys={'boundary','StructName','soma','mainBranch','endFeet','microdomain','calcium_info', ...
        'rawSignal','calcium_signal','video','meanData','rawData','contrastData','contrastRatio', ...
        'curBoundary','mask','cutlength','coreIndices','uncoreIndices','distance_group_micro', ...
        'boundary_center_micro','correlation_group_micro','coredepart_threshold','netcontact_matrix', ...
        'astromicronetResult','network_metrics','analysisBoundaries'};
    save_data=selectFields(handles,keys);
    save_data.source_dfSignal_micro=handles.dfSignal_micro;
    save_data.dfSignal_micro=handles.analysis_dfSignal_micro;
    if isfield(handles,'rawSignal_micro'),save_data.source_rawSignal_micro=handles.rawSignal_micro;end
    if isfield(handles,'analysis_rawSignal_micro'),save_data.rawSignal_micro=handles.analysis_rawSignal_micro;end
    save_data.brounary_micro=handles.analysisBoundaries;
    save(filePath,'save_data','-v7.3');
end
fprintf('Saved %s\n',filePath);
end

function out=selectFields(s,names)
out=struct();
for k=1:numel(names)
    if isfield(s,names{k}),out.(names{k})=s.(names{k});end
end
end
