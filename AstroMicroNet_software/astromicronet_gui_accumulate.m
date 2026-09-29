function handles=astromicronet_gui_accumulate(handles)
% Append one current cell to the historical aggregate tables.
m=handles.network_metrics;
values={[m.C,m.N-m.C,m.core_fraction],m.edge_counts,m.mean_correlation,m.mean_distance_um};
names={'core_para','core_contact','average_cor','average_dis'};
for k=1:numel(names)
    name=names{k};t=handles.(name);
    handles.(name)=[t;array2table(values{k},'VariableNames',t.Properties.VariableNames)];
end
if ~isfield(handles,'connection_density_all'),handles.connection_density_all=zeros(0,3);end
handles.connection_density_all(end+1,:)=m.density;
handles.coredepart_threshold_all(end+1)=handles.coredepart_threshold;
cor={handles.cor_coco;handles.cor_counco;handles.cor_uncounco};
dis={handles.dis_coco;handles.dis_counco;handles.dis_uncounco};
handles.cor_dis_record{1}{end+1}=cor;handles.cor_dis_record{2}{end+1}=dis;
if ~isfield(handles,'analysis_records'),handles.analysis_records={};end
handles.analysis_records{end+1}=struct('metrics',m,'roi_ids',handles.astromicronetResult.roi_ids, ...
    'final',handles.astromicronetResult.final,'summary',handles.astromicronetResult.summary);
end
