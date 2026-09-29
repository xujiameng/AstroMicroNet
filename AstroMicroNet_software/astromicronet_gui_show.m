function figures=astromicronet_gui_show(handles,mode)
% Keep the historical image canvas and red/blue ROI conventions.
figures=gobjects(0);
if nargin<2,mode='core';end
assert(isfield(handles,'astromicronetResult'),'Run node detection first.');
if strcmp(mode,'statistics')
    m=handles.network_metrics;
    figures(1)=figure;bar([m.C,m.N-m.C]);set(gca,'XTick',1:2,'XTickLabel',{'core','periphery'});ylabel('Number of microdomains');
    figures(2)=figure;bar(m.edge_counts);set(gca,'XTick',1:3,'XTickLabel',m.categories);ylabel('Number of connections');
    figures(3)=figure;groupPlot({handles.cor_coco,handles.cor_counco,handles.cor_uncounco});ylabel('Normalized cross-correlation');
    figures(4)=figure;groupPlot({handles.dis_coco,handles.dis_counco,handles.dis_uncounco});ylabel('Distance (μm)');
    return
end
ax=handles.axes1;
cla(ax);imshow(handles.contrastData,'Parent',ax);hold(ax,'on');
if strcmp(mode,'network')
    [i,j]=find(triu(handles.netcontact_matrix,1));c=handles.boundary_center_micro;
    core=logical(handles.astromicronetResult.final(:));colors=[.10 .45 .70;.51 .43 .69;.79 .12 .11];
    for k=1:numel(i)
        kind=double(core(i(k)))+double(core(j(k)))+1;
        plot(ax,c([i(k),j(k)],2),c([i(k),j(k)],1),'Color',colors(kind,:),'LineWidth',0.5);
    end
end
core=logical(handles.astromicronetResult.final(:));
for k=1:numel(core)
    b=handles.analysisBoundaries{k};color=[.086 .447 .694];
    if core(k),color=[.792 .118 .114];end
    plot(ax,b(:,2),b(:,1),'Color',color,'LineWidth',2);
end
hold(ax,'off');
end

function groupPlot(values)
data=[];group=[];
for k=1:3
    v=values{k}(:);v=v(isfinite(v));data=[data;v];group=[group;k*ones(size(v))]; %#ok<AGROW>
end
if ~isempty(data),boxplot(data,group,'Positions',unique(group)','Symbol','');end
set(gca,'XTick',1:3,'XTickLabel',{'CC','CP','PP'},'XLim',[.5 3.5]);
end
