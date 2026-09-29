function paths=astromicronet_gui_export_image(handles,filePath)
% Export visible axes with vector boundaries and a raster microscopy background.
if nargin<2
    [file,folder]=uiputfile('*.pdf','Save AstroMicroNet image','AstroMicroNet_image.pdf');
    if isequal(file,0),paths={};return;end
    filePath=fullfile(folder,file);
end
[folder,name]=fileparts(filePath);base=fullfile(folder,name);
f=figure('Visible','off','Color','w');clean=onCleanup(@()close(f));
ax=axes(f);copyobj(allchild(handles.axes1),ax);
set(ax,'XLim',get(handles.axes1,'XLim'),'YLim',get(handles.axes1,'YLim'), ...
    'YDir',get(handles.axes1,'YDir'),'DataAspectRatio',get(handles.axes1,'DataAspectRatio'));
axis(ax,'off');
exportgraphics(ax,[base '.pdf'],'ContentType','vector');
print(f,[base '.svg'],'-dsvg','-painters');
savefig(f,[base '.fig']);
paths={[base '.pdf'],[base '.svg'],[base '.fig']};
fprintf('Saved image: %s\n',base);
end
