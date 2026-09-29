function handles = astromicronet_gui_analyze(handles)
% Adapt historical GUI data to the integrated analysis without altering controls.
assert(isfield(handles,'dfSignal_micro'),'AstroMicroNet:MissingSignals','Extract ROI signals first.');
X = double(handles.dfSignal_micro);
assert(isfield(handles,'brounary_micro'),'AstroMicroNet:MissingROIs','ROI boundaries are required.');
boundaries = handles.brounary_micro(:);
N = size(X,2);
assert(numel(boundaries)==N && size(X,1)>1 && all(isfinite(X(:))), ...
    'AstroMicroNet:InvalidSignals','Signals and ROI boundaries must be finite and aligned.');
centers=zeros(N,2);
for i=1:N
    b=double(boundaries{i});
    assert(size(b,2)==2 && ~isempty(b) && all(isfinite(b(:))),'Invalid ROI boundary.');
    centers(i,:)=mean(b,1);
end
opts=struct();
if isfield(handles,'astromicronetOptions'),opts=handles.astromicronetOptions;end
if ~isfield(opts,'rho'),opts.rho=handles.core_precent;end
if ~isfield(opts,'seedKey'),opts.seedKey='gui';end
if ~isfield(opts,'pixelSizeUm'),opts.pixelSizeUm=0.187;end
if ~isfield(opts,'roiIDs'),opts.roiIDs=(1:N)';end
assert(numel(opts.roiIDs)==N && all(diff(opts.roiIDs(:))>0),'ROI identifiers must be increasing.');
parent=1:N;
if isfield(handles,'rawSignal_micro')
    raw=handles.rawSignal_micro;
    assert(isequal(size(raw),size(X)),'Raw and processed ROI traces must be aligned.');
    for i=1:N
        for j=i+1:N
            if isequal(raw(:,i),raw(:,j)) && isequal(centers(i,:),centers(j,:)) && isequal(boundaries{i},boundaries{j})
                assert(isequal(X(:,i),X(:,j)),'AstroMicroNet:InconsistentDuplicate','Identical raw ROI records have different processed traces.');
                parent(j)=parent(i);
            end
        end
    end
end
keep=unique(parent,'stable');[~,mapping]=ismember(parent,keep);
opts.roiIDs=opts.roiIDs(keep);
opts.inputKind='X';
% A positive value in the historical connection control remains an input cutoff.
if isfield(handles,'dis_therehold') && handles.dis_therehold>0
    [~,R]=distance_correlation_ROI(boundaries(keep),X(:,keep));
    R(R<handles.dis_therehold)=0;R(1:size(R,1)+1:end)=0;
    opts.inputKind='R';input=R;
else
    input=X(:,keep);
end
result=astromicronet_analyze(input,opts);
result.gui_version='20260927_v11';
result.source_roi_ids=double(reshape(getSourceIDs(handles,N),[],1));
result.keep_indices=keep(:);result.source_to_unique=mapping(:);
result.pixel_size_um=opts.pixelSizeUm;
result.input_cutoff=handles.dis_therehold;
% Original GUI signal records remain available alongside the unique analysis rows.
handles.astromicronetResult=result;
handles.analysisBoundaries=boundaries(keep);
handles.analysis_dfSignal_micro=X(:,keep);
if isfield(handles,'rawSignal_micro'),handles.analysis_rawSignal_micro=handles.rawSignal_micro(:,keep);end
handles.coreIndices=result.coreIndices(:)';
handles.uncoreIndices=result.peripheryIndices(:)';
handles.correlation_group_micro=result.R;
handles.boundary_center_micro=centers(keep,:);
delta=permute(centers(keep,:),[1 3 2])-permute(centers(keep,:),[3 1 2]);
handles.distance_group_micro=sqrt(sum(delta.^2,3))*opts.pixelSizeUm;
handles.coredepart_threshold=double(result.theta);
handles.netcontact_matrix=result.R>=result.theta;
handles.netcontact_matrix(1:numel(keep)+1:end)=false;
handles=astromicronet_gui_statistics(handles);
end

function ids=getSourceIDs(handles,N)
ids=(1:N)';
if isfield(handles,'astromicronetOptions') && isfield(handles.astromicronetOptions,'roiIDs')
    ids=handles.astromicronetOptions.roiIDs;
end
end
