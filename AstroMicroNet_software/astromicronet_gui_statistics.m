function handles=astromicronet_gui_statistics(handles)
% Summarize each unordered ROI pair once, including valid empty-core results.
R=handles.correlation_group_micro;D=handles.distance_group_micro;
core=logical(handles.astromicronetResult.final(:));N=numel(core);
assert(isequal(size(R),[N N]) && isequal(size(D),[N N]),'Network arrays are not aligned.');
[i,j]=find(triu(true(N),1));edge=sub2ind([N N],i,j);
kind=double(core(i))+double(core(j));retained=handles.netcontact_matrix(edge);
names={'coco','counco','uncounco'};codes=[2 1 0];
counts=zeros(1,3);possible=counts;density=nan(1,3);meanCor=density;meanDis=density;
for k=1:3
    selected=kind==codes(k) & retained;
    weights=reshape(R(edge(selected)),1,[]);distances=reshape(D(edge(selected)),1,[]);
    counts(k)=sum(selected);possible(k)=sum(kind==codes(k));
    if possible(k)>0,density(k)=counts(k)/possible(k);end
    if counts(k)>0,meanCor(k)=mean(weights);meanDis(k)=mean(distances);end
    handles.(['num_contact_' names{k}])=counts(k);
    handles.(['cor_' names{k}])=weights;handles.(['dis_' names{k}])=distances;
    handles.(['average_cor_' names{k}])=meanCor(k);handles.(['average_dis_' names{k}])=meanDis(k);
end
handles.num_core=sum(core);handles.num_uncore=N-sum(core);
handles.connection_density=density;
handles.micro_group=cell(sum(core),1);
C=find(core);P=find(~core);
if ~isempty(C)
    for p=reshape(P,1,[])
        values=R(p,C);values(~handles.netcontact_matrix(p,C) | values<=0)=-Inf;
        [value,which]=max(values);
        if isfinite(value),handles.micro_group{which}(end+1)=p;end
    end
end
handles.core_cantact_counts=cellfun(@numel,handles.micro_group);
handles.unique_counts=unique(handles.core_cantact_counts);
handles.core_cantact_frequency=arrayfun(@(x)sum(handles.core_cantact_counts==x),handles.unique_counts);
handles.network_metrics=struct('N',N,'C',sum(core),'core_fraction',mean(core), ...
    'categories',{{'CC','CP','PP'}},'edge_counts',counts,'possible_pairs',possible, ...
    'density',density,'mean_correlation',meanCor,'mean_distance_um',meanDis);
end
