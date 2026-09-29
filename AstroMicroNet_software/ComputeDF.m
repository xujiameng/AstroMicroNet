function [df_f] = ComputeDF(fvalue,base_length,fs,display)
% fvalue: N X 1 vector;
% if  prctile(fvalue,25)<0.5
%     df_f=zeros(size(fvalue));
%     return;
% end
fvalue_comb=[fvalue(end-fs*base_length+1:end);fvalue];
df=fvalue_comb;
if display
    h = waitbar(0,'Please wait...');
end
for i=1:length(fvalue_comb)     
    if i<=(length(fvalue_comb)-fs*base_length+1)        
        f0=fvalue_comb(i:i+fs*base_length-1);
        x = prctile(f0,25);
        df(i)=fvalue_comb(i)/x-1;
    else
        f0=fvalue_comb(i:end);
        x = prctile(f0,25);
        df(i)=fvalue_comb(i)/x-1;
    end
    if display
        waitbar(i/length(fvalue_comb),h,[num2str(round(100*i/length(fvalue_comb))),'% done']);
    end
end
if display
    close(h);
end
df_f=df(fs*base_length+1:end);