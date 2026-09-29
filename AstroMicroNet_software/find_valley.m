function [valley_value,valley_loc]  = find_valley(data, threshold,baseline)
% 找谷值
% 对于谷值，左边点和右边点均比谷值大，认为是谷值
% data:一维数据
% data = data_not(k, peak_tmp{k, 1}(i) : peak_tmp{k, 1}(i + 1));
if ~exist('threshold', 'var')
    threshold = max(data);
end
valley_value = [];
valley_loc = [];
data_len = length(data);
valley_cnt = 0;
for i = 2 : data_len - 1
    if abs(data(i)-baseline) <= threshold && data(i) <= data(i - 1) && data(i) <= data(i + 1)
        valley_cnt = valley_cnt + 1;
        valley_value(valley_cnt) = data(i);
        valley_loc(valley_cnt) = i;
    end
end
