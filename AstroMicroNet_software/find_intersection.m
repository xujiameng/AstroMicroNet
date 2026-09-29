function [valley_loc]  = find_intersection(data, threshold)
% 找谷值
% 对于谷值，左边点和右边点均比谷值大，认为是谷值
% data:一维数据
% data = data_not(k, peak_tmp{k, 1}(i) : peak_tmp{k, 1}(i + 1));
% valley_value = [];
valley_loc = [];
data_len = length(data);
valley_cnt = 0;
for i = 1 : data_len - 1
    if data(i) == threshold || (data(i) < threshold && data(i + 1) > threshold)
        valley_cnt = valley_cnt + 1;
        valley_loc(valley_cnt) = i;
    elseif data(i) > threshold && data(i + 1) < threshold
        valley_cnt = valley_cnt + 1;
        valley_loc(valley_cnt) = i + 1;
    end
end

% 如果没有交点， 选最小值
if valley_cnt == 0
    [~, valley_loc] = min(data);
end
