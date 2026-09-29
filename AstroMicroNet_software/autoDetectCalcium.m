function [start_tmp_merge, peak_tmp_merge, finish_tmp_merge] = autoDetectCalcium(data,  filter_low_start, filter_low_finish)
%   Inputs
%   ------
%   data              : Nx1 vector (signal trace).
%   filter_low_start  : threshold used when searching for start points.
%   filter_low_finish : threshold used when searching for finish points.
%
%   Outputs
%   -------
%   start_tmp_merge   : start indices of detected events.
%   peak_tmp_merge    : peak indices of detected events.
%   finish_tmp_merge  : finish indices of detected events.
%   Notes
%   -----
%   This function depends on helper functions:
%     - find_valley(...)
%     - find_intersection(...)



%% 设置阈值为3SD   基线为均值
   filter_high = mean(data)+ 3 * std(data);
   filter_baseline = mean(data);
   bug_therehold = mean(data)+ 7 * std(data);
%% 找峰值点
    [~, peak_tmp] = findpeaks(data, 'MinPeakHeight', filter_high);
    
    %% 找起始点
    peak_len = length(peak_tmp);    % 峰值点数量
    start_tmp = zeros(peak_len, 1); % 谷值点
    if peak_len ~= 0
        % 检测每个峰值点的起始点和结束点, 利用谷值 和 交点
        for i = 1 : peak_len
            % 可能是起始点和结束点的谷值
            % 如果只有1个峰值点
            if peak_len == 1
                [~, start_loc_tmp] = find_valley(data(1 : peak_tmp(i), 1), filter_low_start,filter_baseline);     
                [start_intersect_tmp] = find_intersection(data(1 : peak_tmp(i), 1), filter_low_start);
            % 如果是第1个峰值点
            elseif i == 1   
                [~, start_loc_tmp] = find_valley(data(1 : peak_tmp(i), 1), filter_low_start,filter_baseline);
                [start_intersect_tmp] = find_intersection(data(1 : peak_tmp(i), 1),  filter_low_start);
            % 如果是最后1个峰值点
            elseif i == peak_len
                [~, start_loc_tmp] = find_valley(data(peak_tmp(i - 1) : peak_tmp(i), 1), filter_low_start,filter_baseline);
                start_loc_tmp = start_loc_tmp + peak_tmp(i - 1) - 1;

                [start_intersect_tmp] = find_intersection(data(peak_tmp(i - 1) : peak_tmp(i), 1),  filter_low_start);
                start_intersect_tmp = start_intersect_tmp + peak_tmp(i - 1) - 1;
            else   %else与else if无差别
                [~, start_loc_tmp] = find_valley(data(peak_tmp(i - 1) : peak_tmp(i), 1), filter_low_start,filter_baseline);
                start_loc_tmp = start_loc_tmp + peak_tmp(i - 1) - 1;

                [start_intersect_tmp] = find_intersection(data(peak_tmp(i - 1) : peak_tmp(i), 1),  filter_low_start);
                start_intersect_tmp = start_intersect_tmp + peak_tmp(i - 1) - 1;
            end
            % 如果 有起始点， 就选择右边第1个
            if ~isempty(start_loc_tmp)
                start_tmp(i, 1) = start_loc_tmp(end);
            elseif i == 1   % 如果没有起始点并且是第一个峰值点
                start_tmp(i, 1) = start_intersect_tmp(end);
            end
        end
    end
    
    %% 找结束点
    peak_len = length(peak_tmp);
    finish_tmp = zeros(peak_len, 1);
    if peak_len ~= 0      
        % 前 n - 1个结束点
        for i = 1 : peak_len - 1
            cnt = i;
            while cnt >= 1 && start_tmp(cnt) == 0
                cnt = cnt - 1;
            end
%             if start_tmp(cnt) ~= 0
            filter_low_finish_end = max(data(start_tmp(cnt), 1) + mean(data(:,1)) - std(data(:, 1)), filter_low_finish);
            [~, finish_tmp_k] = find_valley(data(peak_tmp(i) : peak_tmp(i + 1), 1), filter_low_finish_end,filter_baseline);
            finish_tmp_k = finish_tmp_k + peak_tmp(i) - 1;
%             end
            % 判断 第 i + 1 个起始点
            % 1. i + 1 的起始点为0——没有起始点
            % 2. 第 cnt 和 i＋1 的钙信号的峰值点之间，只有一个谷值
            % 3. 第 cnt 个钙信号的起始点 存在 
            % 4. 第 cnt 个钙信号的起始点 满足 起始点阈值要求，即不是通过此方法得到的
            if start_tmp(i + 1) == 0 && length(finish_tmp_k) == 1 && start_tmp(cnt) ~= 0 && data(start_tmp(cnt), 1) <= filter_low_start(1, 1)
                start_tmp(i + 1) = finish_tmp_k(end);
            end
            % 获取结束点
            if ~isempty(finish_tmp_k)
                finish_tmp(i, 1) = finish_tmp_k(1);
            end
        end
        % 检测最后1个钙信号
        cnt = peak_len;
        while cnt >= 1 && start_tmp(cnt) == 0
            cnt = cnt - 1;
        end
        if start_tmp(cnt) ~= 0
            filter_low_finish_end = max(data(start_tmp(cnt), 1) + mean(data(:,1)) - std(data(:, 1)), filter_low_finish);
            [~, finish_tmp_k] = find_valley(data(peak_tmp(peak_len) : end, 1), filter_low_finish_end,filter_baseline);
            finish_tmp_k = finish_tmp_k + peak_tmp(peak_len) - 1;
        end
        if ~isempty(finish_tmp_k)
            finish_tmp(peak_len, 1) = finish_tmp_k(1);
        end
    end
    
    %% 筛选
    start_tmp_merge = start_tmp;
    peak_tmp_merge = peak_tmp;
    finish_tmp_merge = finish_tmp;
%     %添加起止点差值阈值
%     delete_num = [];
%     for i = 1 : length(peak_tmp_merge)
%         data_start = start_tmp(i);
%         data_finish = finish_tmp(i);
%         if abs(data(data_start+1)-data(data_finish+1)) > filter_diff
%            delete_num = [delete_num;i];
%         end
%     end
%     % 删掉多余点
%     for i = 1 : length(delete_num)
%     start_tmp_merge(delete_num(i)-i+1) = [];
%     peak_tmp_merge(delete_num(i)-i+1) = [];
%     finish_tmp_merge(delete_num(i)-i+1) = [];
%     end

             
    peak_len = length(peak_tmp_merge);
    if peak_len ~= 0
        for i = 2 : peak_len
            % 如果start为0， 证明和前一个信号事同一个信号
            % 将当前的start 和前一个信号的finish 置0
            % 将两个钙信号中较小的置0
            if start_tmp_merge(i) == 0
                finish_tmp_merge(i - 1) = 0;
                if data(peak_tmp_merge(i - 1), 1) > data(peak_tmp_merge(i), 1)
                    peak_tmp_merge(i) = peak_tmp_merge(i - 1);
                end
                peak_tmp_merge(i - 1) = 0;
            end
        end
        % 第1个钙信号 没有起始点——去掉此种情况，没有起始点用交点
        % 最后1个钙信号没有结束点
        if finish_tmp_merge(peak_len) == 0
            peak_tmp_merge(peak_len) = 0;
            cnt = peak_len;
            while start_tmp_merge(cnt) == 0
                cnt = cnt - 1;
            end
            start_tmp_merge(cnt) = 0;
        end
    end
    
    not_ca_len = length(find(peak_tmp_merge == 0));
    peak_tmp_merge = sort(peak_tmp_merge);
    peak_tmp_merge = peak_tmp_merge(not_ca_len + 1: end);
    start_tmp_merge = sort(start_tmp_merge);
    start_tmp_merge = start_tmp_merge(not_ca_len + 1: end);
    finish_tmp_merge = sort(finish_tmp_merge);
    finish_tmp_merge = finish_tmp_merge(not_ca_len + 1: end);

end
    

