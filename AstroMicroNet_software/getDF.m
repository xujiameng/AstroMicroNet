function [smoothData] = getDF(rawSignal,fs)
    if isempty(fs)
    fs = 5;
    end
    step = 1500;%最初f窗口大小
    % 根据原始波形获取df/f,并进行滤波
    %% 处理 取每20帧最小值
    rawSignalLen = length(rawSignal(:, 1));
    df_f = rawSignal;
    for i = 1 : 1 : rawSignalLen
        tmp_right = min(rawSignalLen - i + 1, step / 2);
        tmp_left = min(i, step / 2);
        min200 = prctile(rawSignal(i - tmp_left + 1 : i + tmp_right - 1, :), 20)';   % 每一列最小值

%         for k = 1 :  length(rawSignal(1,:))
%             if i<251
%                 min20(k) = mean(rawSignal(1  : 500, k)); 
%             elseif i<rawSignalLen-250
%                 min20(k) = mean(rawSignal(i-250 : i +250, k));   % 每一列最小值
%             else
%                 min20(k) = mean(rawSignal(rawSignalLen-250 : rawSignalLen, k));
%             end
%         end
%         for k = 1 :  length(rawSignal(1,:))
%             if i<251
%                 min200(k) = mean(rawSignal(1  : 500, k));   % 每一列最小值
%             elseif i<rawSignalLen-250
%                 min200(k) = mean(rawSignal(i-250  : i+250, k));
%             else
%                 min200(k) = mean(rawSignal(rawSignalLen-500:rawSignalLen , k));
%             end
%         end
%         for k = 1 :  length(rawSignal(1,:))
%             if i<101
%                 min200(k) = mean(rawSignal(1  : 200, k));   % 每一列最小值
%             elseif i<rawSignalLen-100
%                 min200(k) = mean(rawSignal(i-100  : i+100, k));
%             else
%                 min200(k) = mean(rawSignal(rawSignalLen-200:rawSignalLen , k));
%             end
%         end

        % df/f0
        f0 = min200;
        for j = 1 : length(f0)
            df_f(i, j) = (rawSignal(i, j) - f0(j)) / f0(j);
        end
    end
    
    %% 平滑滤波(暂时去掉）
    smoothData = df_f;
    for k = 1 : length(rawSignal(1, :))
        smoothData(:, k) = smooth(df_f(:, k), 20);
%         smoothData(:, k) = movmean(df_f(:, k), fs);
%        
    end
end

