function [rawSignal] = rawSignalOnROI(video, boundary)
    % 根据ROI获取一维光强波形
    timeLength = length(video(1, 1, :)); % 时间序列
    ROINum = length(boundary);  % roi的数量
    mask = boundaryToMask(video(:, :, 1), boundary);    % 获取每个连通域的掩膜
    
    rawSignal = zeros(timeLength, ROINum);  % 原始波形，一列为一个ROI
    for i = 1 : length(boundary)
        maskIndex = find(mask(:, :, i) == 1);    % 找掩膜为1的索引
        for j = 1 : timeLength
            frame = video(:, :, j); % 当前帧
            frameSumIntensity = sum(frame(maskIndex));
            frameMeanIntensity = frameSumIntensity / length(maskIndex);
            rawSignal(j, i) = frameMeanIntensity;
        end
        fprintf("正在计算第%d个ROI\n", i)
    end
%     fprintf("计算完成\n")
%     save rawSignal rawSignal
%     figure, plot(rawSignal(:, 5));
end

