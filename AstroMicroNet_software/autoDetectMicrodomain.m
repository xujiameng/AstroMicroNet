function [finalBoundary, mask] = autoDetectMicrodomain(image, minThresholdIntensity, minAreaThreshold, maxAreaThreshold)
    % minThresholdIntensity: 最小强度阈值
    % minAreaThreshold: 最小面积阈值
    % maxAreaThreshold: 最大面积阈值
    % 根据三个参数，获取microdomain
    
    %%
    boundary = cell(0); % 所有不大的连通域
    boundaryArea = [];  % 记录保存的连通域的面积
    maxBoundary = cell(0);  % 所有大的连通域
    flag = 1;   % 第一次不保留低强度区域
    while 1
        if flag == 1
            % 第一次检测
            [boundaryTmp, boundaryAreaTmp] = autoDetectMicrodomainFirst(image, minThresholdIntensity);
            flag = 0;
        else
            % 第二次检测
            tmpBoundary = maxBoundary(1);
            maxBoundary = maxBoundary(2 : end); % 从队列中取出第一个
            [boundaryTmp, boundaryAreaTmp] = autoDetectMicrodomainSecond(image, tmpBoundary);
            % 如果计算出来的是原来的，强制面积为maxAreaThreshold
            if length(boundaryTmp) == 1 && isequal(boundaryTmp{1}, tmpBoundary{1}) && boundaryAreaTmp > maxAreaThreshold
                boundaryAreaTmp = maxAreaThreshold;
            end
        end
        maxIndex = find(boundaryAreaTmp > maxAreaThreshold);
        minIndex = find(boundaryAreaTmp <= maxAreaThreshold);
        boundary = [boundary; boundaryTmp(minIndex)];   % 小于最大面积阈值的放入队列
        boundaryArea = [boundaryArea; boundaryAreaTmp(minIndex)];   
        maxBoundary = [maxBoundary; boundaryTmp(maxIndex)]; % 大于面积阈值的继续分割
        if isempty(maxBoundary)
            break
        end

    end
    
    %%  去掉小面积连通域
    finalBoundary = cell(0);
    for i = 1 : length(boundaryArea)
        if boundaryArea(i) >= minAreaThreshold
            finalBoundary = [finalBoundary; boundary(i)];
        end
    end
    
    %% 每个连通域的掩膜
    mask = boundaryToMask(image, finalBoundary);
    
    
end

