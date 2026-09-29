function [boundary,areaBoundary] = autoDetectMicrodomainFirst(image,minThresholdIntensity)
    % minThresholdIntensity: 最小强度阈值
    % 进行第一次检测
%     figure, imshow(image)
    [m, n] = size(image);
    %% 根据阈值找活跃像素点
    activateImageMinThreshold = image;
    activateImageMinThreshold(activateImageMinThreshold < minThresholdIntensity) = 0;
%     figure, imshow(activateImageMinThreshold)
    
    %% 确定连通域
    [B,~,~,A] = bwboundaries(activateImageMinThreshold,8, 'noholes'); % 8邻域找连通域
    
    %% 去掉子连通域
    [child, ~] = find(A == 1);
    boundary = cell(0); % 第一次检测连通域
    for i = 1 : length(B)
        if isempty(find(child == i, 1))
            boundary = [boundary; B{i}];
        end
    end
    % 计算连通域面积
    areaBoundary = [];
    for i = 1 : length(boundary)
        areaBoundary = [areaBoundary; regionprops(boundary{i} > 0, 'Perimeter').Perimeter];
    end
end

