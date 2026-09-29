function [B,areaBoundary] = autoDetectMicrodomainSecond(image, tmpBoundary)
    % tmpBoundary: 需要进行多次分割的连通域
    % minAreaThreshold: 最小面积阈值
    % maxAreaThreshold: 最大面积阈值
    % 进行第 n 次检测
    mask = boundaryToMask(image, tmpBoundary);  % 根据连通域获取掩膜
    mask = sum(mask, 3);
    mask(mask ~= 0) = 1;
    dataMask = image .* mask;   % 掩膜后的图像
    
    B = cell(0);    % 分割后的连通域
    areaBoundary = [];  % 分割后，每个连通域的面积
    
    % 均值和标准差
    meanData = mean(dataMask(dataMask ~= 0));
    stdData = std(dataMask(dataMask ~= 0));
    % 最用于分割连通域的强度
    minIntensityThreshold = meanData + 0.5 * stdData;
    %% 大于阈值的连通域
    maxDataMask = dataMask;
    maxDataMask(maxDataMask < minIntensityThreshold) = 0;
    maxMask = zeros(size(maxDataMask));
    % 确定连通域
    [BTmp,LTmp,~,A] = bwboundaries(maxDataMask,8, 'noholes'); % 8邻域找连通域
    % 如果只有一个，并且是 tmpBoundary ，保留，最外层函数中，强制使其面积为最大面积，不再参与下次分割
    if length(BTmp) == 1 && isequal(BTmp{1}, tmpBoundary{1})    
        B = BTmp;
    else
        rawIndex = 0;   % 记录和 tmpBoundary 一样的连通域的索引，不保留此连通域，同时保留此连通域下的子连通域
        for i = 1 : length(BTmp)
            if isequal(BTmp{i}, tmpBoundary{1})
                rawIndex = i;
                break;
            end
        end
        % 以 tmpBoundary 为父连通域的，不算做有父连通域
        if rawIndex ~= 0
            A(:, rawIndex) = 0;
        end
        % 去掉子连通域
        % 其中，父连通域不能是 rawIndex
        [child, parent] = find(A == 1);
        for i = 1 : length(BTmp)
            % 如果是最外层的连通域，舍弃
            if i ~= rawIndex && isempty(find(child == i, 1)) 
                B = [B; BTmp{i}];
            end
        end
        % 大于阈值部分的掩膜
        maxMask = boundaryToMask(image, B);
        maxMask = sum(maxMask, 3);
        maxMask(maxMask ~= 0) = 1;
        % 重新寻找一遍连通域，避免有连通域交错的情况
        [BTmp,LTmp,~,A] = bwboundaries(maxMask,8, 'noholes'); % 8邻域找连通域
        B = BTmp;
    end
    
    %%  小于阈值的连通域

    minDataMask = mask - maxMask;   % 原始掩膜 减去 大于阈值连通域的掩膜
    decreaseImage = minDataMask;
    decreaseImage(decreaseImage ~= 0) = 1;
    % 腐蚀
    decrease = imerode(decreaseImage,[1,1,1;1,1,1;1,1,1]);%3x3正方形结构元素的腐蚀
    % 腐蚀完膨胀
    increaseAfterDecrease = imdilate(decrease,[1,1,1;1,1,1;1,1,1]);%3x3正方形结构元素的膨胀

    % 确定连通域
    [BTmp,LTmp,~,A] = bwboundaries(increaseAfterDecrease,8, 'noholes'); % 8邻域找连通域
    % 同上
    if length(BTmp) == 1 && isequal(BTmp{1}, tmpBoundary{1})
        B = BTmp;
    else
        rawIndex = 0;   % 找最外面的连通域
        for i = 1 : length(BTmp)
            if isequal(BTmp{i}, tmpBoundary{1})
                rawIndex = i;
                break;
            end
        end

        % 以最外层连通域为父连通域的，不算做有父连通域
        if rawIndex ~= 0
            A(:, rawIndex) = 0;
        end

        % 去掉子连通域
        % 其中，父连通域不能是 rawIndex
        [child, parent] = find(A == 1);
        for i = 1 : length(BTmp)
            % 如果是最外层的连通域，舍弃
            if i ~= rawIndex && isempty(find(child == i, 1)) 
                B = [B; BTmp{i}];
            end
        end
    end
    
    % 计算连通域面积
    areaBoundary = [];
    for i = 1 : length(B)
        areaBoundary = [areaBoundary; regionprops(B{i} > 0, 'Perimeter').Perimeter];
    end
end

