%% 示例：使用合成数据测试算法
function [coreScores, adjMatrix, coreIndices,uncoreIndices,corepoints,uncorepoints] = coredetect(boundary_center, distance_group,correlation_group,dis_therehold,core_precent)
    % boundary_center: ROI质心
    % distance_group: 距离矩阵
    % coreScores: 聚合核心分数
    % adjMatrix: 构建核心网络
    % corepoints: 核心点集
    % uncorepoints: 非核心点集
    %core_precent: 核心占比
    
    % 构建距离网络
%     adjMatrix = distance_group/max(distance_group(:));
    adjMatrix = buildNetwork(boundary_center,distance_group,correlation_group,dis_therehold,...
         'density', 'correlation');
    % 设置参数
    paramSettings.alphaRange = 0.1:0.1:1;
    paramSettings.betaRange = 0.1:0.1:1;
    paramSettings.saIterations = 500;
    paramSettings.saInitialTemp = 1;           
    paramSettings.saCoolingRate = 0.8;           
    paramSettings.saFinalTemp = 1e-8;            


    % 执行核心-外围检测
    [coreScores, ~, ~] = corePeripheryDetection(adjMatrix, paramSettings);
    
    % 可视化结果
    [coreIndices,uncoreIndices,corepoints,uncorepoints] = visualizeResults(boundary_center, coreScores,core_precent);
end
%% 核心-外围结构检测算法实现
function [coreScores, allCoreVectors, allRscores] = corePeripheryDetection(adjMatrix, paramSettings)
    % adjMatrix: 输入邻接矩阵
    % paramSettings: 参数设置结构体
    % coreScores: 聚合核心分数
    % allCoreVectors: 各参数下的核心度向量
    % allRscores: 各参数下的核心质量值
    
    % 初始化参数
    if nargin < 2
        paramSettings.alphaRange = 0.1:0.1:0.5;     % α参数范围
        paramSettings.betaRange = 0.2:0.1:1;      % β参数范围
        paramSettings.saIterations = 1000;           % 模拟退火迭代次数
        paramSettings.saInitialTemp = 1;             % 初始温度
        paramSettings.saCoolingRate = 0.8;           % 冷却率
        paramSettings.saFinalTemp = 1e-8;            % 最终温度
    end
    
    N = size(adjMatrix, 1);                        % 节点数
    alphaVals = paramSettings.alphaRange;
    betaVals = paramSettings.betaRange;
    numParams = length(alphaVals) * length(betaVals);
    
    % 预分配存储变量
    allCoreVectors = zeros(N, numParams);
    allRscores = zeros(1, numParams);
    paramIndex = 0;
    
    % 遍历所有α和β参数组合
    for a = 1:length(alphaVals)
        alpha = alphaVals(a);
        for b = 1:length(betaVals)
            beta = betaVals(b);
            paramIndex = paramIndex + 1;
            
            % 1. 生成过渡函数向量
            coreVector = generateTransitionVector(N, alpha, beta);
            
            % 2. 使用模拟退火优化核心质量
            [optimizedCore, rScore] = simulateAnnealing(adjMatrix, coreVector, ...
                paramSettings.saIterations, paramSettings.saInitialTemp, ...
                paramSettings.saCoolingRate, paramSettings.saFinalTemp);
            
            % 存储结果
            allCoreVectors(:, paramIndex) = optimizedCore;
            allRscores(paramIndex) = rScore;
            
            % 显示进度
            fprintf('完成参数组合 %d/%d (α=%.2f, β=%.2f)\n', paramIndex, numParams, alpha, beta);
        end
    end
    
    % 3. 计算聚合核心分数
    coreScores = aggregateCoreScores(allCoreVectors, allRscores);
end

%% 生成过渡函数向量
function transitionVec = generateTransitionVector(N, alpha, beta)
    % 实现尖锐双参数过渡函数
    betaN = floor(beta * N);
    transitionVec = zeros(N, 1);
    
    if betaN == 0  % β=0时，所有节点按外围处理
        for i = 1:N
            transitionVec(i) = ((i) * (1-alpha)) / (2*(N - betaN)) + (1+alpha)/2;
        end
    elseif betaN == N  % β=1时，所有节点按核心处理
        for i = 1:N
            transitionVec(i) = (i * (1-alpha)) / (2*betaN);
        end
    else
        % 分两部分生成过渡向量
        for i = 1:betaN
            transitionVec(i) = (i * (1-alpha)) / (2*betaN);
        end
        for i = betaN+1:N
            transitionVec(i) = ((i - betaN) * (1-alpha)) / (2*(N - betaN)) + (1+alpha)/2;
        end
    end
end

%% 模拟退火优化核心质量
function [optimizedCore, rScore] = simulateAnnealing(adjMatrix, initialCore, maxIter, initialTemp, coolingRate, finalTemp)
    % 初始化
    currentCore = initialCore;
    currentTemp = initialTemp;
    bestCore = currentCore;
    bestR = computeCoreQuality(adjMatrix, currentCore);
    
    % 模拟退火主循环
    for iter = 1:maxIter
        % 生成新解：随机交换两个元素
        newCore = currentCore;
        i = randi(size(newCore, 1));
        j = randi(size(newCore, 1));
        temp = newCore(i);
        newCore(i) = newCore(j);
        newCore(j) = temp;
        
        % 计算新旧核心质量
        currentR = computeCoreQuality(adjMatrix, currentCore);
        newR = computeCoreQuality(adjMatrix, newCore);
        
        % 接受准则
        if newR > currentR || rand() < exp((newR - currentR) / currentTemp)
            currentCore = newCore;
            if newR > bestR
                bestCore = newCore;
                bestR = newR;
            end
        end
        
        % 降温
        currentTemp = currentTemp * coolingRate;
        if currentTemp < finalTemp
            break;
        end
    end
    
    optimizedCore = bestCore;
    rScore = bestR;
end

%% 计算核心质量Rγ
function rScore = computeCoreQuality(adjMatrix, coreVector)
    % 实现Rγ = sum_{i,j} A_ij * C_i * C_j
    rScore = adjMatrix * coreVector;
    rScore = coreVector' * rScore;
end

%% 聚合核心分数
function coreScores = aggregateCoreScores(allCoreVectors, allRscores)
    % 计算CS(i) = Z * sum_{γ} C_i(γ) * Rγ
    N = size(allCoreVectors, 1);
    numParams = size(allCoreVectors, 2);
    coreScores = zeros(N, 1);
    
    % 加权求和
    for i = 1:N
        for p = 1:numParams
            coreScores(i) = coreScores(i) + allCoreVectors(i, p) * allRscores(p);
        end
    end
    
    % 归一化
    maxScore = max(coreScores);
    if maxScore > 0
        coreScores = coreScores / maxScore;
    end
end


%% 构建示例网络（以距离矩阵为例）
function adjMatrix = buildNetwork(points, distance,correlation, dis_therehold,weightType,paratype)
    % points: N×2坐标矩阵
    % distanceThresh: 距离阈值
    % weightType: 权重类型 ('binary'或'gaussian')
    % 后续可否改为相关性等计算网络矩阵
    N = size(points, 1);
    adjMatrix = zeros(N, N);
    correlation_used = correlation;
    if strcmp(paratype, 'distance')
       Thresh =  max(distance(:))/2;
        for i = 1:N
            for j = i+1:N
                dist = distance(i,j);%距离改为相关性？
                if dist >= Thresh
                    if strcmp(weightType, 'binary')
                        adjMatrix(i,j) = 1;
                        adjMatrix(j,i) = 1;
                    elseif strcmp(weightType, 'gaussian')
                        adjMatrix(i,j) = exp(-dist^2 / (2 * Thresh^2));
                        adjMatrix(j,i) = adjMatrix(i,j);
                    elseif strcmp(weightType, 'density')
                        % 基于局部密度的连接权重(测试）
                        localDensity = sum(exp(-distance(distance <= Thresh).^2 / (2 * Thresh^2)));
                        adjMatrix(i,j) = localDensity;
                        adjMatrix(j,i) = localDensity;
                    end
                end
            end
        end
    elseif strcmp(paratype, 'correlation')
        Thresh =  dis_therehold;
%         Thresh =  0;
%         Thresh =  max(correlation(:))/1.5;
%          correlation_used(correlation < dis_therehold) = 0;
        for i = 1:N
            for j = i+1:N
                dist = correlation_used(i,j);
                if dist >= Thresh
                    if strcmp(weightType, 'binary')
                        adjMatrix(i,j) = 1;
                        adjMatrix(j,i) = 1;
                    elseif strcmp(weightType, 'gaussian')
                        adjMatrix(i,j) = exp(-dist^2 / (2 * Thresh^2));
                        adjMatrix(j,i) = adjMatrix(i,j);
                    elseif strcmp(weightType, 'density')
                        % 基于局部密度的连接权重(测试）
%                         localDensity = sum(exp(-correlation(correlation <= Thresh).^2 / (2 * Thresh^2)));
                        adjMatrix(i,j) = correlation_used(i,j);
%                         adjMatrix(j,i) = correlation(j,i);
                    end
                end
            end
        end
    end
end

%% 可视化结果
function [coreIndices,uncoreIndices,corepoints,uncorepoints] = visualizeResults(points, coreScores,core_precent)
    % 按核心分数着色显示数据点
    figure
    scatter(points(:,1), points(:,2), 50, coreScores, 'filled');
    colorbar;
    title('核心-外围结构检测结果');
    xlabel('X坐标');
    ylabel('Y坐标');
    caxis([0, 1]);
    colormap(hot);

    % 标记核心点
    [~, sortedIndices] = sort(coreScores, 'descend');
    coreIndices = sortedIndices(1:floor(core_precent*length(coreScores)));
    uncoreIndices = sortedIndices(floor(core_precent*length(coreScores))+1:end);

    corepoints = points(coreIndices,:);
    uncorepoints = points(uncoreIndices,:);
    hold on;
    scatter(points(coreIndices,1), points(coreIndices,2), 100, 'r', 'x', 'LineWidth', 2);
    legend('所有点（核心分数）', '核心点');
end