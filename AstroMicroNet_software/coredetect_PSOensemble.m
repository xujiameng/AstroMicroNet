function [ensembleCore,finalCoreIndices, finalUncoreIndices] = coredetect_PSOensemble(adjMatrix, distance_matrix, dis_therehold, core_precent)
    % ============================================================
    % 多参数集成粒子群算法核心–外围结构检测
    %
    % 输入:
    %   adjMatrix       —— 邻接矩阵（微域相关性矩阵）
    %   distance_matrix —— ROI间距离矩阵
    %   dis_therehold   —— 连接阈值
    %   core_precent    —— 核心节点比例 (0~1)
    %
    % 输出:
    %   finalCoreIndices   —— 最终识别的核心节点索引
    %   finalUncoreIndices —— 最终识别的外围节点索引
    %   ensembleCore       —— 聚合后的核心分数向量
    %
    % 
    % ============================================================

    % 粒子群参数组合空间
    numParticles = [30, 50];
    maxIter = 150;
    c1_list = [1.5, 2.0];
    c2_list = [1.5, 2.0];
    w_list  = [0.7, 0.9]; % 惯性权重（内部动态更新，但初始可不同）

    numConfigs = length(numParticles)*length(c1_list)*length(c2_list)*length(w_list);
    N = size(adjMatrix,1);

    coreMatrix = zeros(N, numConfigs);
    scoreList  = zeros(1, numConfigs);

    configIndex = 0;

    fprintf('🚀 开始执行多参数粒子群算法集成优化...\n');

    for np = numParticles
        for c1 = c1_list
            for c2 = c2_list
                for w_init = w_list
                    configIndex = configIndex + 1;

                    % 调用单次PSO检测
                    [bestCore, coreIndices, uncoreIndices, fitnessVal] = ...
                        run_PSO_once(adjMatrix, distance_matrix, np, maxIter, c1, c2, w_init, dis_therehold, core_precent);

                    % 存储结果
                    coreMatrix(:, configIndex) = bestCore;
                    scoreList(configIndex) = fitnessVal;

                    fprintf('组合 %02d/%02d → np=%d | c1=%.1f | c2=%.1f | w=%.1f | Score=%.4f\n', ...
                        configIndex, numConfigs, np, c1, c2, w_init, fitnessVal);
                end
            end
        end
    end

    % ===== 聚合结果 =====
    scoreWeights = scoreList / sum(scoreList);
    ensembleCore = coreMatrix * scoreWeights';
    ensembleCore = ensembleCore / max(ensembleCore);

    % ===== 提取最终核心节点 =====
    [~, sortedIndices] = sort(ensembleCore, 'descend');
    finalCoreIndices = sortedIndices(1:floor(core_precent*N));
    finalUncoreIndices = sortedIndices(floor(core_precent*N)+1:end);

    fprintf('✅ 多参数PSO集成完成: 核心 %d | 外围 %d\n', ...
        length(finalCoreIndices), length(finalUncoreIndices));
end

%% ================== 子函数：单次PSO运行 ==================
function [bestCore, coreIndices, uncoreIndices, gbest_fitness] = ...
    run_PSO_once(adjMatrix, distance_matrix, numParticles, maxIter, c1, c2, w_init, dis_therehold, core_precent)

    N = size(adjMatrix, 1);
    core_size = floor(N * core_precent);

    % 阈值过滤
    for i = 1:N
        for j = i+1:N
            if adjMatrix(i,j) < dis_therehold
                adjMatrix(i,j) = 0;
            end
        end
    end

    % 初始化粒子位置和速度
    particles = rand(numParticles, N);
    velocities = 0.1 * randn(numParticles, N);

    % 初始化个体与全局最优
    pbest_position = particles;
    pbest_fitness = zeros(numParticles, 1);
    gbest_fitness = -inf;
    gbest_position = zeros(1, N);

    % 迭代优化
    for iter = 1:maxIter
        % 动态调整惯性权重（线性递减）
        w = w_init - 0.5 * (iter / maxIter);

        for i = 1:numParticles
            [~, sorted_idx] = sort(particles(i,:), 'descend');
            core_vector = zeros(N, 1);
            core_vector(sorted_idx(1:core_size)) = 1;

            current_fitness = coreFitness(adjMatrix, core_vector, distance_matrix);

            % 更新个体和全局最优
            if current_fitness > pbest_fitness(i)
                pbest_fitness(i) = current_fitness;
                pbest_position(i,:) = particles(i,:);
            end

            if current_fitness > gbest_fitness
                gbest_fitness = current_fitness;
                gbest_position = particles(i,:);
            end
        end

        % 更新速度与位置
        for i = 1:numParticles
            r1 = rand();
            r2 = rand();

            velocities(i,:) = w * velocities(i,:) ...
                + c1 * r1 * (pbest_position(i,:) - particles(i,:)) ...
                + c2 * r2 * (gbest_position - particles(i,:));

            particles(i,:) = particles(i,:) + velocities(i,:);
            particles(i,:) = max(0, min(1, particles(i,:))); % 保持在[0,1]
        end
    end

    % 输出最终结果
    [~, sorted_idx] = sort(gbest_position, 'descend');
    coreIndices = sorted_idx(1:core_size);
    uncoreIndices = sorted_idx(core_size+1:end);

    bestCore = zeros(N, 1);
    bestCore(coreIndices) = 1;
end

function rScore = coreFitness(adjMatrix, coreVector,distance_matrix)
% 实现Rγ = sum_{i,j} A_ij * C_i * C_j
rScore = adjMatrix * coreVector;
rScore = coreVector' * rScore;
rScore = rScore/sum(coreVector);

% distance_matrix = distance_matrix/max(distance_matrix(:));
% % 计算核心ROI之间的平均距离
% core_indices = find(coreVector);
% 
% % 提取核心ROI之间的距离
% core_distances = distance_matrix(core_indices, core_indices);
% % 计算上三角部分的平均值（排除对角线）
% triu_indices = triu(true(length(core_indices)), 1);
% avg_distance = mean(core_distances(triu_indices));
% 
% 
% rScore = rScore  + avg_distance;
end

