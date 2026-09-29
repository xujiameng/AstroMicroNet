function [bestCore, coreIndices, uncoreIndices] = coredetect_particle(adjMatrix,distance_matrix, numParticles, maxIter,c1, c2, dis_therehold,core_precent)
% 输入参数:
% adjMatrix - 邻接矩阵
% numParticles - 粒子数量
% maxIter - 最大迭代次数
% w - 惯性权重
% c1 - 个体学习因子
% c2 - 社会学习因子
% core_precent - 核心节点比例

N = size(adjMatrix, 1);
core_size = floor(N * core_precent);  % 核心节点数量
for i = 1:N
    for j = i+1:N
        if adjMatrix(i,j) < dis_therehold

            adjMatrix(i,j) = 0;

        end
    end
end

[~, dissortedIndices] = sort(adjMatrix(1,:), 'descend');
gbest_fitness_all = 0;
gbest_position_all = [];



for p = 1:200
%     fprintf('第%d次迭代开始: \n',p);

    % 初始化粒子位置和速度
    particles_1 = rand( numParticles-core_size,N);  % 位置在[0,1]之间
    particles_2 = zeros(numParticles, N);
    for i = 1:numParticles
        for j = 1:i
            particles_2(i,dissortedIndices(j)) = 1;
        end
    end
    particles = [particles_1;particles_2];

%     rng(123);  % 设置固定种子，确保每次初始化相同
%     particles = rand(numParticles, N);
    velocities = 0.1 * randn(numParticles, N);  % 初始速度


    % 初始化个体最优和全局最优
    pbest_position = particles;  % 个体历史最优位置
    pbest_fitness = zeros(numParticles, 1);  % 个体历史最优适应度

    gbest_fitness = -inf;  % 全局最优适应度
    gbest_position = zeros(1, N);  % 全局最优位置

    % 迭代优化
    for iter = 1:maxIter
        % 动态调整惯性权重（线性递减）
        w = 0.9 - 0.5 * (iter / maxIter);

        % 评估每个粒子
        for i = 1:numParticles
            % 将粒子位置转换为核心向量（取概率最高的节点）
            [~, sorted_idx] = sort(particles(i,:), 'descend');
            core_vector = zeros(N, 1);
            core_vector(sorted_idx(1:core_size)) = 1;

            % 计算适应度
            current_fitness = coreFitness(adjMatrix, core_vector,distance_matrix);

            % 更新个体最优
            if current_fitness > pbest_fitness(i)
                pbest_fitness(i) = current_fitness;
                pbest_position(i,:) = particles(i,:);
            end

            % 更新全局最优
            if current_fitness > gbest_fitness
                gbest_fitness = current_fitness;
                gbest_position = particles(i,:);
            end
        end

        % 更新所有粒子的速度和位置
        for i = 1:numParticles
            r1 = rand();
            r2 = rand();

            % 速度更新公式
            velocities(i,:) = w * velocities(i,:) ...
                + c1 * r1 * (pbest_position(i,:) - particles(i,:)) ...
                + c2 * r2 * (gbest_position - particles(i,:));

            % 位置更新
            particles(i,:) = particles(i,:) + velocities(i,:);

            % 边界处理（确保位置在[0,1]范围内）
            particles(i,:) = max(0, min(1, particles(i,:)));
        end

        % 显示进度
%                 if mod(iter, 100) == 0
%                     fprintf('迭代 %d/%d，最佳适应度: %.4f\n', iter, maxIter, gbest_fitness);
%                 end

    end
%     fprintf('最佳适应度: %.4f\n', gbest_fitness);
    if mod(p, 50) == 0
        fprintf('迭代 %d次，\n',p);
    end
    %全局归纳 
    if gbest_fitness > gbest_fitness_all
        gbest_fitness_all = gbest_fitness;
        gbest_position_all = gbest_position;
    end

end

fprintf('多次迭代后，最佳适应度: %.4f\n', gbest_fitness_all);


% 输出最终结果
% [~, sorted_idx] = sort(gbest_position, 'descend');
% coreIndices = sorted_idx(1:core_size);
% uncoreIndices = sorted_idx(core_size+1:end);
[~, sorted_idx] = sort(gbest_position_all, 'descend');
coreIndices = sorted_idx(1:core_size);
uncoreIndices = sorted_idx(core_size+1:end);

bestCore = zeros(N, 1);
bestCore(coreIndices) = 1;
end


%% 计算核心质量Rγ
function rScore = coreFitness(adjMatrix, coreVector,distance_matrix)
% 实现Rγ = sum_{i,j} A_ij * C_i * C_j
rScore = adjMatrix * coreVector;
rScore = coreVector' * rScore;
rScore = rScore/sum(coreVector);

distance_matrix = distance_matrix/max(distance_matrix(:));
% 计算核心ROI之间的平均距离
core_indices = find(coreVector);

% 提取核心ROI之间的距离
core_distances = distance_matrix(core_indices, core_indices);
% 计算上三角部分的平均值（排除对角线）
triu_indices = triu(true(length(core_indices)), 1);
avg_distance = mean(core_distances(triu_indices));


rScore = rScore  + avg_distance;
end





%% 计算核心质量Rγ（初始）
% function rScore = coreFitness(adjMatrix, coreVector)
% % 实现Rγ = sum_{i,j} A_ij * C_i * C_j
% rScore = adjMatrix * coreVector;
% rScore = coreVector' * rScore;
% rScore = rScore/sum(coreVector);
% 
% end