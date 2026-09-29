function [bestCore,coreIndices,uncoreIndices] = coredetect_heredity(adjMatrix,distance_soma_micro,popSize, maxGen, crossoverRate, mutationRate,dis_therehold,core_precent)
% adjMatrix：输入的参数矩阵
% distance_soma_micro:ROI距离矩阵
% popSize：种群大小，每代中包含的个体数量
% maxGen：递归次数
% crossoverRate：交叉概率，控制交叉操作的发生频率
% mutationRate：变异概率，控制变异操作的发生频率
% dis_therehold:矩阵的阈值设置
% core_precent：核心比例
Thresh =  dis_therehold;
N = size(adjMatrix, 1);
for i = 1:N
    for j = i+1:N
        if adjMatrix(i,j) < Thresh

            adjMatrix(i,j) = 0;

        end
    end
end

save_precent = 20;

% 遗传算法主函数
N = size(adjMatrix, 1);


[~, dissortedIndices] = sort(distance_soma_micro(1,:), 'descend');
% 初始化种群（随机生成改为按距离排序）
population = rand(popSize, N);
for i = 1:save_precent
    para_orgin = population(i,:);
    for j = 1:N
        population(i,j) = para_orgin(dissortedIndices(j));
    end
end





% 迭代进化
for gen = 1:maxGen
    % 计算适应度
    fitness = zeros(popSize, 1);
    for i = 1:popSize
        fitness(i) = coreFitness(adjMatrix,population(i,:)');
    end

    [~, sortedIndices] = sort(fitness, 'descend');%对计算的核心总量排序


    for i = 1:popSize

        newPopulation(i,:) = population(sortedIndices(i),:);
        fitness_sort(i,:) = fitness(sortedIndices(i),:);

    end

    % 交叉操作
    for i = save_precent+1:2:popSize-1
        if rand < crossoverRate
            % 单点交叉
            crossPoint = randi(N-1);
            crosslen = randi(N-crossPoint);
            temp = newPopulation(i, crossPoint+1:crossPoint+crosslen);
            newPopulation(i, crossPoint+1:crossPoint+crosslen) = newPopulation(i+1, crossPoint+1:crossPoint+crosslen);
            newPopulation(i+1, crossPoint+1:crossPoint+crosslen) = temp;
        end
    end

    % 变异操作
    for i = save_precent+1:popSize
        for j = 1:N
            if rand < mutationRate
                newPopulation(i,j) = rand;  % 随机变异
            end
        end
    end

    population = newPopulation;

    % 显示进度
    if mod(gen, 10000) == 0
        fprintf('遗传算法迭代 %d/%d，最佳适应度: %.4f\n', gen, maxGen, max(fitness));
    end
end

% 返回最优解
[~, bestIdx] = max(fitness);
bestCore = population(bestIdx, :)';

%计算核心外围分布
[~, sortedIndices] = sort(bestCore, 'descend');
coreIndices = sortedIndices(1:floor(core_precent*length(bestCore)));
uncoreIndices = sortedIndices(floor(core_precent*length(bestCore))+1:end);



end



%% 计算核心质量Rγ
function rScore = coreFitness(adjMatrix, coreVector)
    % 实现Rγ = sum_{i,j} A_ij * C_i * C_j
    rScore = adjMatrix * coreVector;
    rScore = coreVector' * rScore;
end