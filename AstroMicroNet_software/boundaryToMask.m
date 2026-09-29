function [mask] = boundaryToMask(data, boundary)
    % 根据连通域获取掩膜
    [m, n] = size(data);
    mask = zeros(m, n, length(boundary));
    
    for i = 1 : length(boundary)
        tmp_mask = zeros(m, n);

        index = sub2ind([m, n], boundary{i, 1}(:, 1), boundary{i, 1}(:, 2));
        tmp_mask(index) = 1;
        [~,L,~,~] = bwboundaries(tmp_mask,8); % 8邻域找连通域
        mask(:, :, i) = L;
    end
end

