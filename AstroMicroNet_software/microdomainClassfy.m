function [save_microdomain_struct, all_mask] = microdomainClassfy(mask, boundary, boundary_name)
    % 根据连通域和掩膜，获取掩膜内的所有ROI
    % mask:圈定的掩膜
    % boundary：微域的边界
    
    [m, n] = size(mask);
    save_microdomain_struct = {};  % 保留的ROI信息，结构体包括mask和boundary
    all_mask = zeros(m, n);
    for i = 1 : length(boundary)
        tmp_name = boundary_name{i, 1};
        if tmp_name ~= "microdomain"
            continue;
        end
        index = sub2ind([m, n], boundary{i, 1}(:, 1), boundary{i, 1}(:, 2));
        tmp_ROI = zeros(m, n);
        tmp_ROI(index) = 1;
        [~,L,~,~] = bwboundaries(tmp_ROI,8); % 8邻域找连通域
        tmp_mask = L .* mask;
        % 如果该连通域与掩膜有交集 就保留
        if ~isempty(find(tmp_mask == 1, 1))
            tmp.ROIInfo.mask = L;
            tmp.ROIInfo.boundary = boundary{i, 1};
            all_mask = all_mask + L;
            if isempty(save_microdomain_struct)
                save_microdomain_struct = {tmp};
            else
                save_microdomain_struct = [save_microdomain_struct; tmp];
            end
        end
    end
end

