function [distance_group] = distance_ROI(boundary)
% 计算距离
for i = length(boundary)
    ROI_center = Centroid_center(boundary);
    distance_group = zeros(length(ROI_center(1,:)));     % 存储两两之间的距离
end


x = 0;%x为表示计算的点是第x个
for i = 1 : length(ROI_center(:,1))   %i设定为被计算的点的位数
    x = x + 1;
    y = 0;%y为表示算与被计算点相关参数的点是第y个
    for m = 1 : length(ROI_center(:,1))%m为计算其与被计算点相关参数的点，进行遍历
        y = y + 1;
        % 距离
        dis = sqrt((ROI_center(x, 1) - ROI_center(y, 1))^2 +((ROI_center(x, 2) - ROI_center(y, 2))^2));% sqrt是平方根函数，计算距离
        distance_group(x, y) = dis;

    end
end

end



function [ROI_center] = Centroid_center(boundary)
    % 获取每个连通域的质心
    ROI_center = [];
    for i = 1:length(boundary(:,1))
        ROI_loc = boundary{i,1};
        center(1,1) = mean(ROI_loc(:,1));
        center(1,2) = mean(ROI_loc(:,2));
        ROI_center = [ROI_center;center];
    end
end

