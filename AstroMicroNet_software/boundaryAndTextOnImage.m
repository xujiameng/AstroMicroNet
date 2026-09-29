function boundaryAndTextOnImage(B,L,set_flag, structName)
%输入：B：范围坐标合集,L：标志位,set_flag:绘制相同颜色的核心外围时使用 structName：结构体名称
%输出：将B内的坐标集按不同颜色绘制出来

%颜色选集
%     colors=['b' 'g' 'r' 'c' 'm' 'y'];
%     colors = ["#0067C1", "#87CEFA", "#98FB98", "#CFD62F", "#DA67E1", "#00FFFF"]"#FF5733", "#33FF57";
colors = ["#CA1E1D", "#1672B1", "#3357FF", "#F3FF33", "#FF33F3", "#ADD8E6", "#DDA0DD",...
          "#FF9933", "#99FF33", "#33FF99", "#3399FF", "#9933FF", "#FF3399", "#99FF99", "#3333FF",...
          "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFFFFF",...
          "#FFA500", "#800080", "#008000", "#000080", "#800000", "#808000", "#008080",...
          "#C0C0C0", "#808080", "#FFC0CB", "#90EE90", "#ADD8E6", "#F0E68C", "#20B2AA", "#CD5C5C"];


if set_flag == 0
    if ~isempty(structName)

        %  不显示结构名称
        if L == 0
            for k=1 : length(B)
                if structName{k, 2} == 1
                    boundary = B{k};
                    cidx = mod(k,length(colors))+1;%选择颜色
                    rgb_color = hex2rgb(colors(cidx));
                    plot(boundary(:,2), boundary(:,1),'Color',...
                        rgb_color,'LineWidth',2);
                    %              plot(boundary(:,2), boundary(:,1),...
                    %                colors(cidx),'LineWidth',2);
                end
            end
        end
        %显示结构名称
        if L ~= 0
            for k=1 : length(B)
                %         if structName{k, 1} == "soma"||structName{k, 1} == "Main Branch"||structName{k, 1} == "End Feet"
                boundary = B{k};
                cidx = mod(k,length(colors))+1;%选择颜色
                rgb_color = hex2rgb(colors(cidx));
                plot(boundary(:,2), boundary(:,1),'Color',...
                    rgb_color,'LineWidth',2);

                %randomize text position for better visibility（在列表勾选不使用）
                rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
                col = boundary(rndRow,2); row = boundary(rndRow,1);
                %h = text(col+1, row-1, num2str(L(row,col)));
                h = text(col+1, row-1, num2str(structName{k, 1}));
                set(h,'Color',rgb_color,'FontSize',14,'FontWeight','bold');
                %             set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
                %         end
            end
        end

    else
        %忽略判断，圈定层级范围使用
        for k=1 : length(B)

            boundary = B{k};

            rgb_color = hex2rgb("#FF0000");
            plot(boundary(:,2), boundary(:,1),'Color',...
                rgb_color,'LineWidth',2);

        end
    end
else
    if L ==0
        for k=1 : length(B)
            if structName{k, 2} == 1
                boundary = B{k};

                rgb_color = hex2rgb(colors(set_flag));
                plot(boundary(:,2), boundary(:,1),'Color',...
                    rgb_color,'LineWidth',2);

            end
        end
    else
        for k=1 : length(B)
            if structName{k, 2} == 1
                boundary = B{k};

                rgb_color = hex2rgb(colors(set_flag));
                fill(boundary(:,2), boundary(:,1), rgb_color)
                plot(boundary(:,2), boundary(:,1),'Color',...
                    rgb_color,'LineWidth',2);

            end
        end
    end
end
end

function rgb = hex2rgb(hex)
    %  hex 是字符串去掉 #
    if startsWith(hex, '#')
        hex = char(hex);
        hex = hex(2:end); % 移除 #
    end
    rgb = sscanf(hex, '%2x%2x%2x')' / 255;
end