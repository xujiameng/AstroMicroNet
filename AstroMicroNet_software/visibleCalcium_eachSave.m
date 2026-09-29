function visibleCalcium_eachSave(calcium_info, saveDir, baseName)
% visibleCalcium_eachSave  Plot each signal separately and save as vector.
%
% Usage:
%   visibleCalcium_eachSave(calcium_info)
%   visibleCalcium_eachSave(calcium_info, saveDir)
%   visibleCalcium_eachSave(calcium_info, saveDir, baseName)

    if nargin < 2 || isempty(saveDir)
        saveDir = pwd;
    end
    if nargin < 3 || isempty(baseName)
        baseName = 'calcium';
    end
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    n = numel(calcium_info);
    ts = datestr(now,'yyyymmdd_HHMMSS_FFF');  % 同一批导出的时间戳一致

    for i = 1:n

        % ===== 1) 单独画一条信号 =====
        fig = figure('Color','w','Visible','off'); % Visible='off' 不弹窗，直接保存
        ax = axes(fig); hold(ax,'on');

        dfSignal = calcium_info{i}.signal.dfSignal;
        start    = calcium_info{i}.calciumIndex.start;
        finish   = calcium_info{i}.calciumIndex.finish;

        plot(ax, dfSignal, 'k', 'LineWidth', 2);

        for j = 1:numel(start)
            idx = start(j,1):finish(j,1);
            plot(ax, idx, dfSignal(idx,1), 'r', 'LineWidth', 2);
        end

        axis(ax,'tight');
        box(ax,'off');
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        ax.Color = 'none';

        drawnow;

        % ===== 2) 保存为矢量 =====
        set(fig, 'Renderer', 'painters');  % 关键：矢量渲染器

        outBase = fullfile(saveDir, sprintf('%s_%s_%03d', baseName, ts, i));

%         print(fig, [outBase '.eps'], '-depsc', '-vector');
        print(fig, [outBase '.svg'], '-dsvg',  '-vector');

        % （可选）同时保存 PDF（AI 最稳）
        % print(fig, [outBase '.pdf'], '-dpdf', '-painters', '-bestfit');

        % ===== 3) 关闭该 figure（避免开太多）=====
        close(fig);
    end
end