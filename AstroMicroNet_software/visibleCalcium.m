function visibleCalcium(calcium_info)
     figure();
%     fig = figure();
    for i = 1 : length(calcium_info)
        subplot(length(calcium_info), 1, i)
        
        dfSignal = calcium_info{i}.signal.dfSignal;
        start = calcium_info{i}.calciumIndex.start;
        finish = calcium_info{i}.calciumIndex.finish;
        plot(dfSignal, 'k', 'LineWidth', 2); hold on;
        for j = 1 : length(start)
            plot(start(j, 1) : finish(j, 1), ...
                dfSignal(start(j, 1) : finish(j, 1), 1),'r', 'LineWidth', 2); hold on;   
        end

        %È¥³ýÉÏÓÒ±ß¿ò¿Ì¶È
        box off  
        %ÒÆ³ý×ø±êÖá±ß¿ò
        set(gca,'Visible','off');
        %ÉèÖÃ±³¾°Îª°×É«
        set(gcf,'color','w');
    end
end

