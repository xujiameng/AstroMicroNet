function waveParams =  WaveParams(waveData, startIdx, endIdx, fs)
    % 输入：
    %   waveData：原始波形数据
    %   startIdx：所有波形的起始索引
    %   endIdx：所有波形的结束索引
    %   fs：采样率
    % 输出：
    %   waveParams：波形参数结构体数组
    %   featureMatrix：波形参数结构体数组

    if ~isempty(startIdx)

        waveData = waveData(:);  % 转为列向量
        if length(startIdx) ~= length(endIdx)
            error('startIdx 和 endIdx 必须长度相同');
        end
        numWaves = length(startIdx);
        %     waveParams_std(numWaves) = struct();  % 初始化结构体数组
        %     waveParams(numWaves) = struct();
        % 循环计算每个波形的参数
        for i = 1:numWaves
            % 1. 获取当前波形的起止索引
            s = startIdx(i);
            e = endIdx(i);

            % 检查索引有效性
            if s < 1 || e > length(waveData) || s >= e
                warning(['第', num2str(i), '个波形的起止索引无效，跳过计算']);
                continue;
            end

            % 2. 截取当前波形片段
            wave = waveData(s:e);  % 截取波形（包含起止点）
            t = (0:length(wave)-1) / fs;  % 时间轴（单位：秒）
            n = length(wave);

            % 3. 计算参数
            % （1）AMP（幅度：峰峰值）
            maxVal = max(wave);
            minVal = min(wave);
            AMP = maxVal - minVal;

            % （2）Area（面积：时域积分）
            Area = trapz(t, wave);

            % （3）Dura（持续时间）
            Dura = (e - s) / fs;  % 时间 = 索引差 / 采样率

            % （4）主峰检测（阈值待确定，0.1？）
            peakThreshold = 0.1 * AMP;  % 排除噪声峰的阈值
            if peakThreshold == 0  % 波形无波动
                peaks = [];
                mainPeakIdx = [];
                mainPeakVal = 0;
            else
                [peaks, peaks_loc] = findpeaks(wave, 'MinPeakHeight', peakThreshold);
                if isempty(peaks)
                    mainPeakIdx = [];
                    mainPeakVal = 0;
                else
                    % 取最高的峰为主峰
                    %                 [~, maxPeakIdx] = max(peakProps.PeakHeight);
                    %                 mainPeakIdx = peaks(maxPeakIdx);  % 主峰在当前波形中的相对索引
                    %                 mainPeakVal = wave(mainPeakIdx);
                    peakHeights = wave(peaks_loc);
                    [~, maxPeakIdx] = max(peakHeights);
                    mainPeakIdx = peaks_loc(maxPeakIdx);
                    mainPeakVal = peakHeights(maxPeakIdx);
                end
            end

            % （5）HW（半高宽）
            if isempty(mainPeakIdx) || mainPeakVal == 0
                HW = 0;
            else


                half = min(wave) + (wave(mainPeakIdx) - min(wave)) / 2;
                tmp_left = find_intersection(wave(1:mainPeakIdx), half);
                tmp_right = find_intersection(wave(mainPeakIdx:end), half);
                riseIdx = tmp_left(end);
                fallIdx = tmp_right(1)+mainPeakIdx-1;


                %             halfHeight = mainPeakVal / 2;
                %             % 上升沿：最后一个达到半高的点（主峰左侧）
                %             risePart = wave(1:mainPeakIdx);
                %             riseIdx = find(risePart >= halfHeight, 1, 'last');
                %             if isempty(riseIdx), riseIdx = 1; end  % 用起点
                %
                %             % 下降沿：第一个低于半高的点（主峰右侧）
                %             fallPart = wave(mainPeakIdx:end);
                %             fallIdxRel = find(fallPart <= halfHeight, 1, 'first');
                %             if isempty(fallIdxRel)
                %                 fallIdx = n;  % 用终点
                %             else
                %                 fallIdx = mainPeakIdx + fallIdxRel - 1;  % 转换为相对索引
                %             end
                HW = (fallIdx - riseIdx)/ fs;
            end

            % （6）谷值检测
            %        [valleys, peaks_loc] = findpeaks(-wave, 'MinPeakHeight', peakThreshold);  % 谷值是负向峰值,方法错误
            valleys_loc = find_intersection(wave, min(wave));
            valleyVals = wave(valleys_loc);  % 谷值实际值

            % （7）US（上升斜率）、（8）UT（上升时间）
            if isempty(mainPeakIdx) || mainPeakVal == 0
                US = 0;
                UT = 0;
            else
                % 左侧最近的谷值
                leftValleys = valleys_loc(valleys_loc < mainPeakIdx);
                if isempty(leftValleys)
                    leftValleyIdx = 1;  % 用起点
                    leftValleyVal = wave(leftValleyIdx);
                else
                    leftValleyIdx = leftValleys(end);  % 最近左侧谷值
                    leftValleyVal = wave(leftValleyIdx);
                end
                % 计算斜率和时间
                timeDiff = t(mainPeakIdx) - t(leftValleyIdx);
                if timeDiff == 0
                    US = 0;
                    UT = 0;
                else
                    US = (mainPeakVal - leftValleyVal) / timeDiff;
                    UT = timeDiff;
                end
            end

            % （9）DS（下降斜率）、（10）DT（下降时间）
            if isempty(mainPeakIdx) || mainPeakVal == 0
                DS = 0;
                DT = 0;
            else
                % 右侧最近的谷值
                rightValleys = valleys_loc(valleys_loc > mainPeakIdx);
                if isempty(rightValleys)
                    rightValleyIdx = n;  % 用终点
                    rightValleyVal = wave(rightValleyIdx);
                else
                    rightValleyIdx = rightValleys(1);  % 最近右侧谷值
                    rightValleyVal = wave(rightValleyIdx);
                end
                % 计算斜率和时间
                timeDiff = t(rightValleyIdx) - t(mainPeakIdx);
                if timeDiff == 0
                    DS = 0;
                    DT = 0;
                else
                    DS = (rightValleyVal - mainPeakVal) / timeDiff;
                    DT = timeDiff;
                end
            end

            % （11）NP（主峰数量）
            NP = length(peaks);

            % （12）NV（谷值数量）
            NV = length(valleys_loc);

            % （13）Bump（谷值为峰值一半的数量）
            if isempty(valleys_loc) || mainPeakVal == 0
                Bump = 0;
            else
                Bump = sum(valleyVals >= (mainPeakVal+min(wave) / 2));
            end

            % （14）MPDura（主峰持续时间，等于半高宽）
            MPDura = HW;

            % 存储参数
            %          AMP_all = [AMP_all;AMP];
            %          Area_all = [Area_all;Area];
            %          Dura_all = [Dura_all;Dura];
            %          HW_all = [HW_all;HW];
            %          US_all = [US_all;US];
            %          DS_all = [DS_all;DS];
            %          UT_all = [UT_all;UT];
            %          DT_all = [DT_all;DT];
            %          NP_all = [NP_all;NP];
            %          NV_all = [NV_all;NV];
            %          Bump_all = [Bump_all;Bump];
            %          MPDura_all = [MPDura_all;MPDura];

            waveParams.AMP(i) = AMP;
            waveParams.Area(i)  = Area;
            waveParams.Dura(i)  =Dura;
            waveParams.HW(i)  = HW;
            waveParams.US(i)  = US;
            waveParams.DS(i)  = DS;
            waveParams.UT(i)  = UT;
            waveParams.DT(i)  =DT;
            waveParams.NP(i)  =NP;
            waveParams.NV(i)  = NV;
            waveParams.Bump(i)  = Bump;
            waveParams.MPDura(i)  = MPDura;
        end

    else 
        
            waveParams.AMP = [];
            waveParams.Area  = [];
            waveParams.Dura  =[];
            waveParams.HW  = [];
            waveParams.US  = [];
            waveParams.DS  = [];
            waveParams.UT  = [];
            waveParams.DT  =[];
            waveParams.NP  =[];
            waveParams.NV  = [];
            waveParams.Bump  = [];
            waveParams.MPDura  = [];
        

    end
%         featureMatrix(:, 1) = waveParams.AMP;
%         featureMatrix(:, 2) = waveParams.Area;
%         featureMatrix(:, 3) = waveParams.Dura;
%         featureMatrix(:, 4) = waveParams.HW;
%         featureMatrix(:, 5) = waveParams.US;
%         featureMatrix(:, 6) = waveParams.DS;
%         featureMatrix(:, 7) = waveParams.UT;
%         featureMatrix(:, 8) = waveParams.DT;
%         featureMatrix(:, 9) = waveParams.NP;
%         featureMatrix(:, 10) = waveParams.NV;
%         featureMatrix(:, 11) = waveParams.Bump;
%         featureMatrix(:, 12) = waveParams.MPDura;
end