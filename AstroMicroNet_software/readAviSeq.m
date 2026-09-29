function [data,avi_stream] = readAviSeq(fileNmae)
    avi_stream = VideoReader(fileNmae); %读入
    height = avi_stream.Height; % 每帧高
    width = avi_stream.Width;   % 每帧宽
%     frame_rate = avi_stream.FrameRate;  % 帧率
    frame_num = avi_stream.NumFrames;  % 帧数
    
    data = zeros(height, width, frame_num, 'uint8');
    for i = 1 : frame_num
        if mod(i, 100) == 0
            fprintf("%d/%d Frames\n", i, frame_num)
        end
        data(:, :, i) = readFrame(avi_stream);
    end
end

