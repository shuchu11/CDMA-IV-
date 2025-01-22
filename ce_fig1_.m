% 讀取數據
data = readtable('ce_waveform_data.csv');

% 提取數據
iterations = unique(data.Iteration);
colors = lines(length(iterations)); % 為不同迭代次數生成顏色

figure;
hold on;
for i = 1:length(iterations)
    iter = iterations(i);
    iter_data = data(data.Iteration == iter, :);
    % 添加 LineWidth 屬性來調寬軌跡線
    plot(iter_data.Time, iter_data.Phase, 'Color', colors(i, :), ...
         'LineWidth', 2, 'DisplayName', sprintf('Iteration %d', iter));
end

% 圖例與標題
legend('show');
xlabel('Time [T_c]');
ylabel('Phase [radians]');
title('Phase Evolution of Signature Waveforms');
grid on;
hold off;
