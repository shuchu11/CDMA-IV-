clear all;
% MATLAB 程式碼: 生成根升餘弦濾波器數據並輸出為文件
rf = 0.1;      % 滾降係數
span = 16  ;     % 濾波器範圍（符號數）
sps = 32;      % 每符號取樣數

% 計算根升餘弦濾波器
filter = rcosdesign(rf, span, sps, 'sqrt');

% 將濾波器係數存到一個檔案中
fileID = fopen('filter_data.txt', 'w');
fprintf(fileID, '%.8e\n', filter);
fclose(fileID);

% 確認數據輸出完成
disp('Filter data has been written to filter_data.txt');
