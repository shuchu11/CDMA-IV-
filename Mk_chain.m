% 馬可夫鏈模型
% P0 為初始概率向量，P 為轉移概率矩陣，chain_length 為馬可夫鏈長度

P = [0.67 0.33; 0.33 0.67]; % 轉移概率矩陣
P0 = [0.8 0.2];         % 初始概率
chain_length = 32;      % 馬可夫鏈長度

% 初始化馬可夫鏈
chain = zeros(1, chain_length); % 用於存放狀態序列

% 生成初始狀態
prob = cumsum(P0);          % 計算累積概率
choice = rand();            % 生成隨機數
temp = find(prob >= choice);% 找到累積概率大於隨機數的位置
target = temp(1);           % 取出第一個位置
chain(1) = target;          % 保存初始狀態

% 生成後續狀態
for i = 2:chain_length
    P_temp = P(chain(i-1), :); % 用上一個狀態確定轉移概率
    prob = cumsum(P_temp);     % 計算累積概率
    choice = rand();           % 生成隨機數
    temp = find(prob >= choice); % 找到累積概率大於隨機數的位置
    target = temp(1);          % 取出第一個位置
    chain(i) = target;         % 保存此次狀態
end

% 將狀態 1 映射為 1，狀態 2 映射為 -1
chain(chain == 1) = 1;
chain(chain == 2) = -1;

% 輸出結果
disp(chain);

% 將濾波器係數存到一個檔案中
fileID = fopen('chain_data.txt', 'w');
fprintf(fileID, '%.8e\n', chain);
fclose(fileID);

% 確認數據輸出完成
disp('Filter data has been written to chain_data.txt');
