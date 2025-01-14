% 创建 pcamatrix 表格，并添加新的列
pcatable = table();
vars_sleep = { 'slpprdp','slpeffp','timerem','timest1','timest2','timest34',...
    'timeremp','timest1p','timest2p','times34p',...
    'hslptawp','waso','ai_all','ai_nrem',...
    'ai_rem','slplatp','remlaiip',... % PSG睡眠特征
    };

vars_to_add = { 'TST','Efficiency','REM duration','N1 duration','N2 duration','N3/4 duration',...
    'REM proportion','N1 proportion','N2 proportion','N3/4 proportion',...
    'Wake per hour','WASO','Arousal Index','NREM Arousal Index',...
    'REM Arousal Index','Latency','REM Latency',... % PSG睡眠特征
    };

vars_label_num = { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'
    };

for i = 1:numel(vars_sleep)
    pcatable.(vars_sleep{i}) = dataset.(vars_sleep{i});
end

% 将清理后的表格转换为 double 类型的数组
pcamatrix= table2array(pcatable);

%% 主成分分析
%pca标准化
pcamatrix_standardized = zscore(pcamatrix);

% 进行主成分分析
[coeff, score, ~, ~, explained] = pca(pcamatrix_standardized);

%% 查找高度相关的变量
% 计算相关系数矩阵
[R, pValue] = corr(pcamatrix);

% 查找高度相关的变量对
[i, j] = find(triu(abs(R), 1) > 0.80);

% 显示高度相关的变量对和它们的相关系数
highlyCorrelatedPairs = [i, j, R(sub2ind(size(R), i, j))];

%% 定义一个color bar
% 定义颜色
colorNeg = [0 0.4470 0.7410];  % 蓝色，对应 -1
colorPos = [0.8500 0.3250 0.0980];   % 橙色，对应 1
colorMid = [1, 1, 1];                   % 白色，对应 0

% 创建颜色映射
numColors = 256;  % 定义颜色级别数量
halfNum = numColors / 2;
colorMap = zeros(numColors, 3);

% 从 -1 到 0 的颜色渐变
colorMap(1:halfNum, :) = interp1([1, halfNum], [colorNeg; colorMid], 1:halfNum);

% 从 0 到 1 的颜色渐变
colorMap(halfNum+1:numColors, :) = interp1([halfNum+1, numColors], [colorMid; colorPos], halfNum+1:numColors);


%% 描述变量之间相关性的热图
% 计算相关系数矩阵
corrMatrix = corr(pcamatrix, 'Rows', 'pairwise');  % 使用'pairwise'来忽略NaN值

% 绘制热图
figure;
h=heatmap(vars_label_num, vars_label_num, corrMatrix, 'Colormap', colorMap, 'CellLabelColor','none');
h.ColorLimits = [-1, 1];  % 设置颜色条范围
h.Colormap = colorMap;   % 显示颜色条
xlabel('Sleep Index');
ylabel('Sleep Index');

% 获取当前活动的坐标轴句柄并设置图形窗口的位置和大小
fig = gcf;
ax = gca;
fig.Position = [100, 100, 320, 280];  % 设置图形窗口的位置和大小
ax.FontSize = 8;  % 设置坐标轴刻度字体大小
ax.FontName = 'Arial';  % 设置字体为Arial

%% 输出cbar
cbarFigure = figure('Position', [100, 100, 400, 400]);
cbar = colorbar;
colormap(cbarFigure, colorMap);
caxis([-1 1]); % 设置颜色条范围
title(cbar, 'Correlation');

%%
% 计算累积解释的方差百分比
cumulativeExplained = cumsum(explained);

% 使用插值使曲线更平滑
x = 1:length(cumulativeExplained);
xi = linspace(1, length(cumulativeExplained), 300);  % 创建更细的x值
smoothedY = interp1(x, cumulativeExplained, xi, 'spline');  % 使用样条插值

% 绘制 Elbow Plot
figure;
plot(cumulativeExplained, '-o', 'LineWidth', 1,'MarkerSize', 4, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410]);
xlabel('Number of Principal Components', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Cumulative Explained Variance (%)', 'FontName', 'Arial', 'FontSize', 10);
grid off;  % 显示网格以增加可读性

% 调整坐标轴字体和字号
ax = gca;  % 获取当前坐标轴对象
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.TickDir = 'out';  % 设置刻度线朝外
ax.XTick = 1:4:length(cumulativeExplained);  % 每三个主成分一个刻度
ax.Box = 'off';% 移除顶部和右侧的坐标轴线
ax.XLim = [1, 17];  % 设置x轴的上限为17

% 调整图形窗口的大小以确保文字和图形位置合适
fig = gcf;
fig.Position = [100, 100, 200, 180];  % 调整图形窗口的位置和大小

%% 根据clustergram的结果重新排列
pcatable = table();
% 定义原始变量名数组
vars_sleep = { 'slpprdp','slpeffp','timerem','timest1','timest2','timest34',...
    'timeremp','timest1p','timest2p','times34p',...
    'hslptawp','waso','ai_all','ai_nrem',...
    'ai_rem','slplatp','remlaiip'};

vars_to_add = { 'TST','Efficiency','REM duration','N1 duration','N2 duration','N3/4 duration',...
    'REM proportion','N1 proportion','N2 proportion','N3/4 proportion',...
    'Wake per hour','WASO','Arousal Index','NREM Arousal Index',...
    'REM Arousal Index','Latency','REM Latency'};

% 定义新的排序索引
new_order = [13, 14, 11, 15, 4, 8, 12, 5, 9, 17, 16, 1, 2, 3, 7, 6, 10];

% 重新排列变量名
sorted_vars_sleep = vars_sleep(new_order);
sorted_vars_to_add = vars_to_add(new_order);


for i = 1:numel(sorted_vars_sleep)
    pcatable.(sorted_vars_sleep{i}) = dataset.(sorted_vars_sleep{i});
end

% 将清理后的表格转换为 double 类型的数组
pcamatrix= table2array(pcatable);

%% 定义一个color bar
% 定义颜色
colorNeg = [0 0.4470 0.7410];  % 蓝色，对应 -1
colorPos = [0.8500 0.3250 0.0980];   % 橙色，对应 1
colorMid = [1, 1, 1];                   % 白色，对应 0

% 创建颜色映射
numColors = 256;  % 定义颜色级别数量
halfNum = numColors / 2;
colorMap = zeros(numColors, 3);

% 从 -1 到 0 的颜色渐变
colorMap(1:halfNum, :) = interp1([1, halfNum], [colorNeg; colorMid], 1:halfNum);

% 从 0 到 1 的颜色渐变
colorMap(halfNum+1:numColors, :) = interp1([halfNum+1, numColors], [colorMid; colorPos], halfNum+1:numColors);

%% 描述变量之间相关性的热图
% 计算相关系数矩阵
corrMatrix = corr(pcamatrix, 'Rows', 'pairwise');  % 使用'pairwise'来忽略NaN值
vars_label_num = { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'};

% 绘制热图
figure;
h=heatmap(vars_label_num, vars_label_num, corrMatrix, 'Colormap', colorMap, 'CellLabelColor','none');
h.ColorLimits = [-1, 1];  % 设置颜色条范围
h.Colormap = colorMap;   % 显示颜色条
xlabel('Sleep Index');
ylabel('Sleep Index');

% 获取当前活动的坐标轴句柄并设置图形窗口的位置和大小
fig = gcf;
ax = gca;
fig.Position = [100, 100, 320, 280];  % 设置图形窗口的位置和大小
ax.FontSize = 8;  % 设置坐标轴刻度字体大小
ax.FontName = 'Arial';  % 设置字体为Arial

axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.TickDirection = 'out';
% cb.Location = 'manual';


% 添加标题以注释是使用Pearson相关系数
title('Pearson’s Correlation Coefficient');

%%
tree = linkage(corrMatrix,'average','correlation');
H = dendrogram(tree,'Orientation','left');
set(H,'LineWidth',1);

%%
% 获取数据矩阵大小
[rows, cols] = size(corrMatrix);

% 创建图形窗口
fig = figure('Position', [500, 200, 800, 750], 'Name', 'Standard Layout');

% 创建热图区域
placeMat=zeros(7,7);placeMat(2:7,2:7)=1;
axMain=subplot(7,7,find(placeMat'));
heatmap(corrMatrix,'Colormap',colorMap);
axMain.XLim=[1,cols]+[-.5,.5];
axMain.YLim=[1,rows]+[-.5,.5];
axMain.YAxisLocation = 'right';
axMain.YDir = 'reverse';
axMain.XTick = 1:cols;
axMain.YTick = 1:rows;
xlabel('Sleep Index');
ylabel('Sleep Index');
hold on;

% 创建树状图区域
axTree1 = subplot(7, 7, (1:6) * 7 + 1);
[H, ~, ~] = dendrogram(tree, 'Orientation', 'left');
axTree1.XAxis.Visible = 'off';  % 隐藏x轴
axTree1.YAxis.Visible = 'off';  % 隐藏y轴
hold on;

% 创建颜色条区域
axBar = subplot(7, 8, 1:8);  % 在顶部全宽度展示
axBar.Position = [0.1, 0.9, 0.8, 0.05];  % 定位颜色条
axBar.Color = 'none';
axBar.XColor = 'none';
axBar.YColor = 'none';
colorbar(axBar, 'northoutside');  % 将颜色条放置在顶部外侧


%%
% 获取数据矩阵大小
[rows, cols] = size(corrMatrix);

% 创建图形窗口
fig = figure('Position', [500, 200, 800, 750], 'Name', 'Standard Layout');

% 创建热图区域
placeMat=zeros(7,7);placeMat(2:7,2:7)=1;
axMain=subplot(7,7,find(placeMat'));
imagesc(axMain, corrMatrix,[-1 1]);
colormap(axMain, colorMap);  % 设置热图的颜色映射
cb = colorbar(axMain, 'northoutside');  % 将颜色条放置在热图上方
cb.TickDirection = 'out';  % 设置颜色条的刻度方向向外

axMain.XLim=[1,cols]+[-.5,.5];
axMain.YLim=[1,rows]+[-.5,.5];
axMain.YAxisLocation = 'right';
axMain.YDir = 'reverse';
axMain.XTick = 1:cols;
axMain.YTick = 1:rows;
set(axMain, 'YTickLabel', sorted_vars_to_add, 'FontName', 'Arial', 'FontSize', 10);
% axMain.YTickLabel = sorted_vars_to_add;
axMain.Box = 'off'; 
axMain.TickDir = 'out';
xlabel('Sleep Index');
ylabel('Sleep Index');
hold on;

% 添加网格线以分隔单元格
% 绘制黑色网格线 —— 水平线
LineX = repmat([[1, cols] + [-0.5, 0.5], nan], [rows + 1, 1]).';
LineY = repmat((0.5:1:(rows + 0.5)).', [1, 3]).';
plot(axMain, LineX(:), LineY(:), 'Color', 'k', 'LineWidth', 1);  % 使用黑色('k')

% 绘制黑色网格线 —— 垂直线
LineY = repmat([[1, rows] + [-0.5, 0.5], nan], [cols + 1, 1]).';
LineX = repmat((0.5:1:(cols + 0.5)).', [1, 3]).';
plot(axMain, LineX(:), LineY(:), 'Color', 'k', 'LineWidth', 1);  % 使用黑色('k')

% 创建树状图区域
axTree1 = subplot(7, 7, (1:6) * 7 + 1);
[H, ~, ~] = dendrogram(tree, 'Orientation', 'left');
set(H, 'Color', 'k', 'LineWidth', 1);  % 设置树状图的线条颜色为黑色
axTree1.XAxis.Visible = 'off';  % 隐藏x轴
axTree1.YAxis.Visible = 'off';  % 隐藏y轴
hold on;

% 确保所有元素的字体一致
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 10);

title('Pearson’s Correlation Coefficient');


%%
% 创建颜色条区域
axBar = subplot(7, 7, 1:7);  % 在顶部全宽度展示
axBar.Position = [0.1, 0.4, 0.8, 0.5];  % 定位颜色条
axBar.Color = 'none';
axBar.XColor = 'none';
axBar.YColor = 'none';
colorbar(axBar, 'northoutside');  % 将颜色条放置在顶部外侧
colormap(axBar, colorMap);  % 确保ColorBar使用同一颜色映射
hold on;
