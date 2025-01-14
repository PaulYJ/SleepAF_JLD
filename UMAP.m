%% 
addpath('/Users/MATLAB/Path/UMAP path/umapFileExchange (4.4)/umap');

% 导入数据
pcamatrix_reduced = score(:, 1:9);
test = pcamatrix_reduced;

%%
% 使用 UMAP 进行降维（调整参数）
[reduced_data, umap] = run_umap(test, ...
    'n_components', 2, ...
    'n_neighbors',25,...
    'min_dist', 0.05, ...
    'metric', 'mahalanobis');

% 输出数据
writematrix(reduced_data,'umap_embeddings.csv');

%% 根据单一sleep metrics在UMAP图上创造颜色映射
% 修改变量 TST=slpprdp;WASO = waso;Latency = slplatp;
vars_sleep = { 'age','ahi_a0h3','slpprdp','slpeffp','timerem','timest1','timest2','timest34',...
    'timeremp','timest1p','timest2p','times34p',...
    'hslptawp','waso','ai_all','ai_nrem',...
    'ai_rem','slplatp','remlaiip',... % PSG睡眠特征
    };

vars_to_add = {'Age', 'AHI','TST','Efficiency','REM Duration','N1 Duration','N2 Duration','N3 Duration',...
    'REM Proportion','N1 Proportion','N2 Proportion','N3 Proportion',...
    'Wake per hour','WASO','Arousal Index','NREM Arousal Index',...
    'REM Arousal Index','Latency','REM Latency',... % PSG睡眠特征
    };

for i = 1:numel(vars_sleep)
    variable_projection = dataset.(vars_sleep{i});
    variableLabel = vars_to_add{i};

% 计算TST值的范围，用于颜色映射
min_tst = min(variable_projection);
max_tst = max(variable_projection);

% 使用 scatter3 函数创建三维散点图，其中颜色代表TST值
figure;
scatter(reduced_data(:,1), reduced_data(:,2), 2, variable_projection, 'filled');
colormap(jet); % 使用彩色图谱表示不同的TST
colorbar; % 显示颜色条
caxis([min_tst max_tst]); % 标准化颜色映射到TST值的范围

% 添加轴标签和标题
xlabel('UMAP 1', 'FontName', 'Arial', 'FontSize', 12);
ylabel('UMAP 2', 'FontName', 'Arial', 'FontSize', 12);
title(variableLabel,'FontName', 'Arial', 'FontSize', 12);  

% 调整坐标轴字体和字号
ax = gca;  % 获取当前坐标轴对象
ax.FontName = 'Arial';
ax.FontSize = 12;
ax.TickDir = 'out';  % 设置刻度线朝外
ax.Box = 'off';             % 关闭边框
axis([-2.5 5 -7.5 2.5]);

% 调整图形窗口的大小以确保文字和图形位置合适
fig = gcf;
fig.Position = [100, 100, 280, 200];  % 调整图形窗口的位置和大小


