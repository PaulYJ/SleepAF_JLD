
%% Sleep Stage analysis
%% 使用imagesc对齐睡眠数据（未调色版）
% 筛选 NhoodGroup == i 的个体
filteredData = combinedTable(combinedTable.NhoodGroup == 2, :);
% numIndividuals = height(filteredData);
numIndividuals = 100;

maxDuration = 32400; % 设置最大睡眠时间为9小时（32400秒）

% 初始化存储睡眠数据的cell数组
sleepDataCells = cell(numIndividuals, 1);

for i = 1:numIndividuals
    sleepStages = filteredData.sleepStages{i};
    startTime = sleepStages.Start(1); % 假设第一个记录即为入睡时间
    relativeStart = sleepStages.Start - startTime; % 转换为入睡后的相对时间
    relativeEnd = relativeStart + sleepStages.Duration;
    
    % 将 Stage 1, 2, 3, 4 统一修改为 2
    modifiedStages = sleepStages.Stage;
    modifiedStages(ismember(sleepStages.Stage, [1, 2, 3, 4])) = 2;

    % 存储转换后的数据，确保最大时长为9小时
    sleepDataCells{i} = table(relativeStart, relativeEnd, modifiedStages, ...
                              'VariableNames', {'RelativeStart', 'RelativeEnd', 'Stage'});
end

timeStep = 30; % 时间步长为30秒
timeline = 0:timeStep:maxDuration;  % 时间线调整为固定的9小时
numTimePoints = length(timeline);
alignedSleepStages = NaN(numTimePoints, numIndividuals);

for i = 1:numIndividuals
    for j = 1:height(sleepDataCells{i})
        stage = sleepDataCells{i}.Stage(j);
        startIdx = find(timeline >= sleepDataCells{i}.RelativeStart(j), 1, 'first');
        endIdx = find(timeline <= sleepDataCells{i}.RelativeEnd(j), 1, 'last');
        if ~isempty(startIdx) && ~isempty(endIdx)
            alignedSleepStages(startIdx:endIdx, i) = stage;
        end
    end
end


% 将alignedSleepStages中的NaN替换为0以匹配colormap
alignedSleepStages(isnan(alignedSleepStages)) = 0;

figure;
imagesc(timeline / 3600, 1:numIndividuals, alignedSleepStages'); % 转换时间为小时
colormap(jet);
caxis([0 5]); % 设置颜色映射的范围为0到5
% colorbar;
xlabel('Time (hours)');
ylabel('Individuals');
title('Aligned Sleep Stages Across Individuals');

% 调整图窗大小为400x400像素
set(gcf, 'Position', [100, 100, 450, 400]);


%% 根据各个Type的睡眠结构作堆叠柱状图
% 不算sleep onset的时间
% 假设 sleepStages 包含各个阶段的时间数据和 Type 信息
% 初始化每个个体的睡眠阶段时间数组
numIndividuals = height(combinedTable);
individualNREM = zeros(numIndividuals, 1);
individualREM = zeros(numIndividuals, 1);
individualWake = zeros(numIndividuals, 1);

% 遍历每个个体，计算每个睡眠阶段的时间
for i = 1:numIndividuals
    sleepStages = combinedTable.sleepStages{i};

    % 确定第一次入睡的索引（非Wake阶段）
    firstSleepIndex = find(sleepStages.Stage ~= 0, 1, 'first');

    % 计算每个睡眠阶段的时间
    for j = firstSleepIndex:length(sleepStages.Stage)
        stage = sleepStages.Stage(j);
        duration = sleepStages.Duration(j);

        if (stage == 2 || stage == 3 || stage == 4)
            % 累加NREM时间（假设NREM为阶段2, 3, 4）
            individualNREM(i) = individualNREM(i) + duration;
        elseif stage == 5
            % 累加REM时间
            individualREM(i) = individualREM(i) + duration;
        elseif stage == 0
            % 累加Wake时间
            individualWake(i) = individualWake(i) + duration;
        end
    end
end

%%
% 过滤掉 NhoodGroup 中的 NA 值
validGroups = combinedTable.NhoodGroup(~isnan(combinedTable.NhoodGroup));

% 获取每个 NhoodGroup 的唯一值并忽略 NA
types = unique(validGroups);
numTypes = length(types);

nremTimes = zeros(numTypes, 1);
remTimes = zeros(numTypes, 1);
wakeTimes = zeros(numTypes, 1);


for k = 1:numTypes
    typeIndex = combinedTable.NhoodGroup == types(k);
    nremTimes(k) = sum(individualNREM(typeIndex));
    remTimes(k) = sum(individualREM(typeIndex));
    wakeTimes(k) = sum(individualWake(typeIndex));
end

% 计算每种类型的总睡眠时间
totalTimes = nremTimes + remTimes + wakeTimes;

% 转换为百分比
nremPercents = nremTimes ./ totalTimes * 100;
remPercents = remTimes ./ totalTimes * 100;
wakePercents = wakeTimes ./ totalTimes * 100;

%% 绘制堆叠柱形图
figure;
% 注意数据顺序调整为[Wake, NREM, REM]
b = bar(categorical(types), [wakePercents, nremPercents, remPercents], 'stacked');

% 为不同的条形设置颜色，顺序也要对应调整
b(1).FaceColor = [6/255, 8/255, 131/255];  % 浅蓝色
b(2).FaceColor = [27/255, 251/255, 227/255];  % 深蓝色（这里使用类似绿松石的颜色）
b(3).FaceColor = [123/255, 4/255, 4/255];   % 深红色

% 每根柱子不要边框黑线
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(3).EdgeColor = 'none';

xlabel('Type');
ylabel('Percentage of Total Sleep Time');
% title('Proportional Sleep Structure by Type');

% 设置图例并移除图例边框
lgd = legend({'Wake', 'NREM', 'REM'}, 'Location', 'best');
lgd.Box = 'off';

% 设置坐标轴外观
set(gca, 'Box', 'off', 'TickDir', 'out');  % Tick out，Box off

% 调整图窗大小为450x400像素
set(gcf, 'Position', [100, 100, 450, 400]);


%% 单段NREM的中位数
% 初始化存储中位数的数组
numIndividuals = height(combinedTable);
individualMedians = zeros(numIndividuals, 1);
nrem_frequency = zeros(numIndividuals, 1);
nrem_total = zeros(numIndividuals, 1);

   % 遍历每个个体
for i = 1:numIndividuals
    sleepStages = combinedTable.sleepStages{i};
    
    % 将 Stage 1, 2, 3, 4 统一视为 2
    nremStages = ismember(sleepStages.Stage, [1, 2, 3, 4]);
    sleepStages.Stage(nremStages) = 2;

    % 标记 NREM 开始和结束
    nremStarts = [sleepStages.Stage(1) == 2; diff(sleepStages.Stage == 2) == 1];
    nremEnds = [diff(sleepStages.Stage == 2) == -1; sleepStages.Stage(end) == 2];

    % 找到每个 NREM 段的开始和结束索引
    startsIdx = find(nremStarts);
    endsIdx = find(nremEnds);

    % 初始化数组以存储每个 NREM 段的时长
    nremLengths = zeros(length(startsIdx), 1);

    % 计算每个 NREM 段的时长
    for j = 1:length(startsIdx)
        if j <= length(endsIdx)
            nremLengths(j) = sum(sleepStages.Duration(startsIdx(j):endsIdx(j)));
        end
    end

    % 计算并存储中位数，如果没有NREM段，则赋值NaN
    if isempty(nremLengths)
        individualMedians(i) = NaN;
    else
        individualMedians(i) = median(nremLengths);
    end

    nrem_frequency(i) = length(nremLengths);
    nrem_total(i) = sum(nremLengths);
end

% 获取每个 NhoodGroup 的组ID
validGroups = combinedTable.NhoodGroup(~isnan(combinedTable.NhoodGroup));
groups = unique(validGroups);


% 按 NhoodGroup 分组计算每组中位数的平均值
log_individualMedians = log10(individualMedians);
% groupMeans = arrayfun(@(g) mean(individualMedians(combinedTable.NhoodGroup == g)), groups);
groupMeans = arrayfun(@(g) mean(log_individualMedians(combinedTable.NhoodGroup == g)), groups); % 这里显示单段NREM的对数

% 显示结果
disp('Group Averages of Medians for REM:');
disp(table(groups, groupMeans, 'VariableNames', {'NhoodGroup', 'AverageOfMedians'}));

%% 储存单段NREM的中位数并绘图
finalTable = table(combinedTable.visitnumber, combinedTable.nsrrid, ...
                   combinedTable.NhoodGroup, log_individualMedians,...
                    'VariableNames', {'visitnumber', 'nsrrid', 'NhoodGroup', 'IndividualMedians'});


% 使用 innerjoin 合并表格
mergedDataset = outerjoin(dataset, finalTable, ...
                          'Keys', {'nsrrid', 'visitnumber'}, ...
                          'MergeKeys', true, ...
                          'Type', 'left', ...
                          'LeftVariables', dataset.Properties.VariableNames, ...
                          'RightVariables', {'IndividualMedians'});

%% 绘制NREM时长中位数的ROC曲线图，Group 2 vs 其他组，确保AUC > 0.5
groups = unique(mergedDataset.NhoodGroup);
groups(groups == 0 | isnan(groups)) = [];  % 移除组0和NaN

% 定义颜色
colors = {'#E64C4C', '#E6B422', '#4DB24C', '#1A99CC', '#333399', '#B24C98'};

% 获取组2的数据
groupData2 = mergedDataset.IndividualMedians(mergedDataset.NhoodGroup == 2);

% 创建图形窗口
figure;

% 定义子图的布局
nCols = 5;  % 设定每行显示3个图
nRows = ceil((length(groups) - 1) / nCols);  % 计算需要的行数

index = 1;  % 初始化索引，用于跳过Group2的位置
% 为每个组绘制ROC曲线，使用组2作为参考
for i = 1:length(groups)
    if groups(i) == 2
        continue;  % 跳过组2
    end
    subplot(nRows, nCols, index);  % 指定子图位置
    index = index + 1;
    hold on;
    
    % 获取当前组的数据
    currentData = mergedDataset.IndividualMedians(mergedDataset.NhoodGroup == groups(i));
    
    % 组合数据，并调整labels使得AUC > 0.5
    scores = [currentData; groupData2];
    labels = [true(size(currentData)); false(size(groupData2))];  % 调整labels使得AUC > 0.5
    
     % 计算ROC曲线，检查AUC值
    [~, ~, ~, AUC] = perfcurve(labels, scores, true);
    
    % 如果AUC < 0.5，颠倒标签
    if AUC < 0.5
        labels = [false(size(currentData)); true(size(groupData2))];
    end
    
    % 重新计算ROC曲线
    [X, Y, ~, AUC] = perfcurve(labels, scores, true);

    % 绘制ROC曲线
    plot(X, Y, 'LineWidth', 2, 'Color', colors{groups(i)});
    
    % 添加对角线
    plot([0 1], [0 1], '--', 'Color', colors{2}, 'LineWidth', 1.5);  % 调整对角线样式
    
    % 设置图例和图形标签
    legend(['AUC = ' num2str(AUC, '%.3f')], 'Location', 'Best', 'Box', 'off');
    xlabel('1 - Specificity');
    ylabel('Sensitivity');
    % title(['Type ' num2str(groups(i)) ' vs. Type 2']);
    title(['Type ' num2str(groups(i))]);
    grid off;
    set(gca, 'TickDir', 'out');  % 设置刻度线向外
    hold off;
end

% 调整图窗大小和子图间距
set(gcf, 'Position', [100, 100, 700, 120]); % 调整整个图窗大小


%% 绘制单段NREM的提琴图
addpath('/Users/Documents/MATLAB/Path/violinplot');
%% 从 mergedDataset 中提取数据
allData = mergedDataset.IndividualMedians;  % 假设这是需要绘制的数据列
groupLabels = mergedDataset.NhoodGroup;  % 组标签列

% 过滤掉0和NaN数据
validIndices = groupLabels ~= 0 & ~isnan(allData) & ~isnan(groupLabels);
allData = allData(validIndices);
groupLabels = groupLabels(validIndices);

% 定义颜色（将十六进制颜色转换为0到1范围的RGB值）
hexColors = {'#E64C4C', '#E6B422', '#4DB24C', '#1A99CC', '#333399', '#B24C98'};
rgbColors = cellfun(@(x) [hex2dec(x(2:3))/255, hex2dec(x(4:5))/255, hex2dec(x(6:7))/255], hexColors, 'UniformOutput', false);
rgbColors = vertcat(rgbColors{:});  % 转换成矩阵形式

% 绘制小提琴图
figure;
violinplot(allData, groupLabels, 'ViolinColor', rgbColors, ...
    'ViolinAlpha', 1,...
    'ShowData', false, ...         % 不显示数据点
    'QuartileStyle', 'shadow', ... % 设置四分位数样式为 shadow
    'ShowBox', true, ...           % 显示箱形图
    'ShowMean', true, ...          % 显示均值
    'ShowMedian', false, ...       % 不显示中位数
    'ShowWhiskers', false);        % 不显示胡须

title('NREM fragmentation');
xlabel('Types');
ylabel('Individual Medians (s)');
grid off;

% 设置tick方向为外，关闭box
set(gca, 'TickDir', 'out', 'Box', 'off');

% 修改图窗大小
set(gcf, 'Position', [100, 100, 300, 200]); % 调整整个图窗大小


%% 统计Type1和Type2中waso的差别
% 提取NhoodGroup == 1 和 NhoodGroup == 2的数据
group1Data = mergedDataset(mergedDataset.NhoodGroup == 1, :);
group2Data = mergedDataset(mergedDataset.NhoodGroup == 2, :);

% 提取WASO数据
wasoGroup1 = group1Data.waso;
wasoGroup2 = group2Data.waso;

% 创建分组向量，用于标识数据的分组
wasoGroupLabels = [repmat({'Type 1'}, length(wasoGroup1), 1); repmat({'Type 2'}, length(wasoGroup2), 1)];

% 定义颜色（将十六进制颜色转换为0到1范围的RGB值）
hexColors = {'#E64C4C', '#E6B422'};  % 自定义的两个颜色
rgbColors = cellfun(@(x) [hex2dec(x(2:3))/255, hex2dec(x(4:5))/255, hex2dec(x(6:7))/255], hexColors, 'UniformOutput', false);
rgbColors = vertcat(rgbColors{:});  % 转换成矩阵形式

% 绘制WASO的小提琴图
figure;
violinplot([wasoGroup1; wasoGroup2], wasoGroupLabels, 'ViolinColor', rgbColors, ...
    'ViolinAlpha', 1,...
    'ShowData', false, ...         % 不显示数据点
    'QuartileStyle', 'shadow', ... % 设置四分位数样式为 shadow
    'ShowBox', true, ...           % 显示箱形图
    'ShowMean', true, ...          % 显示均值
    'ShowMedian', false, ...       % 不显示中位数
    'ShowWhiskers', false);        % 不显示胡须

% title('WASO Distribution by Type');
xlabel('Type');
ylabel('Wake After Sleep Onset (min)');
grid off;

% 设置tick方向为外，关闭box
set(gca, 'TickDir', 'out', 'Box', 'off');

% 修改图窗大小
set(gcf, 'Position', [100, 100, 150, 200]); % 调整整个图窗大小

%% 对waso进行统计检验
[~, p_waso12] = ttest2(wasoGroup1, wasoGroup2);
fprintf('WASO: p-value = %.5f (t-test)\n',p_waso12);

%% 统计Type1和Type2中Age和ahi_a0h3的差别
% 提取 Type 1 和 Type 2 数据
group1Data = mergedDataset.waso(mergedDataset.NhoodGroup == 1);
group2Data = mergedDataset.waso(mergedDataset.NhoodGroup == 2);

% 组合数据
scores = [group1Data; group2Data];

% 设置标签，Type 1 为 0，Type 2 为 1（颠倒位置）
labels = [zeros(size(group1Data)); ones(size(group2Data))];

% 计算ROC曲线
[X, Y, ~, AUC] = perfcurve(labels, scores, 1);

% 绘制ROC曲线
figure;
hold on;
plot(X, Y, 'LineWidth', 2, 'Color', '#E64C4C');  % 使用第一个颜色

% 添加对角线
plot([0 1], [0 1], '--', 'Color', '#E6B422', 'LineWidth', 1.5);

% 设置图例和图形标签
legend(sprintf('AUC = %.3f', AUC), 'Location', 'Best', 'Box', 'off');
xlabel('1 - Specificity');
ylabel('Sensitivity');
% title('ROC Curve: Type 2 vs Type 1');
grid off;
set(gca, 'TickDir', 'out');  % 设置刻度线向外
hold off;

% 调整图窗大小和子图间距
set(gcf, 'Position', [100, 100, 120, 120]); % 调整整个图窗大小


%%
% 绘制散点图
figure;
scatter(group1Data.ahi_a0h3, group1Data.waso, 'r', 'filled');
hold on;
scatter(group2Data.ahi_a0h3, group2Data.waso, 'b', 'filled');
xlabel('AHI');
ylabel('Wake Count');
legend({'Type 1', 'Type 2'});
title('AHI vs Wake Count by Type');
grid on;

%%
% 去除 Type 1 中包含 NaN 的行
cleanedGroup1Data = group1Data(~isnan(group1Data.ahi_a0h3) & ~isnan(group1Data.waso), :);

% 去除 Type 2 中包含 NaN 的行
cleanedGroup2Data = group2Data(~isnan(group2Data.ahi_a0h3) & ~isnan(group2Data.waso), :);

% 计算 Type 1 的 Pearson 相关系数
[r_type1, p_value_type1] = corr(cleanedGroup1Data.ahi_a0h3, cleanedGroup1Data.waso, 'Type', 'Pearson');
fprintf('Type 1: Pearson correlation coefficient = %.4f, p-value = %.4f\n', r_type1, p_value_type1);

% 计算 Type 2 的 Pearson 相关系数
[r_type2, p_value_type2] = corr(cleanedGroup2Data.ahi_a0h3, cleanedGroup2Data.waso, 'Type', 'Pearson');
fprintf('Type 2: Pearson correlation coefficient = %.4f, p-value = %.4f\n', r_type2, p_value_type2);

%%
% 绘制Type 1的散点图和回归线
figure;
subplot(2, 1, 1); % 2行的子图，第一行
scatter(cleanedGroup1Data.ahi_a0h3, cleanedGroup1Data.waso, 10, 'filled', 'MarkerFaceColor', '#E64C4C'); % 缩小散点大小
hold on;
lsline; % 添加回归线
hline1 = lsline;
set(hline1, 'LineWidth', 2); % 增加回归线粗细
xlabel('AHI');
ylabel('WASO (min)');
title('Type 1');
grid off;

% 在图上显示相关系数和p值
text(min(cleanedGroup1Data.ahi_a0h3), max(cleanedGroup1Data.waso), ...
    sprintf('r = %.2f, p = %.3f', r_type1, p_value_type1), 'FontSize', 12, 'Color', 'k'); % 改为黑色

% 绘制Type 2的散点图和回归线
subplot(2, 1, 2); % 第二行
scatter(cleanedGroup2Data.ahi_a0h3, cleanedGroup2Data.waso, 10, 'filled', 'MarkerFaceColor', '#E6B422'); % 缩小散点大小
hold on;
hline2 = lsline;
set(hline2, 'LineWidth', 2); % 增加回归线粗细
xlabel('AHI');
ylabel('WASO (min)');
title('Type 2');
grid off;

% 在图上显示相关系数和p值
text(min(cleanedGroup2Data.ahi_a0h3), max(cleanedGroup2Data.waso), ...
    sprintf('r = %.2f, p = %.3f', r_type2, p_value_type2), 'FontSize', 12, 'Color', 'k'); % 改为黑色

% 调整图窗大小和子图间距
set(gcf, 'Position', [100, 100, 200, 300]); % 调整整个图窗大小

%% read respiratory
function [WAKE, NREM, REM, HYPO]=parseTable(input)

totaltime = 36000;
sleep_stage = table2array(input.sleepStages{1});
breath_event = table2array(input.respiratoryEvents{1});
WAKE = zeros(totaltime,1);
NREM = zeros(totaltime,1);
REM = zeros(totaltime,1);
HYPO = zeros(totaltime,1);
for i = 1:size(sleep_stage,1)
    statetemp = sleep_stage(i,1);
    if statetemp==0
        WAKE(sleep_stage(i,2)+1:sleep_stage(i,2)+sleep_stage(i,3)+1)=1;
    elseif statetemp == 1 || statetemp == 2 || statetemp == 3
        NREM(sleep_stage(i,2)+1:sleep_stage(i,2)+sleep_stage(i,3)+1)=1;
    elseif statetemp == 5
        REM(sleep_stage(i,2)+1:sleep_stage(i,2)+sleep_stage(i,3)+1)=1;
    end
end
for i = 1:size(breath_event,1)
    HYPO(ceil(breath_event(i,1)):ceil(breath_event(i,1)+round(breath_event(i,2))))=1;
end

end