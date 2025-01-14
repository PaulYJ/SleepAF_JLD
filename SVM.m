%% 仅看第一次横断面的发病数据
dataset.reducedX=reduced_data(:,1);
dataset.reducedY=reduced_data(:,2);

% 分开2次访视的数据
dataset_shhs1 = dataset(dataset.visitnumber==1,:); %5312个病人
dataset_shhs2 = dataset(dataset.visitnumber==2,:); %2529个病人


%% 随机抽样
% 初始化用来存储所有训练集的单元数组
trainingSets = cell(10, 1);

% 从数据集中筛选出 afib==1 和 afib==0 的行
afib1 = dataset(dataset.afib == 1, :);
afib0 = dataset(dataset.afib == 0, :);

% 检查是否有足够的数据进行抽样
if height(afib1) < 100 || height(afib0) < 100
    error('不足以从每个组中抽取 50 个个体');
end

% 进行 10 次抽样
for i = 1:10
     % 每次迭代设置一个新的随机种子
    rng(i);  % 设置种子为当前迭代次数，确保每次迭代可重复且不同
    % 随机抽样 afib==1 的个体
    rows1 = datasample(afib1, 100, 'Replace', false);
    % 随机抽样 afib==0 的个体
    rows0 = datasample(afib0, 100, 'Replace', false);
    
    % 合并两个子集来形成训练集
    trainingSets{i} = [rows1; rows0];
end

% 可以选择输出第一次抽样的结果来查看
% disp(head(trainingSets{1}));

%% 将10个分类器的结果综合预测shhs1的cvdoutcome风险
% 假设 trainingSets 已经包含了之前抽样得到的训练集
% 初始化一个用于存储 SVM 模型的单元数组
svmModels = cell(10, 1);
accuracy = zeros(10, 1); % 初始化准确率数组

for i = 1:10
    % 选择第 6 到第 23 列的数据（第23列是cvdoutcome）
    subset = trainingSets{i}(:, 6:23);
    
    % 提取分类标签
    labels = subset.afib;
    
    % 提取特征（假设除 'afib' 外的所有列都是特征）
    features = subset(:, setdiff(subset.Properties.VariableNames, 'afib'));
    
    % 使用 70% 的数据进行训练，30% 的数据进行测试
    cv = cvpartition(size(features, 1), 'Holdout', 0.3);
    svmModel = fitcsvm(features(cv.training,:), labels(cv.training), 'Standardize', true );

    
    % 存储模型
    svmModels{i} = svmModel;

    % 使用模型对测试集进行预测
    testPredictions = predict(svmModel, features(cv.test,:));
    
    % 计算测试集的准确率
    testLabels = labels(cv.test);
    accuracy(i) = sum(testPredictions == testLabels) / numel(testLabels);
end

% 输出每个模型的准确率
disp('Model accuracies:');
disp(accuracy);

%% 抽样新的测试集并用多数投票法计算综合正确率

%%% 从原始数据集中抽样200人并用多数投票法评估模型综合性能
% 再次使用相同方法抽样200人
rng(15);
rows1_test = datasample(afib1, 100, 'Replace', false);
rows0_test = datasample(afib0, 100, 'Replace', false);
testSet = [rows1_test; rows0_test];

% 测试集特征和标签
testFeatures = testSet(:, 6:22);
testLabels = testSet.afib;

% 使用10个模型进行预测
testPredictions = zeros(height(testSet), 10);
for i = 1:10
    testPredictions(:, i) = predict(svmModels{i}, testFeatures);
end

% 多数投票法
sumPredictions = sum(testPredictions, 2);
finalPredict = sumPredictions > 5;
combinedAccuracy = sum(finalPredict == testLabels) / height(testSet);

disp('Combined accuracy after majority voting:');
disp(combinedAccuracy);


%% 抽100个validation set

validationSets = cell(100, 1);

% 进行 100 次抽样
for i = 1:100
     % 每次迭代设置一个新的随机种子
    rng(i);  % 设置种子为当前迭代次数，确保每次迭代可重复且不同
    % 随机抽样 afib==1 的个体
    rows1 = datasample(afib1, 100, 'Replace', false);
    % 随机抽样 afib==0 的个体
    rows0 = datasample(afib0, 100, 'Replace', false);
    % 合并两个子集来形成训练集
    validationSets{i} = [rows1; rows0];
end

validation_results = zeros(100, 10);

for i = 1:100
% 测试集特征和标签
    testFeatures = validationSets{i}(:, 6:22);
    testLabels = validationSets{i}.afib;

    % 使用10个模型进行预测
    testPredictions = zeros(height(validationSets{i}), 10);
    for m = 1:10
        testPredictions(:, m) = predict(svmModels{m}, testFeatures);
    end

% 多数投票法
    for p = 1:10
    sumPredictions = sum(testPredictions, 2);
    finalPredict = sumPredictions > p-1;
    combinedAccuracy = sum(finalPredict == testLabels) / height(testSet);
    validation_results(i,p)= combinedAccuracy;
    end
end

disp('Combined accuracy after majority voting:');
disp(combinedAccuracy);


%% Shuffle Control 对100个validation set

% 初始化存储打乱结果的变量
shuffleResults = zeros(100, 10);

% 进行 100 次抽样
for i = 1:100
    % 使用相同的数据集
    testSet = validationSets{i};  % 使用未打乱的验证集数据

    % 打乱标签
    shuffledLabels = testSet.afib(randperm(height(testSet)));

    % 测试集特征和打乱后的标签
    testFeatures = testSet(:, 6:22);
    testLabels = shuffledLabels;  % 使用打乱后的标签

    % 使用10个模型进行预测
    testPredictions = zeros(height(testSet), 10);
    for m = 1:10
        testPredictions(:, m) = predict(svmModels{m}, testFeatures);
    end

    % 多数投票法
    for p = 1:10
        sumPredictions = sum(testPredictions, 2);
        finalPredict = sumPredictions > p-1;
        combinedAccuracy = sum(finalPredict == testLabels) / height(testSet);
        shuffleResults(i, p) = combinedAccuracy;
    end
end

% 输出或保存shuffle control结果
disp(shuffleResults);

%% shuffle control validation

shuffle_validation_results = zeros(100, 10);

for i = 1:100
% 测试集特征和标签
    testFeatures = validationSets{i}(:, 6:22);
    testLabels = validationSets{i}.afib;

    % 使用10个模型进行预测
    testPredictions = zeros(height(validationSets{i}), 10);
    for m = 1:10
        testPredictions(:, m) = predict(shuffledSvmModels{m}, testFeatures);
    end

% 多数投票法
    for p = 1:10
    sumPredictions = sum(testPredictions, 2);
    finalPredict = sumPredictions > p-1;
    combinedAccuracy = sum(finalPredict == testLabels) / height(testSet);
    shuffle_validation_results(i,p)= combinedAccuracy;
    end
end

% disp('Shuffled model: Combined accuracy after majority voting');
% disp(combinedAccuracy);

%% 绘制两个图片的折线图
% 计算平均准确率和标准差
meanValidation = mean(validation_results);
stdValidation = std(validation_results);

meanShuffle = mean(shuffle_validation_results);
stdShuffle = std(shuffle_validation_results);

% 创建新图形
figure;

% 为每个阈值绘制标准差区域
hold on;
x = 1:10;

% 绘制标准差区域（灰色背景）
fill([x fliplr(x)], [meanShuffle+stdShuffle fliplr(meanShuffle-stdShuffle)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% 绘制标准差区域（蓝色背景）
fill([x fliplr(x)], [meanValidation+stdValidation fliplr(meanValidation-stdValidation)], [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% 绘制平均准确率的折线图
plot(x, meanValidation, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', [0 0.4470 0.7410]);
plot(x, meanShuffle, '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5]);


% 图形设置
% title('Accuracy for Majority Voting Threshold');
xlabel('Threshold (i)');
ylabel('Accuracy');
legend('Shuffle SD', 'Validation SD', 'Validation Mean', 'Shuffle Mean', 'Box', 'off','Location', 'best');
grid off;
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 12;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XLim = [1 10]; % 设置x轴的范围
hold off;

% 调整图形窗口的大小以确保文字和图形位置合适
fig = gcf;
fig.Position = [100, 100, 240, 150];  % 调整图形窗口的位置和大小


%% 初始化存储所有混淆矩阵的变量
totalConfMat = zeros(2, 2);  % 假设afib是一个二分类问题

for i = 1:100
    % 测试集特征和标签
    testFeatures = validationSets{i}(:, 6:22);
    testLabels = validationSets{i}.afib;

    % 使用10个模型进行预测并应用多数投票法
    testPredictions = zeros(height(validationSets{i}), 10);
    for m = 1:10
        testPredictions(:, m) = predict(svmModels{m}, testFeatures);
    end

    % 多数投票法
    sumPredictions = sum(testPredictions, 2);
    finalPredict = sumPredictions > 5;  % 多数投票法：超过一半即判定为正类
    finalPredict = double(finalPredict);

    % 计算当前验证集的混淆矩阵
    confMat = confusionmat(testLabels, finalPredict);
    totalConfMat = totalConfMat + confMat;  % 累加所有混淆矩阵
end

% 计算平均混淆矩阵
averageConfMat = totalConfMat / 100;

% 显示平均混淆矩阵
disp('Average Confusion Matrix:');
disp(averageConfMat);

%% 定义一个color bar
% 定义颜色
colorPos = [0 0.4470 0.7410];  % 蓝色，对应 -1
colorNeg = [0.8500 0.3250 0.0980];   % 橙色，对应 1
colorMid = [1, 1, 1];                   % 白色，对应 0

% 创建颜色映射
numColors = 256;  % 定义颜色级别数量
halfNum = numColors / 2;
colorMap = zeros(numColors, 3);

% 从 -1 到 0 的颜色渐变
colorMap(1:halfNum, :) = interp1([1, halfNum], [colorNeg; colorMid], 1:halfNum);

% 从 0 到 1 的颜色渐变
colorMap(halfNum+1:numColors, :) = interp1([halfNum+1, numColors], [colorMid; colorPos], halfNum+1:numColors);

%% 计算100个测试集的平均混淆矩阵并作图
% finalPredict = double(finalPredict);
% [confMat, order] = confusionmat(testLabels, finalPredict);

% 显示混淆矩阵
disp('Confusion Matrix:');
disp(array2table(averageConfMat, 'VariableNames', {'Predicted_Negative', 'Predicted_Positive'}, 'RowNames', {'Actual_Negative', 'Actual_Positive'}));

% 计算其他性能指标
TP = averageConfMat(2,2);
TN = averageConfMat(1,1);
FP = averageConfMat(1,2);
FN = averageConfMat(2,1);

accuracy = (TP + TN) / (TP + TN + FP + FN);
precision = TP / (TP + FP);
recall = TP / (TP + FN);
F1 = 2 * (precision * recall) / (precision + recall);

% 显示性能指标
disp('Performance Metrics:');
fprintf('Accuracy: %.4f\n', accuracy);
fprintf('Precision: %.4f\n', precision);
fprintf('Recall: %.4f\n', recall);
fprintf('F1 Score: %.4f\n', F1);

fig=figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 180, 160]);
heatmap({'non-AF', 'AF'}, {'non-AF', 'AF'}, averageConfMat, ...
        'ColorbarVisible', 'off', 'Colormap', colorMap);
title('Average Confusion Matrix');
xlabel('Predicted Class');
ylabel('Actual Class');
ax = gca;
ax.FontSize = 8;  % 设置坐标轴刻度字体大小
ax.FontName = 'Arial';  % 设置字体为Arial

%% 用生成好的模型对shhs所有个体进行svm预测
dataset_svm = dataset;

% 提取数据集中的特征
features = dataset_svm(:, 6:22); % 根据你的数据集实际情况调整列的选择

% 初始化预测结果矩阵
predictions = zeros(height(dataset), numel(svmModels));

% 使用每个模型对整个数据集进行预测
for i = 1:numel(svmModels)
    predictions(:, i) = predict(svmModels{i}, features);
end

% 多数投票法确定最终预测
sumPredictions = sum(predictions, 2);
finalPredict = sumPredictions > 5; % 以一半以上模型同意为标准

% 添加最终预测到数据集中
dataset_svm.afibPredict = finalPredict;

% 显示更新后的数据集的一部分
% disp(head(dataset_svm));

%% 写出SVM结果
writetable(dataset_svm,'afib_predict.csv');

%% 绘制最后结果的umap图

% 确定各组数据
group1_idx = dataset_svm.afibPredict == 1;  % afib_predict 为 1 的患者
group2_idx = dataset_svm.afib == 1  % afib 为 1 但 afib_predict 不为 1 的患者
group3_idx = ~dataset_svm.afibPredict == 1;  % 其余患者

% 创建图形
figure;
hold on;

% 绘制其余患者 (灰色实心点)
scatter(dataset_svm.reducedX(group3_idx), dataset_svm.reducedY(group3_idx), 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'SizeData',8,'MarkerFaceAlpha', 0.4);
% 绘制 afib_predict 为 1 的患者 (红色实心点)
scatter(dataset_svm.reducedX(group1_idx), dataset_svm.reducedY(group1_idx), 'filled', 'MarkerFaceColor', [0.8 0.4 0.4],  'SizeData',8,'MarkerFaceAlpha', 0.4);
% 绘制 afib 为 1 的患者 (红色菱形)
scatter(dataset_svm.reducedX(group2_idx), dataset_svm.reducedY(group2_idx), 'Marker','diamond', 'MarkerEdgeColor', 'none','SizeData',20,'MarkerFaceColor', [0.6 0.3 0.3],'MarkerFaceAlpha', 1);

% 设置图形属性
xlabel('UMAP 1', 'FontName', 'Arial', 'FontSize', 12);
ylabel('UMAP 2', 'FontName', 'Arial', 'FontSize', 12);
legend({'non AF Predict','AF Predict', 'AF Patients'},'Box', 'off', 'Location', 'best'); % 添加图例

axis([-2.5 5 -7.5 2.5]);

%导出图像
ax = gca;
ax.FontName = 'Arial';  % 设置字体
ax.FontSize = 14;  % 设置字体大小
ax.Box = 'off';             % 关闭边框
ax.TickDir = 'out';  % 设置刻度线朝外

fig=gcf;
fig.Position = [100, 100, 420,300];  % 调整图形窗口的位置和大小


%% slide window analysis
x = dataset(:,'reducedX').reducedX;
y = dataset(:,'reducedY').reducedY;

var_names={'reducedX','reducedY'};

x_fos = dataset{dataset.afib == 1, var_names{1}};
y_fos = dataset{dataset.afib == 1, var_names{2}};

% Parameters  
n = 0.5; % Window size, best, 1.5
m = 0.5; %sride, best, 1
x_edges = min(x):m:max(x)-n; % Bin edges for x
y_edges = min(y):m:max(y)-n; % Bin edges for y

% Initialize count matrix
count_af_prob= zeros(length(x_edges)-1, length(y_edges)-1);
count_af= zeros(length(x_edges)-1, length(y_edges)-1);

% Count points in each window
for i = 1:length(x_edges)-1
    for j = 1:length(y_edges)-1
        % Get points in the current window
        fos_prob_in_window = dataset.afibPredict((x >= x_edges(i)) & (x < x_edges(i+1)) & (y >= y_edges(j)) & (y < y_edges(j+1)));
        fos_in_window= (x_fos >= x_edges(i)) & (x_fos < x_edges(i+1)) & (y_fos >= y_edges(j)) & (y_fos < y_edges(j+1));
        count_af_prob(i, j) = sum(fos_prob_in_window);
        count_af(i,j) = sum(fos_in_window);
    end
end

% boxplot(count_fos_prob(:),count_fos(:));

af_prob_by_af_count = {};
for i  = 1:max(count_af(:))
    af_prob_by_af_count{i} = count_af_prob(count_af==i);
end

% plot scatter by fos
figure; hold on
for i = 1:max(count_af(:))
    % Scatter the data points for each vector
    jittered_x = i + 0.05*randn(size(af_prob_by_af_count{i})); % Jitter x-axis to prevent overlap
    scatter(jittered_x, af_prob_by_af_count{i}, 40,[0.4660 0.6740 0.1880], 'filled', 'MarkerFaceAlpha', 0.4);
end

means = cellfun(@mean, af_prob_by_af_count);
% Fit a linear regression line to the means
x = (1:length(af_prob_by_af_count)); % x values (vector indices)

p = polyfit(x(~isnan(means)), means(~isnan(means)), 1); % Fit a line (1st degree polynomial)
fitted_line = polyval(p, x); % Evaluate the fitted line

[R, P] = corrcoef(x(~isnan(means)), means(~isnan(means))); % Calculate correlation coefficient and p-value

% Plot the data
plot(1:length(af_prob_by_af_count), means, 'o', 'LineWidth', 1, 'Color', [0.2 0.2 0.2]); % Means
plot(x, fitted_line, '--', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); % Fitted line

% Axis labels and legend
xlabel('AF count for each sliding window');
ylabel('AF prob');
% legend('Data Points (Mean)', 'Fitted Line', 'Box', 'off', 'Location', 'best');

% Grid and axis settings
grid off;
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 12;
ax.TickDir = 'out';
ax.Box = 'off';
% ax.XLim = [1 length(af_prob_by_af_count)]; % Set x-axis limits dynamically
% ax.YLim = [min(means(~isnan(means))) max(means(~isnan(means)))*1.1]; % Adjust y-axis limits
ax.XLim = [1 length(af_prob_by_af_count)+0.5];
ax.YLim = [min(means(~isnan(means)))-15 max(means(~isnan(means)))*1.1];

% Adjust figure size
fig = gcf;
fig.Position = [100, 100, 240, 150]; % Adjust figure position and size




