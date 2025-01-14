%% Chi square

% 假设 dataset 是已经加载的表，包含 nsrrid 和 visitnumber 列
nsrrid_visit1 = dataset.nsrrid(dataset.visitnumber == 1);
nsrrid_visit2 = dataset.nsrrid(dataset.visitnumber == 2);
common_nsrrid = intersect(nsrrid_visit1, nsrrid_visit2);

% 筛选原始数据集中只包含这些 nsrrid 的行
filtered_dataset= dataset(ismember(dataset.nsrrid, common_nsrrid), :);

% 查看结果 
% disp(filtered_dataset);


%% 根据新发和再现两次访问进行队列研究
% 根据访问次数拆分数据
data_visit1 = filtered_dataset(filtered_dataset.visitnumber == 1, :);
data_visit2 = filtered_dataset(filtered_dataset.visitnumber == 2, :);

% 定义正组的名称
target_groups = {1,2,3,4,5,6};

% 判断每个病人的组状态并创建组别标签
group_labels = cell(size(data_visit1, 1), 1);

for i = 1:length(target_groups)
    target_group = target_groups{i};

    is_positive_visit1 = (data_visit1.NhoodGroup == target_group);
    is_positive_visit2 = (data_visit2.NhoodGroup == target_group);

    % 分类病人到四个组别
    group_labels(is_positive_visit1 & is_positive_visit2) = {'PP'};
    group_labels(is_positive_visit1 & ~is_positive_visit2) = {'PN'};
    group_labels(~is_positive_visit1 & is_positive_visit2) = {'NP'};
    group_labels(~is_positive_visit1 & ~is_positive_visit2) = {'NN'};


data_visit2.group_labels = group_labels;
data_visit2.group_labels = categorical(data_visit2.group_labels);
new_exposure = data_visit2((data_visit2.group_labels == 'NP') | (data_visit2.group_labels == 'NN'), :);
sleep_recovery = data_visit2((data_visit2.group_labels == 'PP') | (data_visit2.group_labels == 'PN'), :);

% 新的组别分类 - New Exposure
new_exposure_labels = cell(size(new_exposure, 1), 1);
new_exposure_labels(new_exposure.group_labels=='NP')  = {'Positive'};
new_exposure_labels(new_exposure.group_labels=='NN')  = {'Negative'};
new_exposure.new_exposure_labels = new_exposure_labels;

% 新的组别分类 - Sleep Recovery
sleep_recovery_labels = cell(size(sleep_recovery, 1), 1);
sleep_recovery_labels(sleep_recovery.group_labels=='PP')  = {'Positive'};
sleep_recovery_labels(sleep_recovery.group_labels=='PN')  = {'Negative'};
sleep_recovery.sleep_recovery_labels = sleep_recovery_labels;

% 卡方检验 - New Exposure
[T_new, chi2_new, p_new, labels_new] = crosstab(new_exposure.afib, new_exposure_labels);
fprintf('Type%d: New Exposure - Chi-squared statistic = %.3f, p-value = %.3f\n', i, chi2_new, p_new);

% 卡方检验 - Sleep Recovery
[T_rec, chi2_rec, p_rec, labels_rec] = crosstab(sleep_recovery.afib, sleep_recovery_labels);
fprintf('Type%d: Sleep Recovery - Chi-squared statistic = %.3f, p-value = %.3f\n', i, chi2_rec, p_rec);

% 绘制New Exposure的热图
    hFig1 = figure;
    hMap1 = heatmap(new_exposure, 'new_exposure_labels', 'afib');
    hMap1.Title = sprintf('New Exposure - Type%d', i);
    hMap1.XLabel = sprintf('p = %.3f', p_new);
    hMap1.YLabel = 'Afib Status';
    hMap1.FontName = 'Arial';
    hMap1.FontSize = 8;
    hMap1.ColorbarVisible = 'off';
    set(hFig1, 'Units', 'pixels', 'Position', [100, 100, 150, 140]);
    % 导出New Exposure热图
    filename_pdf1 = sprintf(file, i);
    exportgraphics(hFig1, filename_pdf1, 'ContentType', 'vector');
    filename_png1 = sprintf(file, i);
    exportgraphics(hFig1, filename_png1, 'Resolution', 300);

    % 绘制Sleep Recovery的热图
    hFig2 = figure;
    hMap2 = heatmap(sleep_recovery,'sleep_recovery_labels','afib');
    hMap2.Title = sprintf('Sleep Recovery - Type%d', i);
    hMap2.XLabel = sprintf('p = %.3f', p_rec);
    hMap2.YLabel = 'Afib Status';
    hMap2.FontName = 'Arial';
    hMap2.FontSize = 8;
    hMap2.ColorbarVisible = 'off';
    set(hFig2, 'Units', 'pixels', 'Position', [300, 100, 150, 140]);
    % 导出Sleep Recovery热图
    filename_pdf2 = sprintf(file,, i);
    exportgraphics(hFig2, filename_pdf2, 'ContentType', 'vector');
    filename_png2 = sprintf(file,, i);
    exportgraphics(hFig2, filename_png2, 'Resolution', 300);
end
