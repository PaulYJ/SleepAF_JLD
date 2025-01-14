%% 导入数据
reduced_data= table2array(readtable('umap_embeddings.csv'));

%%
dataset.NhoodGroup = categorical(dataset.NhoodGroup);
uniqueGroups = categories(dataset.NhoodGroup);

% 对每个组生成一个新的二分类变量
for i = 0:length(uniqueGroups)-1
    currentGroup = uniqueGroups{i+1};  % 从0开始的索引调整为1开始
    newVarName = sprintf('Type%d', i);  % 创建变量名，如 type0, type1, ...
    
    % 生成二分类变量，1 表示当前组，0 表示其他组
    dataset.(newVarName) = (dataset.NhoodGroup == currentGroup);
end

%%  根据两次访视数据绘制纵向结果
% 添加坐标和病人信息 
dataset.reducedX=reduced_data(:,1);
dataset.reducedY=reduced_data(:,2);

% 分开2次访视的数据
dataset_shhs1 = dataset(dataset.visitnumber==1,:); %5312个病人
dataset_shhs2 = dataset(dataset.visitnumber==2,:); %2529个病人


%% 在同一张图上绘制visit1和visit2
figure; hold on; % Open a figure and hold it for multiple plots

scatter(dataset_shhs1.reducedX, dataset_shhs1.reducedY,5, [193/256 183/256 207/256],'filled','DisplayName', 'Visit 1'); % 紫色
scatter(dataset_shhs2.reducedX, dataset_shhs2.reducedY,5, [165/256 203/256 165/256],'filled','DisplayName', 'Visit 2'); %灰绿色


xlabel('UMAP 1', 'FontName', 'Arial', 'FontSize', 12);
ylabel('UMAP 2', 'FontName', 'Arial', 'FontSize', 12);
legend( 'Box', 'off', 'Location', 'best'); % 添加图例

axis([-2.5 5 -7.5 2.5]);

%导出图像
ax = gca;
ax.FontName = 'Arial';  % 设置字体
ax.FontSize = 12;  % 设置字体大小
ax.Box = 'off';             % 关闭边框
ax.TickDir = 'out';  % 设置刻度线朝外
fig=gcf;
fig.Position = [100, 100, 430,300];  % 调整图形窗口的位置和大小

%% 在同一张图上绘制（2次访视都有afib数据的人）
% Example MATLAB Code
figure; hold on; % Open a figure and hold it for multiple plots

% Find common nsrrid
[commonIds, idx1, idx2] = intersect(dataset_shhs1.nsrrid, dataset_shhs2.nsrrid);

% Plot first dataset in blue
 scatter(dataset_shhs1.reducedX, dataset_shhs1.reducedY,5, [0.1 0.6 0.8],'filled','DisplayName', 'Visit 1');

% Plot second dataset in red
 scatter(dataset_shhs2.reducedX, dataset_shhs2.reducedY,5, [0.9 0.7 0.7],'filled','DisplayName', 'Visit 2');

% Draw lines between common nsrrid
for i = 1:length(commonIds)
     if dataset_shhs1.afib(idx1(i)) > 0
         % 特殊标记afib>0的访视
        if dataset_shhs2.afib(idx2(i)) == 0
            scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 30, 'b', 'Marker', 'hexagram','MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)),20, 'r','MarkerEdgeColor', 'r', 'Marker', 'o', 'MarkerFaceColor', 'w','DisplayName', 'Visit 2');
            plot([dataset_shhs1.reducedX(idx1(i)), dataset_shhs2.reducedX(idx2(i))], ...
            [dataset_shhs1.reducedY(idx1(i)), dataset_shhs2.reducedY(idx2(i))], 'k','Color',[0.5 0.5 0.5]);
        end
        if dataset_shhs2.afib(idx2(i)) > 0
            scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 30, 'b', 'Marker', 'hexagram','MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)), 30, 'r', 'Marker', 'hexagram','MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            plot([dataset_shhs1.reducedX(idx1(i)), dataset_shhs2.reducedX(idx2(i))], ...
            [dataset_shhs1.reducedY(idx1(i)), dataset_shhs2.reducedY(idx2(i))], 'k','Color',[0.5 0.5 0.5]);
        end
    elseif dataset_shhs1.afib(idx1(i)) == 0 && dataset_shhs2.afib(idx2(i)) > 0
            scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)), 30, 'r', 'Marker', 'hexagram','MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)),20, 'b','MarkerEdgeColor', 'b', 'Marker', 'o', 'MarkerFaceColor', 'w','DisplayName', 'Visit 1');
            plot([dataset_shhs1.reducedX(idx1(i)), dataset_shhs2.reducedX(idx2(i))], ...
            [dataset_shhs1.reducedY(idx1(i)), dataset_shhs2.reducedY(idx2(i))], 'k','Color',[0.5 0.5 0.5]);
     end
end

% 注意这里要考虑NA值缺失的情况

xlabel('UMAP 1');
ylabel('UMAP 2');
legend off;
grid on;
hold off;

%%  从非type2到type2
figure; hold on;  % 打开新的图形窗口并保持以便添加多个图层

% 查找两个数据集中共有的 nsrrid
[commonIds, idx1, idx2] = intersect(dataset_shhs1.nsrrid, dataset_shhs2.nsrrid);

% 绘制 Non-Type2 的点
nonType2Idx = dataset.Type2 == 0;
h1 = scatter(dataset.reducedX(nonType2Idx), dataset.reducedY(nonType2Idx), 5, 'filled', ...
    'MarkerEdgeColor', [0.827, 0.827, 0.827], 'MarkerFaceColor', [0.827, 0.827, 0.827], 'DisplayName', 'Non-Type2' ,...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);

% 绘制 Type2 的点
type2Idx = dataset.Type2 == 1;
h2 = scatter(dataset.reducedX(type2Idx), dataset.reducedY(type2Idx), 5, 'filled', ...
    'MarkerEdgeColor',  [0.902, 0.706, 0.133], 'MarkerFaceColor',  [0.902, 0.706, 0.133], 'DisplayName', 'Type2', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);

% 连接线和特殊点
h3 = []; h4 = []; h5 = []; h6 = []; h7 = []; h8 = []; % 初始化空数组用于手动设置图例
for i = 1:length(commonIds)
    if dataset_shhs2.Type2(idx2(i)) == 1 && dataset_shhs2.afib(idx2(i)) > 0
        % Type2 with AF in second dataset
        h3(end+1) = scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)), 30, 'yellow', 'Marker', 'hexagram',...
            'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]);

        % 连接线
        h4(end+1) = plot([dataset_shhs1.reducedX(idx1(i)), dataset_shhs2.reducedX(idx2(i))], ...
            [dataset_shhs1.reducedY(idx1(i)), dataset_shhs2.reducedY(idx2(i))], 'k', 'Color', [0.5 0.5 0.5]);

        % Non-Type2 with or without AF in first dataset
        if dataset_shhs1.Type2(idx1(i)) == 0
            if dataset_shhs1.afib(idx1(i)) > 0
                h5(end+1) = scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 30, 'black', 'Marker', 'hexagram',...
                    'MarkerEdgeColor', [0.2, 0.2, 0.2], 'MarkerFaceColor', [0.2, 0.2, 0.2]);
            end
            if dataset_shhs1.afib(idx1(i)) == 0
                h6(end+1) = scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 20, 'black',...
                    'MarkerEdgeColor', [0.2, 0.2, 0.2], 'Marker', 'o', 'MarkerFaceColor', 'w');
            end
        end

        % Type2 with or without AF in first dataset
        if dataset_shhs1.Type2(idx1(i)) == 1
            if dataset_shhs1.afib(idx1(i)) > 0
                h7(end+1) = scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 30, 'black', 'Marker', 'hexagram',...
                    'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
            end
            if dataset_shhs1.afib(idx1(i)) == 0
                h8(end+1) = scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 20, 'black',...
                    'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'Marker', 'o', 'MarkerFaceColor', 'w');
            end
        end
    end
end

xlabel('UMAP 1', 'FontName', 'Arial', 'FontSize', 8);
ylabel('UMAP 2', 'FontName', 'Arial', 'FontSize', 8);
title('New exposure for Type2','FontName', 'Arial', 'FontSize', 8);
axis([-2.5 5 -7 2.5]);

% 手动设置图例，只显示每种点类型一次
legend([h1 h2 h4(1) h3(1)  h5(1) h8(1) h6(1) ], {'Non-Type2', 'Type2', 'Connection',...
    'Type2 patients',  'Non-Type2 patients', 'Type2 healthy ones', 'Non-Type2 healthy ones'},...
    'Box', 'off','Location','best');

grid off;
hold off;

% 导出图像
ax = gca;
ax.FontName = 'Arial';  % 设置字体
ax.FontSize = 8;  % 设置字体大小
ax.Box = 'off';  % 关闭边框
ax.TickDir = 'out';  % 设置刻度线朝外
fig = gcf;
fig.Position = [100, 100, 420, 300];  % 调整图形窗口的位置和大小

%% 从type2到非type2

figure; hold on;  % 打开新的图形窗口并保持以便添加多个图层

% 查找两个数据集中共有的 nsrrid
[commonIds, idx1, idx2] = intersect(dataset_shhs1.nsrrid, dataset_shhs2.nsrrid);

% Plot first dataset in blue
nonType2Idx = dataset.Type2 == 0;
scatter(dataset.reducedX(nonType2Idx), dataset.reducedY(nonType2Idx), 5, 'filled', ...
    'MarkerEdgeColor', [0.827, 0.827, 0.827], 'MarkerFaceColor', [0.827, 0.827, 0.827], 'DisplayName', 'Non-Type2' ,...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);
hold on;
% 绘制Type2的点
type2Idx = dataset.Type2 == 1;
scatter(dataset.reducedX(type2Idx), dataset.reducedY(type2Idx), 5, 'filled', ...
    'MarkerEdgeColor',  [0.902, 0.706, 0.133], 'MarkerFaceColor',  [0.902, 0.706, 0.133], 'DisplayName', 'Type2', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);

% Draw lines between common nsrrid
for i = 1:length(commonIds)
     if dataset_shhs1.Type2(idx1(i))==1 && dataset_shhs1.afib(idx1(i)) > 0
         scatter(dataset_shhs1.reducedX(idx1(i)), dataset_shhs1.reducedY(idx1(i)), 30, 'yellow', 'Marker', 'hexagram','MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
         plot([dataset_shhs1.reducedX(idx1(i)), dataset_shhs2.reducedX(idx2(i))], ...
            [dataset_shhs1.reducedY(idx1(i)), dataset_shhs2.reducedY(idx2(i))], 'k','Color',[0.5 0.5 0.5]);
         if dataset_shhs2.Type2(idx2(i))==0
             if dataset_shhs2.afib(idx2(i)) > 0
                scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)), 30, 'black', 'Marker', 'hexagram','MarkerEdgeColor', [0.2, 0.2, 0.2], 'MarkerFaceColor', [0.2, 0.2, 0.2]);
             end
             if dataset_shhs2.afib(idx2(i)) == 0
                scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)),20, 'black','MarkerEdgeColor', [0.2, 0.2, 0.2], 'Marker', 'o', 'MarkerFaceColor', 'w');
             end
         end
         if dataset_shhs2.Type2(idx2(i))==1
             if dataset_shhs2.afib(idx2(i)) > 0
                scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)), 30, 'black', 'Marker', 'hexagram','MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
             end
             if dataset_shhs2.afib(idx2(i)) == 0
                scatter(dataset_shhs2.reducedX(idx2(i)), dataset_shhs2.reducedY(idx2(i)),20, 'black','MarkerEdgeColor', [0.8500 0.3250 0.0980], 'Marker', 'o', 'MarkerFaceColor', 'w');
             end
         end
     end
end

xlabel('UMAP 1', 'FontName', 'Arial', 'FontSize', 8);
ylabel('UMAP 2', 'FontName', 'Arial', 'FontSize', 8);
title('Post-exposure recovery for Type2' ,'FontName', 'Arial', 'FontSize', 8);
axis([-2.5 5 -7 2.5]);
legend off;
grid off;
hold off;

% 导出图像
ax = gca;
ax.FontName = 'Arial';  % 设置字体
ax.FontSize = 8;  % 设置字体大小
ax.Box = 'off';  % 关闭边框
ax.TickDir = 'out';  % 设置刻度线朝外
fig = gcf;
fig.Position = [100, 100, 420, 300];  % 调整图形窗口的位置和大小
