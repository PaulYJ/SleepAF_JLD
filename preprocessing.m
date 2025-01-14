% 数据预处理
shhs2total= readtable('shhs2-dataset-0.20.0.csv');
shhs1total= readtable('shhs1-dataset-0.20.0.csv');
shhs_cvdoutcome= readtable('shhs-cvd-summary-dataset-0.20.0.csv');
shhs_eegbiomarker= readtable('shhs1-eeg-biomarkers-dataset-0.20.0.csv');
shhs1_ecgbiomarker= readtable('shhs1-hrv-summary-0.20.0.csv');
shhs2_ecgbiomarker= readtable('shhs2-hrv-summary-0.20.0.csv');

%合并表格
shhs1total_updated = outerjoin(shhs1total,shhs_eegbiomarker(:, [1,7:64]), 'Keys', 'nsrrid', 'MergeKeys', true, 'Type', 'left');
shhs1total_updated = outerjoin(shhs1total_updated,shhs_cvdoutcome(:, [1,3:35,38:39]), 'Keys', 'nsrrid', 'MergeKeys', true, 'Type', 'left');
shhs1total_updated = outerjoin(shhs1total_updated,shhs1_ecgbiomarker, 'Keys', 'nsrrid', 'MergeKeys', true, 'Type', 'left');
shhs2total_updated = outerjoin(shhs2total,shhs2_ecgbiomarker, 'Keys', 'nsrrid', 'MergeKeys', true, 'Type', 'left');

% shhs1total_updated为5804*1364的表格
%% 仅清洗掉缺少睡眠特征，及睡眠脑电图质量不好的个体
%1. 清洗SHHS2
shhs2total_updated.visitnumber = repmat(2, height(shhs2total_updated), 1);

vars_demo= {'nsrrid','visitnumber','age_s2','gender','race',... % 基础人口学特征
    };
shhs2total_updated_clean = rmmissing(shhs2total_updated, 'DataVariables', vars_demo);

vars_sleep = { 'slpprdp','slpeffp','timerem','timest1','timest2','timest34',...
    'timeremp','timest1p','timest2p','times34p',...
    'hslptawp','waso','ai_all','ai_nrem','ai_rem','slplatp','remlaiip',... % PSG睡眠特征，其中ai_all删掉了很多人
    };

shhs2total_updated_clean = rmmissing(shhs2total_updated_clean, 'DataVariables', vars_sleep);

% 剩下2602个样本

missingIndices = (shhs2total_updated_clean.overall_shhs2<4); % 删除了信号质量为 1: Unsatisfactory   2: Poor   3: Fair的个体
shhs2total_updated_clean(missingIndices,:)=[];

% 剩下2529个样本(其中51个房颤病人)
%%
%2. 清洗SHHS1
shhs1total_updated.visitnumber = repmat(1, height(shhs1total_updated), 1);
vars_demo= {'nsrrid','visitnumber','age_s1','gender','race',... % 基础人口学特征
    };
shhs1total_updated_clean = rmmissing(shhs1total_updated, 'DataVariables', vars_demo);

vars_sleep = { 'slpprdp','slpeffp','timerem','timest1','timest2','timest34',...
    'timeremp','timest1p','timest2p','times34p',...
    'hslptawp','waso','ai_all','ai_nrem','ai_rem','slplatp','remlaiip',... % PSG睡眠特征
    };

shhs1total_updated_clean = rmmissing(shhs1total_updated_clean, 'DataVariables', vars_sleep);
% 剩下5615个人

% 2. 删除脑电图水平中伪影程度较高的个体（artifact free）
missingIndices = (shhs1total_updated_clean.eeg1qual==1 | shhs1total_updated_clean.eeg2qual==1);
shhs1total_updated_clean(missingIndices,:)=[];
%剩下5312个人 （其中57个房颤病人）

%%
% 3. 汇总SHHS1和SHHS2的数据，根据nsrrid排列
vars_sleep = {'slpprdp', 'slpeffp', 'timerem', 'timest1', 'timest2', 'timest34', ...
    'timeremp', 'timest1p', 'timest2p', 'times34p', ...
    'hslptawp', 'waso', 'ai_all', 'ai_nrem', 'ai_rem', 'slplatp', 'remlaiip'};
vars_to_find = { 'bmi_s1', 'ahi_a0h3', 'cai0p','ahi_c0h3',...
    'afib', 'mi15', 'angina15', 'hf15', 'stroke15', ...
    'srhype', 'parrptdiab', 'smokstat_s1', 'alcoh', 'diasbp',...
    'systbp', 'chol', 'hdl', 'trig'};
shhs1_selected = shhs1total_updated_clean(:, {'nsrrid','visitnumber', 'age_s1', 'gender', 'race',vars_sleep{:},vars_to_find{:}});
shhs2total_updated_clean = shhs2total_updated_clean(:, {'nsrrid','visitnumber', 'age_s2','gender','race', vars_sleep{:},'bmi_s2', 'ahi_a0h3', 'cai0p','afib'});
shhs1_selected.Properties.VariableNames{'age_s1'} = 'age';
shhs2total_updated_clean.Properties.VariableNames{'age_s2'} = 'age';
shhs2total_updated_clean.Properties.VariableNames{'bmi_s2'} = 'bmi_s1';
% 为 shhs2_selected 补全缺失的变量
for var = vars_to_find
    if ~ismember(var{:}, shhs2total_updated_clean.Properties.VariableNames)
        shhs2total_updated_clean.(var{:}) = NaN(height(shhs2total_updated_clean), 1);
    end
end
dataset = vertcat(shhs1_selected, shhs2total_updated_clean);

%dataset共5312+2529=7841个样本
% 根据 nsrrid 列进行排序
dataset = sortrows(dataset, 'nsrrid');
% 这样就整理好啦\@@/


%% 比较2组样本的基本睡眠特征及人口学特征
stats = table(); % 创建一个新表格存储结果

% 1. 计算连续型正态分布变量的Mean(SD)
vars_to_cal = {'age_s1','bmi_s1','ahi_a0h3','height','weight', 'hip','neck20','alcoh',...
    'chol','hdl','trig','diasbp','systbp',...
    };
total_samples = height(shhs1total_updated_clean); % 获取数据集的总样本数
% 初始化表格，预先定义列名和所需的行数
stats = table(cell(length(vars_to_cal), 1), zeros(length(vars_to_cal), 1), zeros(length(vars_to_cal), 1), ...
              'VariableNames', {'Mean_SD', 'NonNaN_Count', 'Total_Length'}, ...
              'RowNames', vars_to_cal);
% removesample = outerjoin(removesample,shhs1total_updated_clean{vars_to_cal}, 'Keys', 'nsrrid', 'MergeKeys', true, 'Type', 'left');

for i = 1:length(vars_to_cal)
    var = vars_to_cal{i};
    data = shhs1total_updated_clean.(var);

    % 计算均值和标准差
    mean_val = mean(data, 'omitnan');
    sd_val = std(data, 'omitnan');

    % 识别异常值
    outlier_indices = abs(data - mean_val) > 3 * sd_val;
    data(outlier_indices) = NaN;

    % 重新计算平均值和标准差
    mean_val_updated = mean(data, 'omitnan');
    sd_val_updated = std(data, 'omitnan');
    non_nan_count = sum(~isnan(data)); % 计算非NaN的数量

    % 直接填充表格对应行
    stats.Mean_SD{i} = sprintf('%.1f (%.1f)', mean_val_updated, sd_val_updated);
    stats.NonNaN_Count(i) = non_nan_count;
    stats.Total_Length(i) = total_samples;
end

disp(stats);
