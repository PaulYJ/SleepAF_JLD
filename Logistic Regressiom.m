%% Logistic regression for 2 cohorts
dataset= readtable('dataset_Nhood_updated.csv');

%%
shhs1total_updated_clean.anar = (shhs1total_updated_clean.anar1a1 == 1) | (shhs1total_updated_clean.anar1b1 == 1) | (shhs1total_updated_clean.anar1c1 == 1) | (shhs1total_updated_clean.anar31 == 1);
shhs2total_updated_clean.anar = (shhs2total_updated_clean.anar1a2 == 1) | (shhs2total_updated_clean.anar1b2 == 1) | (shhs2total_updated_clean.anar1c2 == 1) | (shhs2total_updated_clean.anar32 == 1);
vars_sleep = {'slpprdp', 'slpeffp', 'timerem', 'timest1', 'timest2', 'timest34', ...
    'timeremp', 'timest1p', 'timest2p', 'times34p', ...
    'hslptawp', 'waso', 'ai_all', 'ai_nrem', 'ai_rem', 'slplatp', 'remlaiip'};
vars_to_find1 = {'afib', 'mi15', 'angina15', 'hf15', 'stroke15', ...
    'bmi_s1', 'ahi_a0h3', ...
    'htnderv_s1','insuln1', 'ohga1', 'smokstat_s1','height','weight','diasbp','systbp', 'htnmed1','lipid1','benzod1','anar','alcoh'};
vars_to_find2 = {'afib', ...
    'bmi_s2', 'ahi_a0h3', ...
    'htnderv_s2', 'insuln2', 'ohga2', 'smokstat_s2','pm207','pm202','avg23bpd_s2','avg23bps_s2','htnmed2','lipid2','benzod2','anar', 'sh328','sh329','sh330'};
shhs1_selected = shhs1total_updated_clean(:, {'nsrrid','visitnumber', 'age_s1', 'gender', 'race',vars_sleep{:},vars_to_find1{:}});
shhs2_selected = shhs2total_updated_clean(:, {'nsrrid','visitnumber', 'age_s2', 'gender', 'race', vars_sleep{:},vars_to_find2{:}});
shhs1_selected.Properties.VariableNames{'age_s1'} = 'age';
shhs2_selected.Properties.VariableNames{'age_s2'} = 'age';

shhs2_selected.Properties.VariableNames{'bmi_s2'} = 'bmi_s1';
shhs2_selected.Properties.VariableNames{'htnderv_s2'} = 'htnderv_s1';
shhs2_selected.Properties.VariableNames{'smokstat_s2'} = 'smokstat_s1';
shhs2_selected.Properties.VariableNames{'pm207'} = 'height';
shhs2_selected.Properties.VariableNames{'pm202'} = 'weight';
shhs2_selected.Properties.VariableNames{'avg23bpd_s2'} = 'diasbp';
shhs2_selected.Properties.VariableNames{'avg23bps_s2'} = 'systbp';
shhs2_selected.Properties.VariableNames{'htnmed2'} = 'htnmed1';
shhs2_selected.Properties.VariableNames{'lipid2'} = 'lipid1';
shhs2_selected.Properties.VariableNames{'benzod2'} = 'benzod1';


shhs2_selected.alcoh=shhs2_selected.sh328+shhs2_selected.sh329+shhs2_selected.sh330;
shhs2_selected = removevars(shhs2_selected, 'sh328');
shhs2_selected = removevars(shhs2_selected, 'sh329');
shhs2_selected = removevars(shhs2_selected, 'sh330');

shhs2_selected.diabetes = NaN(size(shhs2_selected.insuln2));
shhs2_selected.diabetes((shhs2_selected.insuln2 == 1) | (shhs2_selected.ohga2 == 1)) = 1;
index_zero_and_nan = ((shhs2_selected.insuln2 == 0) & isnan(shhs2_selected.ohga2)) | ...
                     ((shhs2_selected.ohga2 == 0) & isnan(shhs2_selected.insuln2))|...
                     ((shhs2_selected.ohga2 == 0) & (shhs2_selected.insuln2 == 0));
shhs2_selected.diabetes(index_zero_and_nan) = 0;


shhs1_selected.diabetes = NaN(size(shhs1_selected.insuln1));
shhs1_selected.diabetes((shhs1_selected.insuln1 == 1) | (shhs1_selected.ohga1 == 1)) = 1;
index_zero_and_nan = ((shhs1_selected.insuln1 == 0) & isnan(shhs1_selected.ohga1)) | ...
                     ((shhs1_selected.ohga1 == 0) & isnan(shhs1_selected.insuln1))|...
                     ((shhs1_selected.ohga1 == 0) & (shhs1_selected.insuln1 == 0));
shhs1_selected.diabetes(index_zero_and_nan) = 0;

shhs1_selected = removevars(shhs1_selected, {'mi15','angina15','hf15','stroke15'});
shhs2_selected.Properties.VariableNames{'insuln2'} = 'insuln1';
shhs2_selected.Properties.VariableNames{'ohga2'} = 'ohga1';

dataset_updated = vertcat(shhs1_selected, shhs2_selected);
dataset_updated= sortrows(dataset_updated, 'nsrrid');
dataset_updated.NhoodGroup = dataset.NhoodGroup;

%% 选定要包括在逻辑回归模型中的变量
variablesForRegression_1 = {'age', 'gender', 'race', 'NhoodGroup'};
variablesForRegression_2 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','NhoodGroup'};
variablesForRegression_3 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','alcoh','NhoodGroup'};

% 提取所有变量名称
allVariables = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','alcoh',...'diasbp','systbp',
    'NhoodGroup'};

% 定义变量标签
variableLabels = struct(...
    'age', 'Age', ...
    'gender', 'Gender', ...
    'race', 'Race', ...
    'bmi_s1', 'BMI', ...
    'ahi_a0h3', 'AHI', ...
    'htnderv_s1', 'Self-reported hypertention', ...
    'parrptdiab', 'Self-reported diabetes', ...
    'smokstat_s1', 'Smoking status', ...
    'alcoh', 'Alcohol use', ...
    ...'diasbp', 'Seated diastolic blood pressure', ...
    ...'systbp', 'Seated systolic blood pressure', ...
    'chol', 'Total cholesterol', ...
    'hdl', 'HDL', ...
    'trig', 'Triglycerides', ...
    'NhoodGroup', 'NhoodGroup' ...
);


% 创建逻辑回归的新表格
logisticRegressionTable = dataset_updated(:, allVariables);

% 将分类变量转化为categorical类型
logisticRegressionTable.gender = categorical(logisticRegressionTable.gender);
logisticRegressionTable.race = categorical(logisticRegressionTable.race);
logisticRegressionTable.htnderv_s1 = categorical(logisticRegressionTable.htnderv_s1);
logisticRegressionTable.diabetes = categorical(logisticRegressionTable.diabetes);
logisticRegressionTable.smokstat_s1 = categorical(logisticRegressionTable.smokstat_s1);

logisticRegressionTable.NhoodGroup = categorical(logisticRegressionTable.NhoodGroup);
% 将二进制输出变量 'diseaseStatus' 加入到表格中，此处为假设的变量名
% 如果你的疾病状态已经在表格中了，确保它是逻辑型变量或二进制数值型变量
logisticRegressionTable.diseaseStatus = dataset_updated.afib; % 替换需要的疾病状态
logisticRegressionTable.afib = dataset_updated.afib;

disp(groupsummary(logisticRegressionTable,"NhoodGroup"));

%%
continuousVariables = {'age', 'bmi_s1', 'ahi_a0h3', 'alcoh'};

% 处理极端异常值 (+/- 3 SD)
for i = 1:length(continuousVariables)
    variableName = continuousVariables{i};
    data = logisticRegressionTable.(variableName);
    meanValue = mean(data, 'omitnan');
    stdValue = std(data, 'omitnan');
    
    % 设置极端值的上下限
    lowerLimit = meanValue - 3 * stdValue;
    upperLimit = meanValue + 3 * stdValue;
    
    % 将超出上下限的值设置为NaN
    logisticRegressionTable.(variableName)(data < lowerLimit | data > upperLimit) = NaN;
end

% 删除含有NaN值的行
logisticRegressionTableClean = rmmissing(logisticRegressionTable);

%这样处理完 7841->4566

%% 生成哑变量（0-N）
% 假设 logisticRegressionTableClean 是你的主数据表
uniqueGroups = categories(logisticRegressionTableClean.NhoodGroup);

% 对每个组生成一个新的二分类变量
for i = 0:length(uniqueGroups)-1
    currentGroup = uniqueGroups{i+1};  % 从0开始的索引调整为1开始
    newVarName = sprintf('Type%d', i);  % 创建变量名，如 type0, type1, ...
    
    % 生成二分类变量，1 表示当前组，0 表示其他组
    logisticRegressionTableClean.(newVarName) = logisticRegressionTableClean.NhoodGroup == currentGroup;
end

% 额外步骤：将所有未定义的NhoodGroup设置为NaN
undefinedIndices = isundefined(logisticRegressionTableClean.NhoodGroup);
for i = 0:8
    newVarName = sprintf('Type%d', i);
    % logisticRegressionTableClean.(newVarName)(undefinedIndices) = NaN;  % 设置NaN
end


%% 2. 构建多变量回归模型（变量数目逐渐增加）
logisticRegressionTableClean.diseaseStatus = categorical(logisticRegressionTableClean.afib); % 重要!替换需要的疾病状态
logisticRegressionTableClean.NhoodGroup = categorical(logisticRegressionTableClean.NhoodGroup);
% 通过 statset 增加迭代次数
options = statset('glmfit');
options.MaxIter = 100;  % 增加迭代次数到100次

% 设置模型选项
options = statset('Display', 'final');

% 初始化结果存储
subtypeResults = {};

% 置信水平
alpha = 0.05;

% 变量组
variablesForRegression_0 = {};
variablesForRegression_1 = {'age', 'gender', 'race'};
variablesForRegression_2 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes'};
variablesForRegression_3 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','alcoh'};

% 变量组数组
allVariablesForRegression = {variablesForRegression_0, variablesForRegression_1, variablesForRegression_2, variablesForRegression_3};
typeVariables = {'Type1', 'Type2', 'Type3', 'Type4', 'Type5', 'Type6'};

% 模型名称数组
modelNames = {'模型0', '模型1', '模型2', '模型3'};
allResults={};
allResults_plot={};

for k = 1:length(allVariablesForRegression)
    for j = 1:length(typeVariables)
        % 构建模型公式，包括基础变量和一个特定的Type变量
        currentVariables = [allVariablesForRegression{k}, typeVariables(j)];
        formula = sprintf('diseaseStatus ~ %s', strjoin(currentVariables, ' + '));
        mdl = fitglm(logisticRegressionTableClean, formula, 'Distribution', 'binomial', 'Options', options, 'Link', 'logit');
    
       % 提取当前Type变量的回归系数和置信区间
        typeVariableName = sprintf('%s_1',typeVariables{j});
            beta = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Row, typeVariableName));
            CI = coefCI(mdl, alpha);
            CI_Lower = CI(strcmp(mdl.Coefficients.Row, typeVariableName), 1);
            CI_Upper = CI(strcmp(mdl.Coefficients.Row, typeVariableName), 2);
            pValue = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Row, typeVariableName));
            
            % 计算 OR 和 95% CI
            OR = exp(beta);
            OR_CI_Lower = exp(CI_Lower);
            OR_CI_Upper = exp(CI_Upper);
            
            % 存储结果
            modelName = sprintf('模型%d_%s', k-1, typeVariables{j});
            allResults = [allResults; {modelName, sprintf('%.2f (%.2f, %.2f)', OR, OR_CI_Lower, OR_CI_Upper), pValue}];
            allResults_plot = [allResults_plot; {sprintf('模型%d', k-1),sprintf(typeVariables{j}), sprintf('%.2f (%.2f, %.2f)', OR, OR_CI_Lower, OR_CI_Upper), pValue,sprintf('%.2f', OR), sprintf('%.2f', OR_CI_Lower),sprintf('%.2f', OR_CI_Upper)}];
    end
end

% 将结果转换为表并添加列名
allResultsTable = cell2table(allResults, 'VariableNames', {'Model', 'OR (95% CI)', 'pValue'});
allResultsTable_plot = cell2table(allResults_plot, 'VariableNames', {'Model','Type', 'OR (95% CI)', 'pValue','OR','OR_Lower','OR_Upper'});

% 显示所有模型回归结果
disp('所有模型的回归结果:');
disp(allResultsTable);

% 输出结果
disp('单变量回归结果已保存为CSV文件。');

%% 对visit1的病人进行所有cvdoutcome的回归（附表）
dataset_shhs1 = shhs1_selected; %5312个病人
dataset_shhs1.mi15 = shhs1total_updated_clean.mi15;
dataset_shhs1.angina15 = shhs1total_updated_clean.angina15;
dataset_shhs1.hf15 = shhs1total_updated_clean.hf15;
dataset_shhs1.stroke15 = shhs1total_updated_clean.stroke15;
dataset_shhs1.chol = shhs1total_updated_clean.chol;
dataset_shhs1.hdl = shhs1total_updated_clean.hdl;
dataset_shhs1.trig = shhs1total_updated_clean.trig;

dataset_shhs1_nhood = dataset(dataset.visitnumber==1,:); %5312个病人
dataset_shhs1.NhoodGroup = categorical(dataset_shhs1_nhood.NhoodGroup);
%%
variablesForRegression_1 = {'age', 'gender', 'race', 'NhoodGroup'};
variablesForRegression_2 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','NhoodGroup'};
variablesForRegression_3 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','alcoh',...
    ...'diasbp','systbp',
    'chol','hdl','trig','NhoodGroup'};

% 提取所有变量名称
allVariables = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1','diabetes','smokstat_s1','alcoh',...'diasbp','systbp',
    'chol','hdl','trig','NhoodGroup'};

% 创建逻辑回归的新表格
logisticRegressionTable = dataset_shhs1(:, allVariables);

% 将分类变量转化为categorical类型
logisticRegressionTable.gender = categorical(logisticRegressionTable.gender);
logisticRegressionTable.race = categorical(logisticRegressionTable.race);
logisticRegressionTable.htnderv_s1 = categorical(logisticRegressionTable.htnderv_s1);
logisticRegressionTable.diabetes = categorical(logisticRegressionTable.diabetes);
logisticRegressionTable.smokstat_s1 = categorical(logisticRegressionTable.smokstat_s1);

% 将二进制输出变量 'diseaseStatus' 加入到表格中，此处为假设的变量名
% 如果你的疾病状态已经在表格中了，确保它是逻辑型变量或二进制数值型变量
logisticRegressionTable.diseaseStatus = dataset_shhs1.afib; % 重要!！！！！！！！！！！！！！！！!替换需要的疾病状态

logisticRegressionTable.mi15 = dataset_shhs1.mi15; 
logisticRegressionTable.angina15 = dataset_shhs1.angina15; 
logisticRegressionTable.hf15 = dataset_shhs1.hf15; 
logisticRegressionTable.stroke15 = dataset_shhs1.stroke15; 
logisticRegressionTable.afib = dataset_shhs1.afib; 

disp(groupsummary(logisticRegressionTable,"NhoodGroup"));

% 定义连续型变量列表
continuousVariables = {'age', 'bmi_s1', 'ahi_a0h3', 'alcoh', ...'diasbp', 'systbp',
    'chol', 'hdl', 'trig'};

% 处理极端异常值 (+/- 3 SD)
for i = 1:length(continuousVariables)
    variableName = continuousVariables{i};
    data = logisticRegressionTable.(variableName);
    meanValue = mean(data, 'omitnan');
    stdValue = std(data, 'omitnan');
    
    % 设置极端值的上下限
    lowerLimit = meanValue - 3 * stdValue;
    upperLimit = meanValue + 3 * stdValue;
    
    % 将超出上下限的值设置为NaN
    logisticRegressionTable.(variableName)(data < lowerLimit | data > upperLimit) = NaN;
end

% 删除含有NaN值的行
logisticRegressionTableClean = rmmissing(logisticRegressionTable);

uniqueGroups = categories(logisticRegressionTableClean.NhoodGroup);

% 对每个组生成一个新的二分类变量
for i = 0:length(uniqueGroups)-1
    currentGroup = uniqueGroups{i+1};  % 从0开始的索引调整为1开始
    newVarName = sprintf('Type%d', i);  % 创建变量名，如 type0, type1, ...
    
    % 生成二分类变量，1 表示当前组，0 表示其他组
    logisticRegressionTableClean.(newVarName) = logisticRegressionTableClean.NhoodGroup == currentGroup;
end

% 额外步骤：将所有未定义的NhoodGroup设置为NaN
undefinedIndices = isundefined(logisticRegressionTableClean.NhoodGroup);
for i = 0:8
    newVarName = sprintf('Type%d', i);
    % logisticRegressionTableClean.(newVarName)(undefinedIndices) = NaN;  % 设置NaN
end

%% 3. 对所有cvdoutcome进行模型0-3的检验
% 通过 statset 增加迭代次数
options = statset('glmfit');
options.MaxIter = 100; % 增加迭代次数到1000次
options.Display = 'final'; % 设置模型选项为最终显示

% 初始化结果存储
allSubtypeResults = {};

% 置信水平
alpha = 0.05;

% 定义要分析的心血管疾病结果
cvdOutcomes = {'mi15', 'angina15', 'hf15', 'stroke15', 'afib'};

% 变量组
variablesForRegression_0 = {};
variablesForRegression_1 = {'age', 'gender', 'race'};
variablesForRegression_2 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1', 'diabetes'};
variablesForRegression_3 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'htnderv_s1', 'diabetes', 'smokstat_s1', 'alcoh','chol', 'hdl', 'trig'};

% variablesForRegression_2 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'srhype', 'parrptdiab', 'smokstat_s1'};
% variablesForRegression_3 = {'age', 'gender', 'race', 'bmi_s1', 'ahi_a0h3', 'srhype', 'parrptdiab', 'smokstat_s1', 'alcoh', 'diasbp', 'systbp', 'chol', 'hdl', 'trig'};


% 变量组数组
allVariablesForRegression = {variablesForRegression_0, variablesForRegression_1, variablesForRegression_2, variablesForRegression_3};
typeVariables = {'Type1', 'Type2', 'Type3', 'Type4', 'Type5', 'Type6'};

% 模型名称数组
modelNames = {'模型0', '模型1', '模型2', '模型3'};
allSubtypeResults = [];

for outcomeIdx = 1:length(cvdOutcomes)
    outcome = cvdOutcomes{outcomeIdx};
    % 确保目标变量是二进制的（0 和 1）, 删除包含 "Don't Know" (8) 样本
    logisticRegressionTableClean_filtered = logisticRegressionTableClean(logisticRegressionTableClean.(outcome) ~= 8, :);
    logisticRegressionTableClean_filtered.diseaseStatus = categorical(logisticRegressionTableClean_filtered.(outcome));
     for k = 1:length(allVariablesForRegression)
        for j = 1:length(typeVariables)
        % 构建模型公式，包括基础变量和一个特定的Type变量
        currentVariables = [allVariablesForRegression{k}, typeVariables(j)];
        formula = sprintf('diseaseStatus ~ %s', strjoin(currentVariables, ' + '));
        mdl = fitglm(logisticRegressionTableClean_filtered, formula, 'Distribution', 'binomial', 'Options', options, 'Link', 'logit');
    
       % 提取当前Type变量的回归系数和置信区间
            typeVariableName = sprintf('%s_1',typeVariables{j});
            beta = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Row, typeVariableName));
            CI = coefCI(mdl, alpha);
            CI_Lower = CI(strcmp(mdl.Coefficients.Row, typeVariableName), 1);
            CI_Upper = CI(strcmp(mdl.Coefficients.Row, typeVariableName), 2);
            pValue = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Row, typeVariableName));
            
            % 计算 OR 和 95% CI
            OR = exp(beta);
            OR_CI_Lower = exp(CI_Lower);
            OR_CI_Upper = exp(CI_Upper);
            
             % 存储结果
            modelName = sprintf('模型%d_%s', k-1, typeVariables{j});
            variableLabel = sprintf('Type%d', j);
            allSubtypeResults = [allSubtypeResults; {outcome, modelName, variableLabel, sprintf('%.2f (%.2f, %.2f)', OR, OR_CI_Lower, OR_CI_Upper), pValue}];
         end
     end
end
           
% 将结果转换为表并添加列名
allSubtypeResultsTable = cell2table(allSubtypeResults, 'VariableNames', {'Outcome', 'Model', 'Type', 'OR (95% CI)', 'pValue'});

% 显示子类型回归结果
disp('子类型回归结果:');
disp(allSubtypeResultsTable);

%% 将结果写入CSV文件



