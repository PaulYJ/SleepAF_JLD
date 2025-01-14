%% 定义要读取XML文件的目录
directoryPath = '/Users/shhs/polysomnography/annotations-events-nsrr/shhs1';

% 获取目录中所有的XML文件
files = dir(fullfile(directoryPath, '*.xml'));

% 初始化一个结构体数组来存储所有文件的解析结果
scoredEventsStructs_1 = struct('Data', [], 'visitnumber', [], 'nsrrid', []);

% 循环处理每个文件
for i = 1:length(files)
    % 构造完整的文件路径
    filename = fullfile(files(i).folder, files(i).name);
    
    % 从文件名解析 visit 和 nsrrid
    tokens = regexp(files(i).name, 'shhs(\d+)-(\d+)-nsrr\.xml', 'tokens');
    visit = str2double(tokens{1}{1});
    nsrrid = str2double(tokens{1}{2});
    
    % 解析 XML 文件
    try
        scoredEventData = parseXML(filename);
        % 存储解析结果和编号信息
        scoredEventsStructs_1(i).Data = scoredEventData;
        scoredEventsStructs_1(i).visitnumber = visit;
        scoredEventsStructs_1(i).nsrrid = nsrrid;
    catch ME
        warning('Failed to parse %s due to error: %s', filename, ME.message);
        % 如果解析失败，保留空值
        scoredEventsStructs_1(i).Data = [];
        scoredEventsStructs_1(i).visitnumber = visit;
        scoredEventsStructs_1(i).nsrrid = nsrrid;
    end
    
    % 显示进度
    fprintf('Processed %d of %d files (%.2f%% complete).\n', i, length(files), (i / length(files) * 100));
end

% 现在 scoredEventsStructs 包含了所有文件的解析结果及其编号


%% 读取shhs2
% 定义要读取XML文件的目录
directoryPath = '/Users/shhs/polysomnography/annotations-events-nsrr/shhs2';

% 获取目录中所有的XML文件
files = dir(fullfile(directoryPath, '*.xml'));

% 初始化一个结构体数组来存储所有文件的解析结果
scoredEventsStructs_2 = struct('Data', [], 'visitnumber', [], 'nsrrid', []);

% 循环处理每个文
for i = 1:length(files)
    % 构造完整的文件路径
    filename = fullfile(files(i).folder, files(i).name);
    
    % 从文件名解析 visit 和 nsrrid
    tokens = regexp(files(i).name, 'shhs(\d+)-(\d+)-nsrr\.xml', 'tokens');
    visit = str2double(tokens{1}{1});
    nsrrid = str2double(tokens{1}{2});
    
    % 解析 XML 文件
    try
        scoredEventData = parseXML(filename);
        % 存储解析结果和编号信息
        scoredEventsStructs_2(i).Data = scoredEventData;
        scoredEventsStructs_2(i).visitnumber = visit;
        scoredEventsStructs_2(i).nsrrid = nsrrid;
    catch ME
        warning('Failed to parse %s due to error: %s', filename, ME.message);
        % 如果解析失败，保留空值
        scoredEventsStructs_2(i).Data = [];
        scoredEventsStructs_2(i).visitnumber = visit;
        scoredEventsStructs_2(i).nsrrid = nsrrid;
    end
    
    % 显示进度
    fprintf('Processed %d of %d files (%.2f%% complete).\n', i, length(files), (i / length(files) * 100));
end

%%
% 初始化 cell 数组，四列分别为 visitnumber, nsrrid, sleepStages 和 respiratoryEvents
allEvents = cell(numel(scoredEventsStructs_1) + numel(scoredEventsStructs_2), 4);
columnNames = {'visitnumber', 'nsrrid', 'sleepStages', 'respiratoryEvents'}; % 列名
index = 1;  % 初始化索引用于存储位置

% 处理每个结构体数组
for structArray = {scoredEventsStructs_1, scoredEventsStructs_2}
    scoredEventsStructs = structArray{1};
    for i = 1:numel(scoredEventsStructs)
        currentStruct = scoredEventsStructs(i);
        
        % 过滤出 Respiratory|Respiratory 和 Stages|Stages 的数据
        respiratoryEvents = currentStruct.Data(strcmp({currentStruct.Data.EventType}, 'Respiratory|Respiratory'));
        stageEvents = currentStruct.Data(strcmp({currentStruct.Data.EventType}, 'Stages|Stages'));

        % 初始化 Stage, Start, Duration
        Stage = [];  % 空的 Stage 列表，初始化为空
        Start = [];
        Duration = [];

        % 提取 Stage, Start, Duration (Stages|Stages)
        for j = 1:numel(stageEvents)
            concept = stageEvents(j).EventConcept;
            numStr = regexp(concept, '\d+', 'match');
            Stage = [Stage; str2double(numStr{1})];  % 转换为数字
            Start = [Start; stageEvents(j).Start];
            Duration = [Duration; stageEvents(j).Duration];
        end
        
        % 创建表格保存所有的 Stage, Start, Duration 数据
        sleepStages = table(Stage, Start, Duration);

       % 提取 Hypopnea 和 Apnea 的 Start 和 Duration
        respStart = [];
        respDuration = [];
        respEvent = {}; % 初始化空 cell 数组，用于存储事件类型
        
        for j = 1:numel(respiratoryEvents)
            concept = respiratoryEvents(j).EventConcept;
            if contains(concept, 'Hypopnea') || contains(concept, 'Apnea')  % 选取 Hypopnea 或 Apnea 事件
                respStart = [respStart; respiratoryEvents(j).Start];
                respDuration = [respDuration; respiratoryEvents(j).Duration];
                respEvent = [respEvent; respiratoryEvents(j).EventConcept];  % 将符合条件的事件添加到 respEvent 中
            end
        end
        
        respiratoryEventsTable = table(respStart, respDuration, respEvent,'VariableNames', {'Start', 'Duration', 'Event'}); % 仅保存呼吸事件的 Start 和 Duration

        % 存储到 cell 数组中
        allEvents{index, 1} = currentStruct.visitnumber;
        allEvents{index, 2} = currentStruct.nsrrid;
        allEvents{index, 3} = sleepStages;  % 直接存储 sleepStages 表格
        allEvents{index, 4} = respiratoryEventsTable;  % 直接存储呼吸事件表格
        index = index + 1;  % 更新索引
    end
end

% 将 allEvents 转换为 table 并添加列名
allEventsTable = cell2table(allEvents, 'VariableNames', columnNames);


%% 将dataset中的NhoodGroup信息添加到allEventsTable中
% 确保 dataset 是表格并包含必要的列
if istable(dataset) && all(ismember({'visitnumber', 'nsrrid', 'NhoodGroup'}, dataset.Properties.VariableNames))
    
    % 连接表格，使用 'left' 方法保留 allEventsTable 中所有行
    combinedTable = outerjoin(allEventsTable, dataset, ...
                              'Keys', {'visitnumber', 'nsrrid'}, ...
                              'MergeKeys', true, ...
                              'Type', 'left', ...
                              'LeftVariables', setdiff(allEventsTable.Properties.VariableNames, 'NhoodGroup'), ...
                              'RightVariables', 'NhoodGroup');

    % 处理缺失值
    % 对于 outerjoin，默认会在结果中填充 NaN 或者空字符串，如果 NhoodGroup 是数值，NaN 就可以；如果是类别或字符串，你可能需要手动替换
    combinedTable.NhoodGroup(ismissing(combinedTable.NhoodGroup)) = NaN;
else
    error('Dataset must be a table and contain columns: visitnumber, nsrrid, and NhoodGroup.');
end

% 显示合并后的表格，检查一部分
% disp(combinedTable(1:10, :));
% allEventsTable = combinedTable.NhoodGroup;


%% 脚本函数部分
function scoredEvents = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure for Scored Events.
try
   xmlDoc = xmlread(filename);
catch
   error('Failed to read XML file %s.', filename);
end

% 获取根节点下的所有 ScoredEvent 子节点
scoredEventNodes = xmlDoc.getElementsByTagName('ScoredEvent');
numEvents = scoredEventNodes.getLength;
scoredEvents = struct('EventType', {}, 'EventConcept', {}, 'Start', {}, 'Duration', {}, 'SignalLocation', {});

% 遍历每个 ScoredEvent 节点
for i = 1:numEvents
    currentNode = scoredEventNodes.item(i-1);
    
    % 解析EventType
    eventType = parseNodeText(currentNode, 'EventType');
    
    % 解析EventConcept
    eventConcept = parseNodeText(currentNode, 'EventConcept');
    
    % 解析Start
    start = parseNodeText(currentNode, 'Start');
    
    % 解析Duration
    duration = parseNodeText(currentNode, 'Duration');
    
    % 解析SignalLocation
    signalLocation = parseNodeText(currentNode, 'SignalLocation');
    
    % 创建结构体
    scoredEvents(i).EventType = eventType;
    scoredEvents(i).EventConcept = eventConcept;
    scoredEvents(i).Start = str2double(start);
    scoredEvents(i).Duration = str2double(duration);
    scoredEvents(i).SignalLocation = signalLocation;
end
end

function text = parseNodeText(parentNode, tagName)
    % 此函数用于安全地解析文本内容，避免空节点引起的错误
    node = parentNode.getElementsByTagName(tagName);
    if node.getLength > 0 && node.item(0).hasChildNodes
        textNode = node.item(0).getFirstChild;
        if ~isempty(textNode) && textNode.getNodeType == textNode.TEXT_NODE
            text = char(textNode.getData);
        else
            text = '';  % 如果没有文本节点，返回空字符串
        end
    else
        text = '';  % 如果没有找到节点，返回空字符串
    end
end


