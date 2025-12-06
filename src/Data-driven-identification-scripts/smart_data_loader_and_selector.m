function [selected_data, data_info] = smart_data_loader_and_selector()
% 智能数据加载与选择

fprintf('   扫描数据文件夹...\n');

% 查找最新的数据文件夹
data_folders = dir('MCG_Data_Stable_*');
if isempty(data_folders)
    error('未找到数据文件夹');
end

% 选择最新的文件夹
[~, latest_idx] = max([data_folders.datenum]);
data_folder = data_folders(latest_idx).name;
fprintf('     使用数据文件夹: %s\n', data_folder);

% 扫描数据文件
data_files = dir(fullfile(data_folder, '*_data.mat'));
fprintf('     发现 %d 个数据文件\n', length(data_files));

%% ====== 分类数据文件 ======
recommended_files = {};
stable_files = {};
all_files = {};

for i = 1:length(data_files)
    filename = data_files(i).name;
    all_files{end+1} = fullfile(data_folder, filename);
    
    if contains(filename, '_RECOMMENDED_')
        recommended_files{end+1} = fullfile(data_folder, filename);
    elseif contains(filename, '_STABLE_') || contains(filename, '_OK_')
        stable_files{end+1} = fullfile(data_folder, filename);
    end
end

fprintf('    推荐数据: %d 个\n', length(recommended_files));
fprintf('     稳定数据: %d 个\n', length(stable_files));

%% ====== 智能选择策略 ======
fprintf('   执行智能选择策略...\n');

% 策略1：优先使用推荐数据
if length(recommended_files) >= 8
    fprintf('    策略1: 使用推荐数据 (足够的高质量数据)\n');
    selected_files = recommended_files;
    selection_strategy = 'recommended_priority';
    
elseif length(recommended_files) >= 4
    fprintf('    策略2: 混合推荐+稳定数据\n');
    selected_files = [recommended_files, stable_files(1:min(8, length(stable_files)))];
    selection_strategy = 'mixed_quality';
    
else
    fprintf('    策略3: 使用所有可用稳定数据\n');
    selected_files = [recommended_files, stable_files];
    selection_strategy = 'all_stable';
end

fprintf('     最终选择: %d 个数据文件\n', length(selected_files));

%% ====== 按激励类型分组选择 ======
fprintf('   按激励类型优化选择...\n');
[optimized_files, type_distribution] = optimize_selection_by_excitation_type(selected_files);

%% ====== 加载选定的数据 ======
fprintf('  加载选定的数据文件...\n');
selected_data = {};
load_summary = struct();

for i = 1:length(optimized_files)
    try
        fprintf('    加载: %s\n', optimized_files{i});
        loaded_data = load(optimized_files{i});
        
        if isfield(loaded_data, 'save_data')
            selected_data{end+1} = loaded_data.save_data;
            
            % 提取关键信息用于显示
            quick_stats = loaded_data.save_data.quick_stats;
            fprintf('       %s I%d F%d (质量: %.2f, 稳定性: %s)\n', ...
                    quick_stats.excitation_name, quick_stats.intensity, ...
                    quick_stats.frequency, quick_stats.quality_score, ...
                    quick_stats.stability_level);
        end
        
    catch ME
        fprintf('       加载失败: %s\n', ME.message);
    end
end

%% ====== 编译数据信息 ======
data_info = struct();
data_info.total_files = length(selected_data);
data_info.data_folder = data_folder;
data_info.selection_strategy = selection_strategy;
data_info.type_distribution = type_distribution;
data_info.load_summary = load_summary;
data_info.load_time = now;

% 统计信息
if ~isempty(selected_data)
    qualities = zeros(length(selected_data), 1);
    excitation_types = {};
    
    for i = 1:length(selected_data)
        qualities(i) = selected_data{i}.quick_stats.quality_score;
        excitation_types{i} = selected_data{i}.quick_stats.excitation_name;
    end
    
    data_info.quality_stats = struct();
    data_info.quality_stats.mean_quality = mean(qualities);
    data_info.quality_stats.min_quality = min(qualities);
    data_info.quality_stats.max_quality = max(qualities);
    data_info.quality_stats.std_quality = std(qualities);
    
    [unique_types, ~, idx] = unique(excitation_types);
    type_counts = accumarray(idx, 1);
    data_info.excitation_summary = struct();
    for i = 1:length(unique_types)
        data_info.excitation_summary.(matlab.lang.makeValidName(unique_types{i})) = type_counts(i);
    end
end

fprintf('  数据加载完成: %d 组数据, 平均质量 %.2f\n', ...
        data_info.total_files, data_info.quality_stats.mean_quality);

end

function [optimized_files, type_distribution] = optimize_selection_by_excitation_type(selected_files)
% 按激励类型优化选择

% 按激励类型分组
type_groups = struct();
for i = 1:length(selected_files)
    try
        % 从文件名提取激励类型
        [~, filename, ~] = fileparts(selected_files{i});
        
        if contains(filename, 'Z_Axis')
            type_name = 'Z_Axis';
        elseif contains(filename, 'Pitch')
            type_name = 'Pitch';
        elseif contains(filename, 'Armz_Joint')
            type_name = 'Armz_Joint';
        elseif contains(filename, 'Army_Joint')
            type_name = 'Army_Joint';
        else
            type_name = 'Unknown';
        end
        
        if ~isfield(type_groups, type_name)
            type_groups.(type_name) = {};
        end
        type_groups.(type_name){end+1} = selected_files{i};
        
    catch
        continue;
    end
end

% 每种类型选择最多6个，保证多样性
optimized_files = {};
type_distribution = struct();

type_names = fieldnames(type_groups);
for i = 1:length(type_names)
    type_name = type_names{i};
    files_of_type = type_groups.(type_name);
    
    % 每种类型最多选6个
    selected_count = min(6, length(files_of_type));
    
    % 优先选择推荐文件
    recommended_of_type = {};
    stable_of_type = {};
    
    for j = 1:length(files_of_type)
        if contains(files_of_type{j}, '_RECOMMENDED_')
            recommended_of_type{end+1} = files_of_type{j};
        else
            stable_of_type{end+1} = files_of_type{j};
        end
    end
    
    % 组合选择
    type_selected = {};
    if length(recommended_of_type) >= selected_count
        type_selected = recommended_of_type(1:selected_count);
    else
        type_selected = [recommended_of_type, ...
                        stable_of_type(1:min(selected_count-length(recommended_of_type), length(stable_of_type)))];
    end
    
    optimized_files = [optimized_files, type_selected];
    type_distribution.(type_name) = length(type_selected);
    
    fprintf('    %s: 选择 %d/%d 个文件\n', type_name, length(type_selected), length(files_of_type));
end

end