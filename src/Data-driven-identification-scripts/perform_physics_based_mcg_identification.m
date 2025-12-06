function [M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results] = perform_physics_based_mcg_identification(selected_data, data_info)
% 基于物理正确性的MCG辨识

fprintf('  开始基于物理正确性的MCG系统辨识...\n');

if isempty(selected_data)
    error('没有可用的数据进行辨识');
end

%% ====== 数据预处理与合并 ======
fprintf('    数据预处理与合并...\n');
[combined_data, processing_info] = preprocess_and_combine_collected_data(selected_data, data_info);

%% ====== FLU坐标系数据提取 ======
fprintf('   FLU坐标系数据提取...\n');
[q_data, q_dot_data, q_ddot_data, u_data, coordinate_info] = extract_flu_coordinate_data(combined_data, processing_info);

%% ====== 数据质量评估 ======
fprintf('     数据质量评估...\n');
quality_assessment = assess_combined_data_quality(q_data, q_dot_data, q_ddot_data, u_data);

%% ====== 分组辨识策略 ======
fprintf('     执行分组辨识策略...\n');
[M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED] = perform_grouped_mcg_estimation(q_data, q_dot_data, q_ddot_data, u_data, quality_assessment);

%% ====== 结果验证 ======
fprintf('     辨识结果验证...\n');
validation_results = validate_physics_identification_results(M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, q_data, q_dot_data, q_ddot_data, u_data);

%% ====== 编译辨识结果 ======
id_results = struct();
id_results.coordinate_system = 'FLU';
id_results.n_dof = size(M_IDENTIFIED, 1);
id_results.n_samples = size(q_data, 1);
id_results.data_sources = length(selected_data);
id_results.processing_info = processing_info;
id_results.coordinate_info = coordinate_info;
id_results.quality_assessment = quality_assessment;
id_results.validation_results = validation_results;
id_results.analysis = validation_results; % 保持兼容性

fprintf('   物理MCG辨识完成\n');
fprintf('     平均相关性: %.3f\n', validation_results.mean_correlation);
fprintf('     整体性能: %s\n', validation_results.overall_performance);

end

function [combined_data, processing_info] = preprocess_and_combine_collected_data(selected_data, data_info)
% 预处理并合并采集的数据

processing_info = struct();
processing_info.original_count = length(selected_data);

% 提取所有数据的公共时间范围
min_duration = inf;
sample_rates = [];

for i = 1:length(selected_data)
    data = selected_data{i};
    if isfield(data, 'simulation') && isfield(data.simulation, 'time')
        duration = data.simulation.time(end) - data.simulation.time(1);
        min_duration = min(min_duration, duration);
        dt = mean(diff(data.simulation.time));
        sample_rates(end+1) = 1/dt;
    end
end

% 统一采样率和时间范围
target_sample_rate = median(sample_rates);
target_duration = min(min_duration, 20); % 最多使用20秒数据
processing_info.target_sample_rate = target_sample_rate;
processing_info.target_duration = target_duration;

fprintf('      目标采样率: %.1f Hz, 目标时长: %.1f 秒\n', target_sample_rate, target_duration);

% 合并数据
combined_data = struct();
combined_data.datasets = {};

valid_count = 0;
for i = 1:length(selected_data)
    try
        data = selected_data{i};
        
        % 重采样到统一时间基准
        resampled_data = resample_dataset_to_standard(data, target_sample_rate, target_duration);
        
        if ~isempty(resampled_data)
            combined_data.datasets{end+1} = resampled_data;
            valid_count = valid_count + 1;
        end
        
    catch ME
        fprintf('       数据集 %d 处理失败: %s\n', i, ME.message);
    end
end

processing_info.valid_count = valid_count;
processing_info.success_rate = valid_count / processing_info.original_count;

fprintf('       成功处理 %d/%d 个数据集\n', valid_count, processing_info.original_count);

end

function resampled_data = resample_dataset_to_standard(data, target_fs, target_duration)
% 将数据集重采样到标准格式

resampled_data = struct();

if ~isfield(data, 'simulation')
    resampled_data = [];
    return;
end

sim_data = data.simulation;
original_time = sim_data.time;

% 计算目标时间向量
target_dt = 1/target_fs;
target_time = 0:target_dt:target_duration;

% 确保不超出原始数据范围
max_time = min(target_duration, original_time(end));
target_time = target_time(target_time <= max_time);

try
    % 重采样所有信号
    resampled_data.time = target_time';
    resampled_data.meta = data.meta;
    resampled_data.quick_stats = data.quick_stats;
    
    % 位置和速度数据
    if isfield(sim_data, 'imu_pos') && size(sim_data.imu_pos, 1) == length(original_time)
        resampled_data.imu_pos = interp1(original_time, sim_data.imu_pos, target_time, 'linear');
        resampled_data.world_vel = interp1(original_time, sim_data.world_vel, target_time, 'linear');
    else
        resampled_data = [];
        return;
    end
    
    % 姿态数据
    if isfield(sim_data, 'body_quat') && size(sim_data.body_quat, 1) == length(original_time)
        resampled_data.body_quat = interp1(original_time, sim_data.body_quat, target_time, 'linear');
        resampled_data.body_gyro = interp1(original_time, sim_data.body_gyro, target_time, 'linear');
    end
    
    % 关节数据
    if isfield(sim_data, 'armz_pos') && length(sim_data.armz_pos) == length(original_time)
        resampled_data.armz_pos = interp1(original_time, sim_data.armz_pos, target_time, 'linear');
        resampled_data.army_pos = interp1(original_time, sim_data.army_pos, target_time, 'linear');
        resampled_data.armz_vel = interp1(original_time, sim_data.armz_vel, target_time, 'linear');
        resampled_data.army_vel = interp1(original_time, sim_data.army_vel, target_time, 'linear');
    end
    
    % 控制信号
    if isfield(sim_data, 'f_total') && length(sim_data.f_total) == length(original_time)
        resampled_data.f_total = interp1(original_time, sim_data.f_total, target_time, 'linear');
        resampled_data.tau_roll = interp1(original_time, sim_data.tau_roll, target_time, 'linear');
        resampled_data.tau_pitch = interp1(original_time, sim_data.tau_pitch, target_time, 'linear');
        resampled_data.tau_yaw = interp1(original_time, sim_data.tau_yaw, target_time, 'linear');
        resampled_data.tau_armz = interp1(original_time, sim_data.tau_armz, target_time, 'linear');
        resampled_data.tau_army = interp1(original_time, sim_data.tau_army, target_time, 'linear');
    end
    
catch ME
    fprintf('        重采样失败: %s\n', ME.message);
    resampled_data = [];
end

end

function [q_data, q_dot_data, q_ddot_data, u_data, coordinate_info] = extract_flu_coordinate_data(combined_data, processing_info)
% 提取FLU坐标系数据

if isempty(combined_data.datasets)
    error(' 没有有效的数据集');
end

% 选择质量最好的数据集或合并多个数据集
if length(combined_data.datasets) == 1
    fprintf('      使用单一数据集\n');
    selected_dataset = combined_data.datasets{1};
    combination_method = 'single';
else
    fprintf('      合并 %d 个数据集\n', length(combined_data.datasets));
    selected_dataset = combine_multiple_datasets(combined_data.datasets);
    combination_method = 'merged';
end

dt = mean(diff(selected_dataset.time));
N = length(selected_dataset.time);

fprintf('      数据长度: %d 样本, 采样间隔: %.4f 秒\n', N, dt);

%% ====== 构建状态向量 ======
q_data = zeros(N, 10);

% 位置 (FLU坐标系: Forward, Left, Up)
q_data(:, 1:3) = selected_dataset.imu_pos;

% 姿态 (从四元数转换为欧拉角)
if isfield(selected_dataset, 'body_quat') && ~isempty(selected_dataset.body_quat)
    euler_angles = quat2euler_flu_corrected(selected_dataset.body_quat);
    q_data(:, 4:6) = euler_angles;
end

% 关节角度
if isfield(selected_dataset, 'armz_pos') && ~isempty(selected_dataset.armz_pos)
    q_data(:, 7) = selected_dataset.armz_pos;
    q_data(:, 8) = selected_dataset.army_pos;
end

% 绳索摆角 (简化处理，设为零或从其他数据估计)
q_data(:, 9:10) = zeros(N, 2);

%% ====== 构建速度向量 ======
q_dot_data = zeros(N, 10);

% 位置速度
q_dot_data(:, 1:3) = selected_dataset.world_vel;

% 角速度
if isfield(selected_dataset, 'body_gyro') && ~isempty(selected_dataset.body_gyro)
    q_dot_data(:, 4:6) = selected_dataset.body_gyro;
else
    % 从姿态数值微分
    for i = 4:6
        q_dot_data(:, i) = gradient(q_data(:, i), dt);
    end
end

% 关节角速度
if isfield(selected_dataset, 'armz_vel') && ~isempty(selected_dataset.armz_vel)
    q_dot_data(:, 7) = selected_dataset.armz_vel;
    q_dot_data(:, 8) = selected_dataset.army_vel;
else
    q_dot_data(:, 7) = gradient(q_data(:, 7), dt);
    q_dot_data(:, 8) = gradient(q_data(:, 8), dt);
end

% 绳索角速度
q_dot_data(:, 9:10) = zeros(N, 2);

%% ====== 计算加速度 ======
q_ddot_data = zeros(N, 10);
for i = 1:10
    q_ddot_data(:, i) = gradient(q_dot_data(:, i), dt);
end

%% ====== 计算控制输入 (基于物理模型) ======
u_data = compute_control_forces_flu_from_collected_data(selected_dataset, q_data, dt);

%% ====== 坐标信息 ======
coordinate_info = struct();
coordinate_info.coordinate_system = 'FLU';
coordinate_info.combination_method = combination_method;
coordinate_info.sample_rate = 1/dt;
coordinate_info.duration = selected_dataset.time(end);
coordinate_info.n_samples = N;
coordinate_info.data_quality = evaluate_extracted_data_quality(q_data, q_dot_data, q_ddot_data, u_data);

fprintf('       FLU坐标数据提取完成\n');

end

function euler_angles = quat2euler_flu_corrected(quat_data)
% 四元数到欧拉角转换 (FLU坐标系)

N = size(quat_data, 1);
euler_angles = zeros(N, 3);

for i = 1:N
    q = quat_data(i, :);
    if size(q, 2) >= 4
        w = q(1); x = q(2); y = q(3); z = q(4);
        
        % 归一化
        norm_q = sqrt(w^2 + x^2 + y^2 + z^2);
        if norm_q > 1e-6
            w = w/norm_q; x = x/norm_q; y = y/norm_q; z = z/norm_q;
        end
        
        % FLU坐标系欧拉角转换
        % Roll (绕X轴)
        sinr_cosp = 2 * (w * x + y * z);
        cosr_cosp = 1 - 2 * (x * x + y * y);
        roll = atan2(sinr_cosp, cosr_cosp);
        
        % Pitch (绕Y轴)
        sinp = 2 * (w * y - z * x);
        sinp = max(-1, min(1, sinp));
        pitch = asin(sinp);
        
        % Yaw (绕Z轴)
        siny_cosp = 2 * (w * z + x * y);
        cosy_cosp = 1 - 2 * (y * y + z * z);
        yaw = atan2(siny_cosp, cosy_cosp);
        
        euler_angles(i, :) = [roll, pitch, yaw];
    end
end
end

function u_data = compute_control_forces_flu_from_collected_data(dataset, q_data, dt)
% 从采集数据计算FLU坐标系控制力

N = size(q_data, 1);
u_data = zeros(N, 8);

% 提取姿态角
euler_angles = q_data(:, 4:6);
roll = euler_angles(:, 1);
pitch = euler_angles(:, 2);
yaw = euler_angles(:, 3);

% 提取控制命令
if isfield(dataset, 'f_total') && ~isempty(dataset.f_total)
    f_total = dataset.f_total;
    tau_roll = dataset.tau_roll;
    tau_pitch = dataset.tau_pitch;
    tau_yaw = dataset.tau_yaw;
    tau_armz = dataset.tau_armz;
    tau_army = dataset.tau_army;
else
    % 如果没有控制数据，使用估计值
    f_total = ones(N, 1) * 19.8; % 约2kg的悬停推力
    tau_roll = zeros(N, 1);
    tau_pitch = zeros(N, 1);
    tau_yaw = zeros(N, 1);
    tau_armz = zeros(N, 1);
    tau_army = zeros(N, 1);
end

% 计算世界坐标系力
for i = 1:N
    % FLU坐标系旋转矩阵
    cy = cos(yaw(i)); sy = sin(yaw(i));
    cp = cos(pitch(i)); sp = sin(pitch(i));
    cr = cos(roll(i)); sr = sin(roll(i));
    
    R_body_to_world = [
        cy*cp,  -sy*cr + cy*sp*sr,   sy*sr + cy*sp*cr;
        sy*cp,   cy*cr + sy*sp*sr,  -cy*sr + sy*sp*cr;
        -sp,     cp*sr,              cp*cr
    ];
    
    % 机体系推力
    f_body = [0; 0; f_total(i)];
    
    % 转换到世界系
    f_world = R_body_to_world * f_body;
    
    % 基于之前的符号修正经验，使用负号
    u_data(i, 1) = -f_world(1);    % X力（前向）
    u_data(i, 2) = -f_world(2);    % Y力（左向）
    u_data(i, 3) = -f_world(3);    % Z力（向上）
    u_data(i, 4) = -tau_roll(i);   % Roll力矩
    u_data(i, 5) = -tau_pitch(i);  % Pitch力矩
    u_data(i, 6) = -tau_yaw(i);    % Yaw力矩
    u_data(i, 7) = -tau_armz(i);   % Armz关节力矩
    u_data(i, 8) = -tau_army(i);   % Army关节力矩
end

end

function combined_dataset = combine_multiple_datasets(datasets)
% 合并多个数据集

if length(datasets) == 1
    combined_dataset = datasets{1};
    return;
end

% 找到公共的最短时间
min_length = inf;
for i = 1:length(datasets)
    min_length = min(min_length, length(datasets{i}.time));
end

% 截取到相同长度
for i = 1:length(datasets)
    datasets{i} = truncate_dataset(datasets{i}, min_length);
end

% 合并策略：使用加权平均
combined_dataset = datasets{1}; % 基础数据集

% 对数值数据进行加权平均
numeric_fields = {'imu_pos', 'world_vel', 'body_quat', 'body_gyro', 'armz_pos', 'army_pos', 'armz_vel', 'army_vel', ...
                  'f_total', 'tau_roll', 'tau_pitch', 'tau_yaw', 'tau_armz', 'tau_army'};

for field_idx = 1:length(numeric_fields)
    field_name = numeric_fields{field_idx};
    
    if isfield(combined_dataset, field_name)
        data_sum = combined_dataset.(field_name);
        weight_sum = 1;
        
        for i = 2:length(datasets)
            if isfield(datasets{i}, field_name) && ~isempty(datasets{i}.(field_name))
                data_sum = data_sum + datasets{i}.(field_name);
                weight_sum = weight_sum + 1;
            end
        end
        
        combined_dataset.(field_name) = data_sum / weight_sum;
    end
end

fprintf('         合并了 %d 个数据集\n', length(datasets));

end

function truncated_dataset = truncate_dataset(dataset, target_length)
% 截取数据集到指定长度

truncated_dataset = dataset;
fields = fieldnames(dataset);

for i = 1:length(fields)
    field_name = fields{i};
    data = dataset.(field_name);
    
    if isnumeric(data) && size(data, 1) > target_length
        truncated_dataset.(field_name) = data(1:target_length, :);
    end
end
end

function quality_info = evaluate_extracted_data_quality(q_data, q_dot_data, q_ddot_data, u_data)
% 评估提取数据的质量

quality_info = struct();

% 基本统计
quality_info.n_samples = size(q_data, 1);
quality_info.n_dof = size(q_data, 2);
quality_info.n_inputs = size(u_data, 2);

% 数据完整性检查
quality_info.has_nan = any(isnan(q_data(:))) || any(isnan(u_data(:)));
quality_info.has_inf = any(isinf(q_data(:))) || any(isinf(u_data(:)));

% 激励充分性评估
quality_info.excitation_adequacy = zeros(8, 1);
for i = 1:min(8, size(u_data, 2))
    if std(u_data(:, i)) > 1e-6
        quality_info.excitation_adequacy(i) = min(1.0, std(u_data(:, i)) / mean(abs(u_data(:, i)) + 1e-6));
    end
end

quality_info.overall_quality = mean(quality_info.excitation_adequacy);

end

function quality_assessment = assess_combined_data_quality(q_data, q_dot_data, q_ddot_data, u_data)
% 评估合并数据的质量

quality_assessment = struct();

[N, n_dof] = size(q_data);
n_inputs = size(u_data, 2);

fprintf('      评估 %d 样本, %d DOF, %d 输入的数据质量\n', N, n_dof, n_inputs);

% DOF质量评估
dof_names = {'X(前)', 'Y(左)', 'Z(上)', 'Roll', 'Pitch', 'Yaw', 'Armz', 'Army'};
quality_assessment.dof_quality = struct();

for i = 1:min(n_dof, n_inputs)
    % 信噪比估计
    signal_var = var(q_ddot_data(:, i));
    noise_var = var(diff(q_ddot_data(:, i))) / 2;
    snr = 10 * log10(max(signal_var / max(noise_var, 1e-12), 1e-6));
    
    % 控制激励评估
    control_std = std(u_data(:, i));
    response_std = std(q_ddot_data(:, i));
    excitation_ratio = control_std / max(response_std, 1e-12);
    
    dof_quality = struct();
    dof_quality.snr = snr;
    dof_quality.excitation_ratio = excitation_ratio;
    dof_quality.control_std = control_std;
    dof_quality.response_std = response_std;
    
    % 质量评分 (0-1)
    snr_score = min(max((snr + 10) / 30, 0), 1);  % -10dB到20dB映射到0-1
    excitation_score = min(excitation_ratio / 2, 1);
    dof_quality.overall_score = 0.6 * snr_score + 0.4 * excitation_score;
    
    quality_assessment.dof_quality.(sprintf('dof_%d', i)) = dof_quality;
    
    if i <= length(dof_names)
        fprintf('        %s: SNR=%.1fdB, 激励比=%.3f, 质量=%.2f\n', ...
                dof_names{i}, snr, excitation_ratio, dof_quality.overall_score);
    end
end

% 整体质量
scores = zeros(min(n_dof, n_inputs), 1);
for i = 1:min(n_dof, n_inputs)
    scores(i) = quality_assessment.dof_quality.(sprintf('dof_%d', i)).overall_score;
end

quality_assessment.overall_score = mean(scores);
quality_assessment.min_score = min(scores);
quality_assessment.max_score = max(scores);

if quality_assessment.overall_score > 0.7
    quality_assessment.overall_rating = 'excellent';
elseif quality_assessment.overall_score > 0.5
    quality_assessment.overall_rating = 'good';
elseif quality_assessment.overall_score > 0.3
    quality_assessment.overall_rating = 'fair';
else
    quality_assessment.overall_rating = 'poor';
end

fprintf('       数据质量评估: %.2f (%s)\n', quality_assessment.overall_score, quality_assessment.overall_rating);

end