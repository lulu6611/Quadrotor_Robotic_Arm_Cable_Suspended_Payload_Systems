function [M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED] = perform_grouped_mcg_estimation(q_data, q_dot_data, q_ddot_data, u_data, quality_assessment)
% 分组MCG估计策略

[N, n_dof] = size(q_data);
n_inputs = size(u_data, 2);

fprintf('      基于 %d 样本执行分组MCG估计\n', N);

%% ====== 初始化MCG矩阵 ======
M_IDENTIFIED = zeros(n_dof, n_dof);
G_IDENTIFIED = zeros(n_dof, 1);
C_IDENTIFIED = zeros(n_dof, n_dof);

%% ====== 组1：位置动力学 (DOF 1-3) ======
fprintf('        组1: FLU位置动力学...\n');
[M_pos, G_pos] = estimate_position_dynamics_improved(q_data(:,1:3), q_dot_data(:,1:3), q_ddot_data(:,1:3), u_data(:,1:3), quality_assessment);
M_IDENTIFIED(1:3, 1:3) = M_pos;
G_IDENTIFIED(1:3) = G_pos;

%% ====== 组2：姿态动力学 (DOF 4-6) ======
fprintf('        组2: 姿态动力学...\n');
[M_att, G_att] = estimate_attitude_dynamics_improved(q_data(:,4:6), q_dot_data(:,4:6), q_ddot_data(:,4:6), u_data(:,4:6), quality_assessment);
M_IDENTIFIED(4:6, 4:6) = M_att;
G_IDENTIFIED(4:6) = G_att;

%% ====== 组3：关节动力学 (DOF 7-8) ======
fprintf('        组3: 关节动力学...\n');
[M_joint, G_joint] = estimate_joint_dynamics_improved(q_data(:,7:8), q_dot_data(:,7:8), q_ddot_data(:,7:8), u_data(:,7:8), quality_assessment);
M_IDENTIFIED(7:8, 7:8) = M_joint;
G_IDENTIFIED(7:8) = G_joint;

%% ====== 组4：绳索动力学 (DOF 9-10) ======
if n_dof >= 10
    fprintf('        组4: 绳索动力学...\n');
    [M_rope, G_rope] = estimate_rope_dynamics_simplified(q_data, q_dot_data, q_ddot_data);
    M_IDENTIFIED(9:10, 9:10) = M_rope;
    G_IDENTIFIED(9:10) = G_rope;
end

%% ====== 主要耦合项估计 ======
fprintf('        耦合项估计...\n');
M_IDENTIFIED = add_main_coupling_terms(M_IDENTIFIED, q_data, q_dot_data, q_ddot_data, u_data, quality_assessment);

%% ====== 确保物理性质 ======
M_IDENTIFIED = ensure_physical_properties(M_IDENTIFIED);

%% ====== 科氏矩阵估计 ======
fprintf('        科氏矩阵估计...\n');
C_IDENTIFIED = estimate_coriolis_matrix_from_data(M_IDENTIFIED, q_data, q_dot_data, q_ddot_data);

fprintf('       分组MCG估计完成\n');

end

function [M_pos, G_pos] = estimate_position_dynamics_improved(q_pos, q_dot_pos, q_ddot_pos, u_pos, quality_assessment)
% 改进的位置动力学估计

[N, ~] = size(q_pos);
M_pos = zeros(3, 3);
G_pos = zeros(3, 1);

% 物理参考值
m_total_ref = 2.019; % kg
g = 9.81; % m/s²

dof_names = {'X(前)', 'Y(左)', 'Z(上)'};

for i = 1:3
    try
        % 获取DOF质量评估
        if isfield(quality_assessment, 'dof_quality') && isfield(quality_assessment.dof_quality, sprintf('dof_%d', i))
            dof_quality = quality_assessment.dof_quality.(sprintf('dof_%d', i));
            data_weight = dof_quality.overall_score;
        else
            data_weight = 0.5;
        end
        
        % 数据有效性检查
        valid_mask = ~(isnan(q_ddot_pos(:, i)) | isnan(u_pos(:, i)) | isinf(q_ddot_pos(:, i)) | isinf(u_pos(:, i)));
        valid_ratio = sum(valid_mask) / N;
        
        if valid_ratio < 0.6
            fprintf('          %s: 数据不足(%.1f%%)，使用物理默认值\n', dof_names{i}, valid_ratio*100);
            M_pos(i, i) = m_total_ref;
            if i == 3
                G_pos(i) = -m_total_ref * g;
            end
            continue;
        end
        
        % 有效数据
        q_ddot_valid = q_ddot_pos(valid_mask, i);
        u_valid = u_pos(valid_mask, i);
        N_valid = length(q_ddot_valid);
        
        % 构建设计矩阵 (包含重力项)
        X = [q_ddot_valid, ones(N_valid, 1)];
        y = u_valid;
        
        % 数据预处理：去除异常值
        outlier_threshold = 3;
        residuals = abs(y - median(y));
        outlier_mask = residuals < outlier_threshold * mad(y, 1);
        
        X_clean = X(outlier_mask, :);
        y_clean = y(outlier_mask);
        
        if length(y_clean) < 20
            fprintf('          %s: 清洁数据不足，使用默认值\n', dof_names{i});
            M_pos(i, i) = m_total_ref;
            if i == 3
                G_pos(i) = -m_total_ref * g;
            end
            continue;
        end
        
        % 加权最小二乘估计
        weights = 1 ./ (1 + 0.2 * abs(y_clean - median(y_clean)));
        W = diag(weights);
        
        % 正则化参数
        lambda = 1e-4 * data_weight;
        
        % 求解
        XtWX = X_clean' * W * X_clean + lambda * eye(2);
        XtWy = X_clean' * W * y_clean;
        
        if rcond(XtWX) > 1e-12
            theta = XtWX \ XtWy;
        else
            theta = pinv(XtWX) * XtWy;
        end
        
        % 提取参数
        mass_estimate = theta(1);
        gravity_estimate = theta(2);
        
        % 物理约束和验证
        mass_lower = 0.5 * m_total_ref;
        mass_upper = 3.0 * m_total_ref;
        
        if mass_estimate >= mass_lower && mass_estimate <= mass_upper
            M_pos(i, i) = mass_estimate;
        else
            M_pos(i, i) = m_total_ref;
            fprintf('          %s: 质量估计超出范围(%.3f)，使用默认值\n', dof_names{i}, mass_estimate);
        end
        
        % 重力项验证 (Z轴应该有显著重力分量)
        if i == 3
            expected_gravity = -m_total_ref * g;
            if abs(gravity_estimate - expected_gravity) < abs(expected_gravity) * 0.3
                G_pos(i) = gravity_estimate;
            else
                G_pos(i) = expected_gravity;
                fprintf('          %s: 重力估计偏差过大(%.3f vs %.3f)，使用默认值\n', dof_names{i}, gravity_estimate, expected_gravity);
            end
        else
            % X, Y方向重力分量应该很小
            if abs(gravity_estimate) < 5.0
                G_pos(i) = gravity_estimate;
            else
                G_pos(i) = 0;
            end
        end
        
        % 拟合质量评估
        y_pred = X_clean * theta;
        ssr = sum((y_clean - y_pred).^2);
        sst = sum((y_clean - mean(y_clean)).^2);
        r_squared = 1 - ssr/sst;
        
        correlation = corrcoef(y_pred, y_clean);
        corr_val = correlation(1, 2);
        
        fprintf('          %s: M=%.3f kg, G=%.3f N, R²=%.3f, 相关性=%.3f\n', ...
                dof_names{i}, M_pos(i,i), G_pos(i), r_squared, corr_val);
        
    catch ME
        fprintf('          %s: 估计失败(%s)，使用默认值\n', dof_names{i}, ME.message);
        M_pos(i, i) = m_total_ref;
        if i == 3
            G_pos(i) = -m_total_ref * g;
        end
    end
end

end

function [M_att, G_att] = estimate_attitude_dynamics_improved(q_att, q_dot_att, q_ddot_att, u_att, quality_assessment)
% 改进的姿态动力学估计

[N, ~] = size(q_att);
M_att = zeros(3, 3);
G_att = zeros(3, 1);

% 物理参考惯量
I_ref = [0.02054, 0.01050, 0.01039]; % kg⋅m² [Ixx, Iyy, Izz]
att_names = {'Roll', 'Pitch', 'Yaw'};

for i = 1:3
    try
        % 获取DOF质量评估
        dof_idx = i + 3; % 姿态DOF从4开始
        if isfield(quality_assessment, 'dof_quality') && isfield(quality_assessment.dof_quality, sprintf('dof_%d', dof_idx))
            dof_quality = quality_assessment.dof_quality.(sprintf('dof_%d', dof_idx));
            data_weight = dof_quality.overall_score;
        else
            data_weight = 0.3; % 姿态通常难以辨识
        end
        
        % 数据预处理
        valid_mask = ~(isnan(q_ddot_att(:, i)) | isnan(u_att(:, i)) | isinf(q_ddot_att(:, i)) | isinf(u_att(:, i)));
        valid_ratio = sum(valid_mask) / N;
        
        if valid_ratio < 0.5 || data_weight < 0.2
            fprintf('          %s: 数据质量不足，使用物理默认值\n', att_names{i});
            M_att(i, i) = I_ref(i);
            continue;
        end
        
        % 有效数据
        q_ddot_valid = q_ddot_att(valid_mask, i);
        u_valid = u_att(valid_mask, i);
        
        % 去除异常值
        outlier_threshold = 3;
        residuals = abs(u_valid - median(u_valid));
        outlier_mask = residuals < outlier_threshold * mad(u_valid, 1);
        
        q_ddot_clean = q_ddot_valid(outlier_mask);
        u_clean = u_valid(outlier_mask);
        
        if length(u_clean) < 15
            fprintf('          %s: 清洁数据不足，使用默认值\n', att_names{i});
            M_att(i, i) = I_ref(i);
            continue;
        end
        
        % 构建设计矩阵
        X = [q_ddot_clean, ones(length(q_ddot_clean), 1)];
        y = u_clean;
        
        % 加权最小二乘 (姿态控制精度较低，权重较均匀)
        weights = ones(length(y), 1);
        W = diag(weights);
        
        % 较强的正则化
        lambda = 1e-3;
        
        % 求解
        XtWX = X' * W * X + lambda * eye(2);
        XtWy = X' * W * y;
        
        theta = XtWX \ XtWy;
        
        % 提取参数
        inertia_estimate = theta(1);
        gravity_torque = theta(2);
        
        % 物理约束 (惯量范围)
        inertia_lower = 0.1 * I_ref(i);
        inertia_upper = 10.0 * I_ref(i);
        
        if inertia_estimate >= inertia_lower && inertia_estimate <= inertia_upper
            M_att(i, i) = inertia_estimate;
        else
            M_att(i, i) = I_ref(i);
            fprintf('          %s: 惯量估计超出范围(%.6f)，使用默认值\n', att_names{i}, inertia_estimate);
        end
        
        % 重力力矩 (姿态系统应该很小)
        if abs(gravity_torque) < 0.5
            G_att(i) = gravity_torque;
        else
            G_att(i) = 0;
        end
        
        % 性能评估
        y_pred = X * theta;
        correlation = corrcoef(y_pred, y);
        corr_val = correlation(1, 2);
        
        fprintf('          %s: I=%.6f kg⋅m², G=%.6f N⋅m, 相关性=%.3f\n', ...
                att_names{i}, M_att(i,i), G_att(i), corr_val);
        
    catch ME
        fprintf('          %s: 估计失败(%s)，使用默认值\n', att_names{i}, ME.message);
        M_att(i, i) = I_ref(i);
    end
end

end

function [M_joint, G_joint] = estimate_joint_dynamics_improved(q_joint, q_dot_joint, q_ddot_joint, u_joint, quality_assessment)
% 改进的关节动力学估计

[N, ~] = size(q_joint);
M_joint = zeros(2, 2);
G_joint = zeros(2, 1);

% 关节惯量参考值
I_joint_ref = [0.005, 0.008]; % kg⋅m² 估计值
joint_names = {'Armz', 'Army'};

for i = 1:2
    try
        % 获取DOF质量评估
        dof_idx = i + 6; % 关节DOF从7开始
        if isfield(quality_assessment, 'dof_quality') && isfield(quality_assessment.dof_quality, sprintf('dof_%d', dof_idx))
            dof_quality = quality_assessment.dof_quality.(sprintf('dof_%d', dof_idx));
            data_weight = dof_quality.overall_score;
        else
            data_weight = 0.6; % 关节通常有较好的辨识性
        end
        
        % 数据预处理
        valid_mask = ~(isnan(q_ddot_joint(:, i)) | isnan(u_joint(:, i)) | isinf(q_ddot_joint(:, i)) | isinf(u_joint(:, i)));
        valid_ratio = sum(valid_mask) / N;
        
        if valid_ratio < 0.4
            fprintf('          %s: 数据不足(%.1f%%)，使用默认值\n', joint_names{i}, valid_ratio*100);
            M_joint(i, i) = I_joint_ref(i);
            continue;
        end
        
        % 有效数据
        q_ddot_valid = q_ddot_joint(valid_mask, i);
        u_valid = u_joint(valid_mask, i);
        
        % 数据质量检查
        if std(q_ddot_valid) < 1e-6 || std(u_valid) < 1e-6
            fprintf('          %s: 数据方差过小，使用默认值\n', joint_names{i});
            M_joint(i, i) = I_joint_ref(i);
            continue;
        end
        
        % 去除异常值
        outlier_threshold = 3;
        combined_data = [q_ddot_valid, u_valid];
        outlier_mask = all(abs(combined_data - median(combined_data)) < outlier_threshold * mad(combined_data, 1), 2);
        
        q_ddot_clean = q_ddot_valid(outlier_mask);
        u_clean = u_valid(outlier_mask);
        
        if length(u_clean) < 10
            fprintf('          %s: 清洁数据不足，使用默认值\n', joint_names{i});
            M_joint(i, i) = I_joint_ref(i);
            continue;
        end
        
        % 构建设计矩阵
        X = [q_ddot_clean, ones(length(q_ddot_clean), 1)];
        y = u_clean;
        
        % 加权最小二乘
        weights = 1 ./ (1 + 0.1 * abs(y - median(y)));
        W = diag(weights);
        
        % 轻度正则化
        lambda = 1e-5 * data_weight;
        
        % 求解
        XtWX = X' * W * X + lambda * eye(2);
        XtWy = X' * W * y;
        
        theta = XtWX \ XtWy;
        
        % 提取参数
        inertia_estimate = theta(1);
        gravity_torque = theta(2);
        
        % 物理约束
        inertia_lower = 1e-6;
        inertia_upper = 0.1;
        
        if inertia_estimate >= inertia_lower && inertia_estimate <= inertia_upper
            M_joint(i, i) = inertia_estimate;
        else
            M_joint(i, i) = max(I_joint_ref(i), 1e-6);
            fprintf('          %s: 惯量估计超出范围(%.6f)，使用约束值\n', joint_names{i}, inertia_estimate);
        end
        
        % 重力力矩 (可能有一定的重力影响)
        if abs(gravity_torque) < 1.0
            G_joint(i) = gravity_torque;
        else
            G_joint(i) = 0;
        end
        
        % 性能评估
        y_pred = X * theta;
        correlation = corrcoef(y_pred, y);
        corr_val = correlation(1, 2);
        
        fprintf('          %s: I=%.6f kg⋅m², G=%.6f N⋅m, 相关性=%.3f\n', ...
                joint_names{i}, M_joint(i,i), G_joint(i), corr_val);
        
    catch ME
        fprintf('          %s: 估计失败(%s)，使用默认值\n', joint_names{i}, ME.message);
        M_joint(i, i) = I_joint_ref(i);
    end
end

end

function [M_rope, G_rope] = estimate_rope_dynamics_simplified(q_data, q_dot_data, q_ddot_data)
% 简化的绳索动力学估计

M_rope = zeros(2, 2);
G_rope = zeros(2, 1);

% 绳索系统参数估计
m_load = 0.24; % kg 载荷质量
L_rope = 1.15; % m 绳长

% 简化的绳索惯量
I_rope_equiv = m_load * L_rope^2;

M_rope(1, 1) = I_rope_equiv;
M_rope(2, 2) = I_rope_equiv;

% 绳索重力项设为零（简化处理）
G_rope(1) = 0;
G_rope(2) = 0;

fprintf('          绳索: I=%.6f kg⋅m² (估计值)\n', I_rope_equiv);

end

function M_coupled = add_main_coupling_terms(M, q_data, q_dot_data, q_ddot_data, u_data, quality_assessment)
% 添加主要耦合项

M_coupled = M;

% 基于物理直觉的主要耦合项
fprintf('          添加物理耦合项...\n');

% 1. Z位置与Roll/Pitch的耦合（重心偏移效应）
coupling_scale_pos_att = 0.08;
M_coupled(3, 4) = coupling_scale_pos_att * sqrt(M(3, 3) * M(4, 4));
M_coupled(3, 5) = coupling_scale_pos_att * sqrt(M(3, 3) * M(5, 5));
M_coupled(4, 3) = M_coupled(3, 4);
M_coupled(5, 3) = M_coupled(3, 5);

% 2. Roll与Armz关节的耦合
coupling_scale_roll_armz = 0.05;
M_coupled(4, 7) = coupling_scale_roll_armz * sqrt(M(4, 4) * M(7, 7));
M_coupled(7, 4) = M_coupled(4, 7);

% 3. Pitch与Army关节的耦合
coupling_scale_pitch_army = 0.05;
M_coupled(5, 8) = coupling_scale_pitch_army * sqrt(M(5, 5) * M(8, 8));
M_coupled(8, 5) = M_coupled(5, 8);

% 4. 关节间的弱耦合
coupling_scale_joint = 0.02;
M_coupled(7, 8) = coupling_scale_joint * sqrt(M(7, 7) * M(8, 8));
M_coupled(8, 7) = M_coupled(7, 8);

% 确保对称性
M_coupled = 0.5 * (M_coupled + M_coupled');

coupling_count = sum(sum(abs(M_coupled - M) > 1e-8)) / 2; % 除以2因为对称
fprintf('           添加了 %d 个耦合项\n', coupling_count);

end

function M_proper = ensure_physical_properties(M)
% 确保质量矩阵的物理性质

% 1. 强制对称性
M_proper = 0.5 * (M + M');

% 2. 确保正定性
[V, D] = eig(M_proper);
eigenvals = diag(D);

% 记录原始特征值
original_min = min(eigenvals);
original_max = max(eigenvals);

% 约束特征值
min_eigenval = 1e-8;
max_eigenval = 1e6;

corrected_count = 0;
if any(eigenvals < min_eigenval)
    corrected_count = sum(eigenvals < min_eigenval);
    eigenvals(eigenvals < min_eigenval) = min_eigenval;
end

if any(eigenvals > max_eigenval)
    corrected_count = corrected_count + sum(eigenvals > max_eigenval);
    eigenvals(eigenvals > max_eigenval) = max_eigenval;
end

% 重构矩阵
M_proper = V * diag(eigenvals) * V';

% 再次强制对称性
M_proper = 0.5 * (M_proper + M_proper');

if corrected_count > 0
    fprintf('           修正了 %d 个特征值 [%.2e, %.2e] → [%.2e, %.2e]\n', ...
            corrected_count, original_min, original_max, min(eigenvals), max(eigenvals));
end

% 条件数检查
cond_num = cond(M_proper);
if cond_num > 1e10
    fprintf('           警告：条件数过大 %.2e\n', cond_num);
end

end

function C_est = estimate_coriolis_matrix_from_data(M, q_data, q_dot_data, q_ddot_data)
% 从数据估计科氏矩阵

[N, n_dof] = size(q_data);
C_est = zeros(n_dof, n_dof);

% 简化的科氏矩阵估计
% 基于角速度和角加速度的经验关系

fprintf('          科氏矩阵基于角速度效应...\n');

% 主要来自姿态运动的科氏力
for i = 4:6  % 姿态DOF
    for j = 4:6
        if i ~= j
            % 基于平均角速度的科氏项
            omega_i = mean(abs(q_dot_data(:, i)));
            omega_j = mean(abs(q_dot_data(:, j)));
            
            if omega_i > 1e-3 && omega_j > 1e-3
                % 科氏系数
                coriolis_coeff = 0.05 * sqrt(M(i, i) * M(j, j)) * min(omega_i, omega_j);
                C_est(i, j) = coriolis_coeff;
            end
        end
    end
end

% 确保反对称性
C_est = 0.5 * (C_est - C_est');

symmetry_error = norm(C_est + C_est', 'fro');
if symmetry_error < 1e-8
    fprintf('           科氏矩阵反对称性验证通过\n');
else
    fprintf('           科氏矩阵反对称性误差: %.2e\n', symmetry_error);
end

end