function [M_HYBRID, G_HYBRID, C_HYBRID, fusion_results] = perform_intelligent_mcg_fusion(...
    M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results, ...
    M_THEORY, G_THEORY, C_THEORY, theory_info)
% 智能混合MCG辨识

fprintf('  执行智能混合MCG辨识...\n');

%% ====== 分析辨识质量 ======
fprintf('    分析辨识质量与理论置信度...\n');
fusion_strategy = analyze_identification_quality_and_theory_confidence(id_results, theory_info);

%% ====== 计算自适应权重 ======
fprintf('    计算自适应融合权重...\n');
fusion_weights = compute_adaptive_fusion_weights(M_IDENTIFIED, G_IDENTIFIED, M_THEORY, G_THEORY, ...
                                                id_results, theory_info, fusion_strategy);

%% ====== 执行分层融合 ======
fprintf('    执行分层融合策略...\n');
[M_HYBRID, fusion_details_M] = fuse_mass_matrix_intelligent(M_IDENTIFIED, M_THEORY, fusion_weights, fusion_strategy);
[G_HYBRID, fusion_details_G] = fuse_gravity_vector_intelligent(G_IDENTIFIED, G_THEORY, fusion_weights, fusion_strategy);
[C_HYBRID, fusion_details_C] = fuse_coriolis_matrix_intelligent(C_IDENTIFIED, C_THEORY, fusion_weights, fusion_strategy);

%% ====== 物理约束校正 ======
fprintf('    应用物理约束校正...\n');
[M_HYBRID, G_HYBRID, C_HYBRID] = apply_physical_constraints_correction(M_HYBRID, G_HYBRID, C_HYBRID, theory_info);

%% ====== 编译融合结果 ======
fusion_results = struct();
fusion_results.strategy = fusion_strategy;
fusion_results.weights = fusion_weights;
fusion_results.fusion_details = struct('M', fusion_details_M, 'G', fusion_details_G, 'C', fusion_details_C);
fusion_results.theory_info = theory_info;
fusion_results.id_results = id_results;
fusion_results.fusion_method = 'intelligent_adaptive';
fusion_results.fusion_timestamp = now;

fprintf('   智能混合辨识完成\n');

end

function fusion_strategy = analyze_identification_quality_and_theory_confidence(id_results, theory_info)
% 分析辨识质量与理论置信度

fusion_strategy = struct();

%% ====== 评估辨识数据质量 ======
if isfield(id_results, 'validation_results') && isfield(id_results.validation_results, 'mean_correlation')
    id_quality = id_results.validation_results.mean_correlation;
elseif isfield(id_results, 'analysis') && isfield(id_results.analysis, 'mean_correlation')
    id_quality = id_results.analysis.mean_correlation;
else
    id_quality = 0.5; % 默认中等质量
end

%% ====== 评估理论模型置信度 ======
if isfield(theory_info, 'confidence_weights') && isfield(theory_info.confidence_weights, 'overall')
    theory_confidence = theory_info.confidence_weights.overall;
else
    theory_confidence = 0.8; % MuJoCo模型默认高置信度
end

fusion_strategy.id_quality = id_quality;
fusion_strategy.theory_confidence = theory_confidence;

%% ====== 确定融合策略 ======
if id_quality > 0.7 && theory_confidence > 0.8
    fusion_strategy.strategy_type = 'balanced_high_quality';
    fusion_strategy.description = '高质量辨识+高置信度理论，平衡融合';
    
elseif id_quality > 0.6 && theory_confidence > 0.7
    fusion_strategy.strategy_type = 'balanced_medium_quality';
    fusion_strategy.description = '中高质量辨识+中高置信度理论，平衡融合';
    
elseif id_quality < 0.4 && theory_confidence > 0.8
    fusion_strategy.strategy_type = 'theory_dominant';
    fusion_strategy.description = '低质量辨识+高置信度理论，理论主导';
    
elseif id_quality > 0.7 && theory_confidence < 0.6
    fusion_strategy.strategy_type = 'identification_dominant';
    fusion_strategy.description = '高质量辨识+低置信度理论，辨识主导';
    
else
    fusion_strategy.strategy_type = 'conservative_balanced';
    fusion_strategy.description = '保守平衡融合';
end

fprintf('      融合策略: %s\n', fusion_strategy.description);
fprintf('      辨识质量: %.3f, 理论置信度: %.3f\n', id_quality, theory_confidence);

end

function fusion_weights = compute_adaptive_fusion_weights(M_ID, G_ID, M_THEORY, G_THEORY, ...
                                                         id_results, theory_info, fusion_strategy)
% 计算自适应融合权重

fusion_weights = struct();

%% ====== 基于策略的基础权重 ======
switch fusion_strategy.strategy_type
    case 'balanced_high_quality'
        base_id_weight = 0.6;
        
    case 'balanced_medium_quality'
        base_id_weight = 0.5;
        
    case 'theory_dominant'
        base_id_weight = 0.25;
        
    case 'identification_dominant'
        base_id_weight = 0.75;
        
    case 'conservative_balanced'
        base_id_weight = 0.45;
        
    otherwise
        base_id_weight = 0.5;
end

%% ====== DOF级别的自适应权重 ======
n_dof = min(size(M_ID, 1), size(M_THEORY, 1));
fusion_weights.dof_weights = zeros(n_dof, 1);

dof_names = {'X位置', 'Y位置', 'Z位置', 'Roll', 'Pitch', 'Yaw', 'Armz关节', 'Army关节', '绳索X', '绳索Y'};

for i = 1:n_dof
    % 获取理论置信度
    if isfield(theory_info, 'confidence_weights') && isfield(theory_info.confidence_weights, 'dof_confidence') ...
            && length(theory_info.confidence_weights.dof_confidence) >= i
        theory_conf_dof = theory_info.confidence_weights.dof_confidence(i);
    else
        theory_conf_dof = 0.7; % 默认值
    end
    
    % 获取辨识质量
    if isfield(id_results, 'quality_assessment') && isfield(id_results.quality_assessment, 'dof_quality') ...
            && isfield(id_results.quality_assessment.dof_quality, sprintf('dof_%d', i))
        id_quality_dof = id_results.quality_assessment.dof_quality.(sprintf('dof_%d', i)).overall_score;
    else
        id_quality_dof = 0.5; % 默认值
    end
    
    % 基于质量差异调整权重
    quality_ratio = id_quality_dof / max(theory_conf_dof, 0.1);
    
    if quality_ratio > 1.5
        % 辨识质量显著优于理论置信度
        dof_id_weight = min(base_id_weight + 0.2, 0.9);
    elseif quality_ratio < 0.7
        % 理论置信度显著优于辨识质量
        dof_id_weight = max(base_id_weight - 0.2, 0.1);
    else
        % 质量相当，使用基础权重
        dof_id_weight = base_id_weight;
    end
    
    % 特殊DOF调整
    if i == 3 % Z位置，重力影响大
        % 检查重力辨识精度
        mujoco_gravity = -theory_info.params.m_total * theory_info.params.g;
        gravity_error = abs(G_ID(3) - mujoco_gravity) / abs(mujoco_gravity);
        
        if gravity_error < 0.05
            dof_id_weight = min(dof_id_weight + 0.15, 0.85);
        elseif gravity_error > 0.2
            dof_id_weight = max(dof_id_weight - 0.2, 0.2);
        end
    end
    
    fusion_weights.dof_weights(i) = dof_id_weight;
    
    if i <= length(dof_names)
        fprintf('        %s: %.0f%%辨识 + %.0f%%理论 (质量比%.2f)\n', ...
                dof_names{i}, dof_id_weight*100, (1-dof_id_weight)*100, quality_ratio);
    end
end

%% ====== 特殊项权重 ======
% 重力项权重
mujoco_gravity = -theory_info.params.m_total * theory_info.params.g;
gravity_agreement = 1 - abs(G_ID(3) - mujoco_gravity) / abs(mujoco_gravity);
fusion_weights.gravity_weight = min(base_id_weight + 0.3 * gravity_agreement, 0.9);

% 耦合项权重（理论通常更可靠）
fusion_weights.coupling_weight = max(base_id_weight - 0.2, 0.2);

% 科氏项权重（理论更可靠）
fusion_weights.coriolis_weight = max(base_id_weight - 0.3, 0.15);

fprintf('        重力项: %.0f%%辨识 (一致性: %.1f%%)\n', fusion_weights.gravity_weight*100, gravity_agreement*100);
fprintf('        耦合项: %.0f%%辨识\n', fusion_weights.coupling_weight*100);
fprintf('        科氏项: %.0f%%辨识\n', fusion_weights.coriolis_weight*100);

end

function [M_HYBRID, fusion_details] = fuse_mass_matrix_intelligent(M_ID, M_THEORY, weights, strategy)
% 智能融合质量矩阵

M_HYBRID = zeros(size(M_THEORY));
fusion_details = struct();

fprintf('        质量矩阵智能融合...\n');

%% ====== 对角线项融合 ======
n_dof = size(M_HYBRID, 1);
diagonal_fusion_info = struct();

for i = 1:n_dof
    if i <= length(weights.dof_weights)
        w_id = weights.dof_weights(i);
    else
        w_id = 0.5;
    end
    w_theory = 1 - w_id;
    
    % 物理约束检查
    id_value = M_ID(i, i);
    theory_value = M_THEORY(i, i);
    
    % 异常值检测
    if id_value <= 0 || id_value > 10 * theory_value
        % 辨识值异常，更多依赖理论
        w_id = min(w_id, 0.2);
        w_theory = 1 - w_id;
        fprintf('          M(%d,%d): 辨识值异常，降低权重\n', i, i);
    end
    
    M_HYBRID(i, i) = w_id * id_value + w_theory * theory_value;
    
    diagonal_fusion_info.(sprintf('dof_%d', i)) = struct(...
        'w_id', w_id, 'w_theory', w_theory, ...
        'id_value', id_value, 'theory_value', theory_value, ...
        'fused_value', M_HYBRID(i, i));
    
    if i <= 8
        fprintf('          M(%d,%d): %.6f = %.2f×%.6f + %.2f×%.6f\n', ...
                i, i, M_HYBRID(i, i), w_id, id_value, w_theory, theory_value);
    end
end

%% ====== 耦合项融合 ======
w_coupling = weights.coupling_weight;
coupling_count = 0;

for i = 1:n_dof
    for j = 1:n_dof
        if i ~= j
            theory_coupling = M_THEORY(i, j);
            id_coupling = M_ID(i, j);
            
            % 如果理论有显著耦合项
            if abs(theory_coupling) > 1e-6
                M_HYBRID(i, j) = (1-w_coupling) * theory_coupling + w_coupling * id_coupling;
                coupling_count = coupling_count + 1;
            else
                % 理论无耦合，谨慎使用辨识结果
                M_HYBRID(i, j) = w_coupling * 0.3 * id_coupling;
            end
        end
    end
end

fprintf('          融合了 %d 个耦合项\n', coupling_count);

%% ====== 确保物理性质 ======
M_HYBRID = 0.5 * (M_HYBRID + M_HYBRID');
M_HYBRID = ensure_positive_definite_fusion(M_HYBRID);

%% ====== 融合质量评估 ======
[quality_score, quality_details] = assess_fusion_quality(M_HYBRID, M_ID, M_THEORY, weights);

fusion_details.diagonal_fusion = diagonal_fusion_info;
fusion_details.coupling_count = coupling_count;
fusion_details.quality_score = quality_score;
fusion_details.quality_details = quality_details;

fprintf('        质量矩阵融合质量: %.3f\n', quality_score);

end

function [G_HYBRID, fusion_details] = fuse_gravity_vector_intelligent(G_ID, G_THEORY, weights, strategy)
% 智能融合重力向量

G_HYBRID = zeros(size(G_THEORY));
fusion_details = struct();

fprintf('        重力向量智能融合...\n');

%% ====== 主重力项融合 ======
w_gravity = weights.gravity_weight;
G_HYBRID(3) = w_gravity * G_ID(3) + (1-w_gravity) * G_THEORY(3);

fprintf('          G_z: %.3f = %.2f×%.3f + %.2f×%.3f\n', ...
        G_HYBRID(3), w_gravity, G_ID(3), (1-w_gravity), G_THEORY(3));

%% ====== 其他项融合 ======
n_dof = length(G_THEORY);
for i = [1,2,4:n_dof]
    if i <= length(weights.dof_weights)
        w_local = weights.dof_weights(i) * 0.6; % 重力项的辨识通常不如主对角线项可靠
    else
        w_local = 0.3;
    end
    
    % 小重力项的融合
    G_HYBRID(i) = w_local * G_ID(i) + (1-w_local) * G_THEORY(i);
end

%% ====== 物理合理性检查 ======
% Z方向重力应该是负值且接近理论值
if G_HYBRID(3) > -5 || G_HYBRID(3) < -50
    fprintf('          Z重力值异常(%.3f)，使用理论值\n', G_HYBRID(3));
    G_HYBRID(3) = G_THEORY(3);
end

fusion_details.gravity_fusion = struct(...
    'w_gravity', w_gravity, ...
    'id_gravity', G_ID(3), 'theory_gravity', G_THEORY(3), ...
    'fused_gravity', G_HYBRID(3));

fprintf('        重力向量融合完成: G_z=%.3f N\n', G_HYBRID(3));

end

function [C_HYBRID, fusion_details] = fuse_coriolis_matrix_intelligent(C_ID, C_THEORY, weights, strategy)
% 智能融合科氏矩阵

w_coriolis = weights.coriolis_weight;

fprintf('        科氏矩阵智能融合...\n');

%% ====== 基本融合 ======
C_HYBRID = (1-w_coriolis) * C_THEORY + w_coriolis * C_ID;

%% ====== 确保反对称性 ======
C_HYBRID = 0.5 * (C_HYBRID - C_HYBRID');

%% ====== 验证反对称性 ======
symmetry_error = norm(C_HYBRID + C_HYBRID', 'fro');

fusion_details.coriolis_fusion = struct(...
    'w_coriolis', w_coriolis, ...
    'symmetry_error', symmetry_error, ...
    'norm_theory', norm(C_THEORY, 'fro'), ...
    'norm_identified', norm(C_ID, 'fro'), ...
    'norm_hybrid', norm(C_HYBRID, 'fro'));

if symmetry_error < 1e-10
    fprintf('         科氏矩阵反对称性验证通过\n');
else
    fprintf('        科氏矩阵反对称性误差: %.2e\n', symmetry_error);
end

fprintf('        科氏矩阵融合完成\n');

end

function [M_corrected, G_corrected, C_corrected] = apply_physical_constraints_correction(M_HYBRID, G_HYBRID, C_HYBRID, theory_info)
% 应用物理约束校正

fprintf('        应用物理约束校正...\n');

M_corrected = M_HYBRID;
G_corrected = G_HYBRID;
C_corrected = C_HYBRID;

%% ====== 质量矩阵约束 ======
% 总质量检查
total_mass_from_matrix = trace(M_corrected(1:3, 1:3)) / 3;
expected_mass = theory_info.params.m_total;
mass_error = abs(total_mass_from_matrix - expected_mass) / expected_mass;

if mass_error > 0.2
    fprintf('           总质量偏差过大(%.1f%%)，进行校正\n', mass_error*100);
    scale_factor = expected_mass / total_mass_from_matrix;
    M_corrected(1:3, 1:3) = M_corrected(1:3, 1:3) * scale_factor;
    fprintf('          校正后总质量: %.3f kg\n', expected_mass);
end

%% ====== 重力向量约束 ======
expected_gravity = -theory_info.params.m_total * theory_info.params.g;
gravity_error = abs(G_corrected(3) - expected_gravity) / abs(expected_gravity);

if gravity_error > 0.15
    fprintf('           重力项偏差过大(%.1f%%)，进行校正\n', gravity_error*100);
    G_corrected(3) = expected_gravity;
    fprintf('          校正后重力: %.3f N\n', expected_gravity);
end

%% ====== 正定性保证 ======
M_corrected = ensure_positive_definite_fusion(M_corrected);

%% ====== 条件数检查 ======
cond_num = cond(M_corrected);
if cond_num > 1e8
    fprintf('          条件数过大(%.2e)，应用正则化\n', cond_num);
    reg_factor = 1e-6;
    M_corrected = M_corrected + reg_factor * eye(size(M_corrected));
    fprintf('          正则化后条件数: %.2e\n', cond(M_corrected));
end

fprintf('        物理约束校正完成\n');

end

function M_pd = ensure_positive_definite_fusion(M)
% 确保正定性（融合专用）

M_pd = 0.5 * (M + M');

[V, D] = eig(M_pd);
eigenvals = diag(D);

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

M_pd = V * diag(eigenvals) * V';
M_pd = 0.5 * (M_pd + M_pd');

if corrected_count > 0
    fprintf('          修正了 %d 个特征值\n', corrected_count);
end

end

function [quality_score, quality_details] = assess_fusion_quality(M_HYBRID, M_ID, M_THEORY, weights)
% 评估融合质量

quality_details = struct();

%% ====== 对角线项质量 ======
n_dof = size(M_HYBRID, 1);
diagonal_scores = zeros(n_dof, 1);

for i = 1:n_dof
    hybrid_val = M_HYBRID(i, i);
    id_val = M_ID(i, i);
    theory_val = M_THEORY(i, i);
    
    % 计算混合值相对于输入值的合理性
    if i <= length(weights.dof_weights)
        w_id = weights.dof_weights(i);
        expected_val = w_id * id_val + (1-w_id) * theory_val;
        relative_error = abs(hybrid_val - expected_val) / max(abs(expected_val), 1e-6);
        diagonal_scores(i) = max(0, 1 - relative_error);
    else
        diagonal_scores(i) = 0.8; % 默认分数
    end
end

quality_details.diagonal_scores = diagonal_scores;
quality_details.mean_diagonal_score = mean(diagonal_scores);

%% ====== 物理性质检查 ======
eigenvals = eig(M_HYBRID);
quality_details.min_eigenvalue = min(eigenvals);
quality_details.condition_number = cond(M_HYBRID);
quality_details.is_positive_definite = all(eigenvals > 1e-10);

%% ====== 综合评分 ======
physical_score = 0;
if quality_details.is_positive_definite
    physical_score = physical_score + 0.4;
end
if quality_details.condition_number < 1e6
    physical_score = physical_score + 0.3;
end
if quality_details.min_eigenvalue > 1e-8
    physical_score = physical_score + 0.3;
end

quality_score = 0.7 * quality_details.mean_diagonal_score + 0.3 * physical_score;
quality_details.physical_score = physical_score;
quality_details.overall_score = quality_score;

end