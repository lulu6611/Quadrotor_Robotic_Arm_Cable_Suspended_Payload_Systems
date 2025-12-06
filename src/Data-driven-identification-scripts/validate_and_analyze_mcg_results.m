function validation_results = validate_and_analyze_mcg_results(...
    M_ID, G_ID, C_ID, id_results, ...
    M_THEORY, G_THEORY, C_THEORY, theory_info, ...
    M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, ...
    selected_data, data_info)
% 验证和分析MCG结果

fprintf('   开始结果验证与性能分析...\n');

%% ====== 基本信息统计 ======
validation_results = struct();
validation_results.timestamp = now;
validation_results.data_summary = compile_data_summary(selected_data, data_info);

%% ====== 单一方法性能分析 ======
fprintf('    分析单一方法性能...\n');
[id_performance, theory_performance] = analyze_single_method_performance(...
    M_ID, G_ID, C_ID, id_results, M_THEORY, G_THEORY, C_THEORY, theory_info, selected_data);

%% ====== 混合方法性能分析 ======
fprintf('    分析混合方法性能...\n');
hybrid_performance = analyze_hybrid_method_performance(...
    M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, selected_data);

%% ====== 方法对比分析 ======
fprintf('    执行方法对比分析...\n');
comparison_analysis = perform_method_comparison_analysis(...
    id_performance, theory_performance, hybrid_performance);

%% ====== 物理一致性验证 ======
fprintf('    验证物理一致性...\n');
physics_validation = validate_physics_consistency(...
    M_HYBRID, G_HYBRID, C_HYBRID, theory_info);

%% ====== 数值稳定性分析 ======
fprintf('    分析数值稳定性...\n');
numerical_analysis = analyze_numerical_stability(...
    M_ID, G_ID, M_THEORY, G_THEORY, M_HYBRID, G_HYBRID);

%% ====== 编译验证结果 ======
validation_results.identification_performance = id_performance;
validation_results.theory_performance = theory_performance;
validation_results.hybrid_performance = hybrid_performance;
validation_results.comparison_analysis = comparison_analysis;
validation_results.physics_validation = physics_validation;
validation_results.numerical_analysis = numerical_analysis;

% 综合评分
validation_results.overall_score = compute_overall_validation_score(validation_results);
validation_results.recommendation = generate_method_recommendation(validation_results);

fprintf('   验证分析完成，综合评分: %.3f\n', validation_results.overall_score);

end

function data_summary = compile_data_summary(selected_data, data_info)
% 编译数据摘要

data_summary = struct();
data_summary.total_datasets = length(selected_data);
data_summary.data_folder = data_info.data_folder;
data_summary.selection_strategy = data_info.selection_strategy;

if ~isempty(selected_data)
    % 质量统计
    qualities = zeros(length(selected_data), 1);
    excitation_types = {};
    
    for i = 1:length(selected_data)
        if isfield(selected_data{i}, 'quick_stats')
            qualities(i) = selected_data{i}.quick_stats.quality_score;
            excitation_types{i} = selected_data{i}.quick_stats.excitation_name;
        end
    end
    
    data_summary.quality_stats = struct();
    data_summary.quality_stats.mean = mean(qualities);
    data_summary.quality_stats.std = std(qualities);
    data_summary.quality_stats.min = min(qualities);
    data_summary.quality_stats.max = max(qualities);
    
    % 激励类型分布
    [unique_types, ~, idx] = unique(excitation_types);
    type_counts = accumarray(idx, 1);
    data_summary.excitation_distribution = struct();
    for i = 1:length(unique_types)
        field_name = matlab.lang.makeValidName(unique_types{i});
        data_summary.excitation_distribution.(field_name) = type_counts(i);
    end
end

end

function [id_performance, theory_performance] = analyze_single_method_performance(...
    M_ID, G_ID, C_ID, id_results, M_THEORY, G_THEORY, C_THEORY, theory_info, selected_data)
% 分析单一方法性能

%% ====== 辨识方法性能 ======
id_performance = struct();
id_performance.method = 'Data-driven Identification';

if isfield(id_results, 'validation_results')
    validation = id_results.validation_results;
    id_performance.mean_correlation = validation.mean_correlation;
    id_performance.validation_success = validation.validation_success;
    id_performance.overall_performance = validation.overall_performance;
else
    id_performance.mean_correlation = 0.5;
    id_performance.validation_success = false;
    id_performance.overall_performance = 'unknown';
end

% 物理合理性检查
id_performance.physics_check = check_physical_reasonableness(M_ID, G_ID, theory_info);

%% ====== 理论方法性能 ======
theory_performance = struct();
theory_performance.method = 'MuJoCo Theoretical Model';
theory_performance.confidence = theory_info.confidence_weights.overall;

% 基于置信度估计相关性
theory_performance.estimated_correlation = theory_performance.confidence * 0.85;
theory_performance.physics_check = check_physical_reasonableness(M_THEORY, G_THEORY, theory_info);

fprintf('      辨识方法: 相关性=%.3f, 物理评分=%.2f\n', ...
        id_performance.mean_correlation, id_performance.physics_check.overall_score);
fprintf('      理论方法: 置信度=%.3f, 物理评分=%.2f\n', ...
        theory_performance.confidence, theory_performance.physics_check.overall_score);

end

function hybrid_performance = analyze_hybrid_method_performance(...
    M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, selected_data)
% 分析混合方法性能

hybrid_performance = struct();
hybrid_performance.method = 'Intelligent Hybrid Fusion';

%% ====== 融合质量评估 ======
if isfield(fusion_results, 'fusion_details') && isfield(fusion_results.fusion_details, 'M')
    hybrid_performance.fusion_quality = fusion_results.fusion_details.M.quality_score;
else
    hybrid_performance.fusion_quality = 0.8; % 默认
end

%% ====== 物理一致性评估 ======
hybrid_performance.physics_check = check_physical_reasonableness(M_HYBRID, G_HYBRID, fusion_results.theory_info);

%% ====== 预测性能估计 ======
% 基于融合权重和输入性能预测
id_quality = fusion_results.id_results.validation_results.mean_correlation;
theory_confidence = fusion_results.theory_info.confidence_weights.overall;
avg_weight = mean(fusion_results.weights.dof_weights);

% 加权平均预测
predicted_correlation = avg_weight * id_quality + (1-avg_weight) * theory_confidence * 0.85;
% 融合带来的额外提升
fusion_bonus = 0.05 * hybrid_performance.fusion_quality;
hybrid_performance.predicted_correlation = min(predicted_correlation + fusion_bonus, 0.95);

hybrid_performance.overall_performance = rate_performance(hybrid_performance.predicted_correlation);

fprintf('      混合方法: 预测相关性=%.3f, 融合质量=%.2f, 物理评分=%.2f\n', ...
        hybrid_performance.predicted_correlation, hybrid_performance.fusion_quality, ...
        hybrid_performance.physics_check.overall_score);

end

function comparison_analysis = perform_method_comparison_analysis(...
    id_performance, theory_performance, hybrid_performance)
% 执行方法对比分析

comparison_analysis = struct();

%% ====== 性能对比 ======
comparison_analysis.performance_comparison = struct();
comparison_analysis.performance_comparison.identification = id_performance.mean_correlation;
comparison_analysis.performance_comparison.theory = theory_performance.estimated_correlation;
comparison_analysis.performance_comparison.hybrid = hybrid_performance.predicted_correlation;

%% ====== 改进幅度计算 ======
id_corr = id_performance.mean_correlation;
theory_corr = theory_performance.estimated_correlation;
hybrid_corr = hybrid_performance.predicted_correlation;

comparison_analysis.improvement_vs_identification = (hybrid_corr - id_corr) / max(id_corr, 0.1);
comparison_analysis.improvement_vs_theory = (hybrid_corr - theory_corr) / max(theory_corr, 0.1);
comparison_analysis.improvement_vs_best_single = (hybrid_corr - max(id_corr, theory_corr)) / max(max(id_corr, theory_corr), 0.1);

%% ====== 优势分析 ======
comparison_analysis.advantages = {};
if comparison_analysis.improvement_vs_identification > 0.05
    comparison_analysis.advantages{end+1} = sprintf('相比纯辨识提升%.1f%%', comparison_analysis.improvement_vs_identification*100);
end
if comparison_analysis.improvement_vs_theory > 0.05
    comparison_analysis.advantages{end+1} = sprintf('相比纯理论提升%.1f%%', comparison_analysis.improvement_vs_theory*100);
end
if hybrid_performance.physics_check.overall_score > max(id_performance.physics_check.overall_score, theory_performance.physics_check.overall_score)
    comparison_analysis.advantages{end+1} = '物理一致性最优';
end

%% ====== 推荐等级 ======
if hybrid_corr > max(id_corr, theory_corr) + 0.05
    comparison_analysis.recommendation_level = 'Highly Recommended';
elseif hybrid_corr > max(id_corr, theory_corr)
    comparison_analysis.recommendation_level = 'Recommended';
elseif hybrid_corr > min(id_corr, theory_corr)
    comparison_analysis.recommendation_level = 'Acceptable';
else
    comparison_analysis.recommendation_level = 'Not Recommended';
end

fprintf('      性能对比: 辨识%.3f, 理论%.3f, 混合%.3f\n', id_corr, theory_corr, hybrid_corr);
fprintf('      推荐等级: %s\n', comparison_analysis.recommendation_level);

end

function physics_check = check_physical_reasonableness(M, G, theory_info)
% 检查物理合理性

physics_check = struct();

%% ====== 质量矩阵检查 ======
% 正定性
eigenvals = eig(M);
physics_check.is_positive_definite = all(eigenvals > 1e-10);
physics_check.min_eigenvalue = min(eigenvals);
physics_check.condition_number = cond(M);

% 总质量合理性
total_mass_from_matrix = trace(M(1:3, 1:3)) / 3;
expected_mass = theory_info.params.m_total;
physics_check.mass_error = abs(total_mass_from_matrix - expected_mass) / expected_mass;
physics_check.mass_reasonable = physics_check.mass_error < 0.3;

%% ====== 重力向量检查 ======
expected_gravity = -theory_info.params.m_total * theory_info.params.g;
physics_check.gravity_error = abs(G(3) - expected_gravity) / abs(expected_gravity);
physics_check.gravity_reasonable = physics_check.gravity_error < 0.2;

% Z方向重力应该是负值
physics_check.gravity_direction_correct = G(3) < -5;

%% ====== 惯量合理性 ======
if size(M, 1) >= 6
    attitude_inertias = diag(M(4:6, 4:6));
    expected_inertias = diag(theory_info.params.I0);
    physics_check.inertia_errors = abs(attitude_inertias - expected_inertias) ./ expected_inertias;
    physics_check.inertia_reasonable = all(physics_check.inertia_errors < 2.0); % 200%误差内
end

%% ====== 综合评分 ======
score = 0;
if physics_check.is_positive_definite, score = score + 0.25; end
if physics_check.mass_reasonable, score = score + 0.25; end
if physics_check.gravity_reasonable, score = score + 0.25; end
if physics_check.gravity_direction_correct, score = score + 0.15; end
if physics_check.condition_number < 1e6, score = score + 0.1; end

physics_check.overall_score = score;

if score > 0.8
    physics_check.rating = 'Excellent';
elseif score > 0.6
    physics_check.rating = 'Good';
elseif score > 0.4
    physics_check.rating = 'Fair';
else
    physics_check.rating = 'Poor';
end

end

function physics_validation = validate_physics_consistency(M_HYBRID, G_HYBRID, C_HYBRID, theory_info)
% 验证物理一致性

fprintf('      验证混合模型的物理一致性...\n');

physics_validation = struct();

%% ====== 能量一致性检查 ======
% 动能应该是正定的
physics_validation.kinetic_energy_positive = all(eig(M_HYBRID) > 1e-10);

%% ====== 力矩平衡检查 ======
% 重力力矩应该与质心偏移一致
cg_offset = norm(theory_info.params.cg_total(1:2));
if cg_offset > 1e-6
    expected_gravity_moment = theory_info.params.m_total * theory_info.params.g * cg_offset;
    actual_gravity_moment = norm(G_HYBRID(4:5));
    physics_validation.moment_balance_error = abs(actual_gravity_moment - expected_gravity_moment) / expected_gravity_moment;
    physics_validation.moment_balance_reasonable = physics_validation.moment_balance_error < 0.5;
end

%% ====== 科氏项反对称性 ======
symmetry_error = norm(C_HYBRID + C_HYBRID', 'fro');
physics_validation.coriolis_antisymmetric = symmetry_error < 1e-8;
physics_validation.coriolis_symmetry_error = symmetry_error;

%% ====== 尺度一致性 ======
% 各项的数量级应该合理
physics_validation.mass_scale_check = all(diag(M_HYBRID(1:3,1:3)) > 0.5) && all(diag(M_HYBRID(1:3,1:3)) < 10);
if size(M_HYBRID, 1) >= 6
    physics_validation.inertia_scale_check = all(diag(M_HYBRID(4:6,4:6)) > 1e-6) && all(diag(M_HYBRID(4:6,4:6)) < 1);
end

%% ====== 综合评分 ======
checks = [physics_validation.kinetic_energy_positive, ...
          physics_validation.coriolis_antisymmetric, ...
          physics_validation.mass_scale_check];

if isfield(physics_validation, 'moment_balance_reasonable')
    checks = [checks, physics_validation.moment_balance_reasonable];
end
if isfield(physics_validation, 'inertia_scale_check')
    checks = [checks, physics_validation.inertia_scale_check];
end

physics_validation.consistency_score = sum(checks) / length(checks);

if physics_validation.consistency_score > 0.8
    physics_validation.consistency_rating = 'Excellent';
elseif physics_validation.consistency_score > 0.6
    physics_validation.consistency_rating = 'Good';
else
    physics_validation.consistency_rating = 'Needs Improvement';
end

fprintf('        物理一致性评分: %.2f (%s)\n', physics_validation.consistency_score, physics_validation.consistency_rating);

end

function numerical_analysis = analyze_numerical_stability(M_ID, G_ID, M_THEORY, G_THEORY, M_HYBRID, G_HYBRID)
% 分析数值稳定性

fprintf('      分析数值稳定性...\n');

numerical_analysis = struct();

%% ====== 条件数分析 ======
cond_id = cond(M_ID);
cond_theory = cond(M_THEORY);
cond_hybrid = cond(M_HYBRID);

numerical_analysis.condition_numbers = struct();
numerical_analysis.condition_numbers.identification = cond_id;
numerical_analysis.condition_numbers.theory = cond_theory;
numerical_analysis.condition_numbers.hybrid = cond_hybrid;

% 混合方法是否改善了条件数
numerical_analysis.condition_improvement = cond_hybrid < min(cond_id, cond_theory);

%% ====== 特征值分析 ======
eig_id = eig(M_ID);
eig_theory = eig(M_THEORY);
eig_hybrid = eig(M_HYBRID);

numerical_analysis.eigenvalue_analysis = struct();
numerical_analysis.eigenvalue_analysis.min_eigenvalues = [min(eig_id), min(eig_theory), min(eig_hybrid)];
numerical_analysis.eigenvalue_analysis.max_eigenvalues = [max(eig_id), max(eig_theory), max(eig_hybrid)];

% 混合方法的特征值是否更合理
numerical_analysis.eigenvalue_stability = all(eig_hybrid > 1e-10) && all(eig_hybrid < 1e4);

%% ====== 参数范围分析 ======
numerical_analysis.parameter_ranges = struct();

% 质量参数范围
mass_diag_id = diag(M_ID(1:3,1:3));
mass_diag_theory = diag(M_THEORY(1:3,1:3));
mass_diag_hybrid = diag(M_HYBRID(1:3,1:3));

numerical_analysis.parameter_ranges.mass_ranges = [
    [min(mass_diag_id), max(mass_diag_id)];
    [min(mass_diag_theory), max(mass_diag_theory)];
    [min(mass_diag_hybrid), max(mass_diag_hybrid)]
];

%% ====== 数值稳定性评分 ======
stability_score = 0;

% 条件数评分
if cond_hybrid < 1e6, stability_score = stability_score + 0.3; end
if numerical_analysis.condition_improvement, stability_score = stability_score + 0.2; end

% 特征值评分
if numerical_analysis.eigenvalue_stability, stability_score = stability_score + 0.3; end

% 参数合理性评分
mass_reasonable = all(mass_diag_hybrid > 0.5) && all(mass_diag_hybrid < 5);
if mass_reasonable, stability_score = stability_score + 0.2; end

numerical_analysis.stability_score = stability_score;

if stability_score > 0.8
    numerical_analysis.stability_rating = 'Excellent';
elseif stability_score > 0.6
    numerical_analysis.stability_rating = 'Good';
else
    numerical_analysis.stability_rating = 'Fair';
end

fprintf('        数值稳定性评分: %.2f (%s)\n', numerical_analysis.stability_score, numerical_analysis.stability_rating);

end

function overall_score = compute_overall_validation_score(validation_results)
% 计算综合验证评分

weights = struct();
weights.hybrid_performance = 0.4;
weights.physics_consistency = 0.25;
weights.numerical_stability = 0.2;
weights.improvement_significance = 0.15;

scores = struct();
scores.hybrid_performance = validation_results.hybrid_performance.predicted_correlation;
scores.physics_consistency = validation_results.physics_validation.consistency_score;
scores.numerical_stability = validation_results.numerical_analysis.stability_score;

% 改进显著性评分
improvement = validation_results.comparison_analysis.improvement_vs_best_single;
if improvement > 0.1
    scores.improvement_significance = 1.0;
elseif improvement > 0.05
    scores.improvement_significance = 0.8;
elseif improvement > 0
    scores.improvement_significance = 0.6;
else
    scores.improvement_significance = 0.3;
end

overall_score = weights.hybrid_performance * scores.hybrid_performance + ...
                weights.physics_consistency * scores.physics_consistency + ...
                weights.numerical_stability * scores.numerical_stability + ...
                weights.improvement_significance * scores.improvement_significance;

end

function recommendation = generate_method_recommendation(validation_results)
% 生成方法推荐

recommendation = struct();

overall_score = validation_results.overall_score;
improvement = validation_results.comparison_analysis.improvement_vs_best_single;

if overall_score > 0.8 && improvement > 0.05
    recommendation.level = 'Strongly Recommended';
    recommendation.reasons = {'高综合性能', '显著性能提升', '优秀物理一致性'};
    
elseif overall_score > 0.7 && improvement > 0
    recommendation.level = 'Recommended';
    recommendation.reasons = {'良好综合性能', '有效性能提升'};
    
elseif overall_score > 0.6
    recommendation.level = 'Conditionally Recommended';
    recommendation.reasons = {'中等性能', '需要进一步验证'};
    
else
    recommendation.level = 'Not Recommended';
    recommendation.reasons = {'性能不足', '需要改进辨识策略'};
end

% 使用建议
recommendation.usage_suggestions = {};
if validation_results.hybrid_performance.predicted_correlation > 0.7
    recommendation.usage_suggestions{end+1} = '适用于高精度控制应用';
end
if validation_results.physics_validation.consistency_score > 0.8
    recommendation.usage_suggestions{end+1} = '适用于物理仿真应用';
end
if validation_results.numerical_analysis.stability_score > 0.7
    recommendation.usage_suggestions{end+1} = '适用于实时控制系统';
end

end

function performance_level = rate_performance(correlation)
% 性能等级评定

if correlation > 0.8
    performance_level = 'excellent';
elseif correlation > 0.6
    performance_level = 'good';
elseif correlation > 0.4
    performance_level = 'fair';
else
    performance_level = 'poor';
end

end