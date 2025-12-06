function save_mcg_processing_results(M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results, ...
                                   M_THEORY, G_THEORY, C_THEORY, theory_info, ...
                                   M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, ...
                                   validation_results, selected_data, data_info)
% 保存MCG处理结果

fprintf('  保存MCG处理结果...\n');

%% ====== 创建结果文件夹 ======
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_folder = sprintf('MCG_Processing_Results_%s', timestamp);

if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

fprintf('     结果文件夹: %s\n', results_folder);

%% ====== 保存主要结果 ======
main_results_file = fullfile(results_folder, 'mcg_main_results.mat');

try
    save(main_results_file, ...
         'M_IDENTIFIED', 'G_IDENTIFIED', 'C_IDENTIFIED', 'id_results', ...
         'M_THEORY', 'G_THEORY', 'C_THEORY', 'theory_info', ...
         'M_HYBRID', 'G_HYBRID', 'C_HYBRID', 'fusion_results', ...
         'validation_results', 'selected_data', 'data_info', '-v7.3');
    
    fprintf('    主要结果已保存: %s (%.1f MB)\n', main_results_file, ...
            get_file_size_mb(main_results_file));
catch ME
    fprintf('    主要结果保存失败: %s\n', ME.message);
end

%% ====== 保存到工作空间 ======
try
    assignin('base', 'M_IDENTIFIED_FINAL', M_IDENTIFIED);
    assignin('base', 'G_IDENTIFIED_FINAL', G_IDENTIFIED);
    assignin('base', 'C_IDENTIFIED_FINAL', C_IDENTIFIED);
    
    assignin('base', 'M_THEORY_FINAL', M_THEORY);
    assignin('base', 'G_THEORY_FINAL', G_THEORY);
    assignin('base', 'C_THEORY_FINAL', C_THEORY);
    
    assignin('base', 'M_HYBRID_FINAL', M_HYBRID);
    assignin('base', 'G_HYBRID_FINAL', G_HYBRID);
    assignin('base', 'C_HYBRID_FINAL', C_HYBRID);
    
    assignin('base', 'mcg_validation_results', validation_results);
    assignin('base', 'mcg_fusion_results', fusion_results);
    
    fprintf('    结果已保存到工作空间\n');
catch ME
    fprintf('    工作空间保存部分失败: %s\n', ME.message);
end

%% ====== 生成详细报告 ======
fprintf('     生成详细分析报告...\n');
generate_comprehensive_analysis_report(results_folder, validation_results, fusion_results, data_info);

%% ====== 生成MCG矩阵文件 ======
fprintf('    导出MCG矩阵...\n');
export_mcg_matrices(results_folder, M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, ...
                    M_THEORY, G_THEORY, C_THEORY, M_HYBRID, G_HYBRID, C_HYBRID);

%% ====== 生成可视化图表 ======
fprintf('     生成性能对比图表...\n');
generate_performance_visualization(results_folder, validation_results);

%% ====== 生成使用说明 ======
fprintf('     生成使用说明...\n');
generate_usage_instructions(results_folder, validation_results);

fprintf('   结果保存完成\n');

end

function generate_comprehensive_analysis_report(results_folder, validation_results, fusion_results, data_info)
% 生成综合分析报告

report_file = fullfile(results_folder, 'MCG_Comprehensive_Analysis_Report.txt');

fid = fopen(report_file, 'w');
if fid == -1
    fprintf('     无法创建分析报告文件\n');
    return;
end

try
    %% ====== 报告头部 ======
    fprintf(fid, '==================== MCG系统辨识综合分析报告 ====================\n');
    fprintf(fid, '生成时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '======================================================================\n\n');
    
    %% ====== 数据概况 ======
    fprintf(fid, ' 数据概况\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    fprintf(fid, '数据源文件夹: %s\n', data_info.data_folder);
    fprintf(fid, '数据集总数: %d\n', validation_results.data_summary.total_datasets);
    fprintf(fid, '选择策略: %s\n', data_info.selection_strategy);
    
    if isfield(validation_results.data_summary, 'quality_stats')
        stats = validation_results.data_summary.quality_stats;
        fprintf(fid, '数据质量统计:\n');
        fprintf(fid, '  平均质量: %.3f ± %.3f\n', stats.mean, stats.std);
        fprintf(fid, '  质量范围: [%.3f, %.3f]\n', stats.min, stats.max);
    end
    
    if isfield(validation_results.data_summary, 'excitation_distribution')
        fprintf(fid, '激励类型分布:\n');
        exc_fields = fieldnames(validation_results.data_summary.excitation_distribution);
        for i = 1:length(exc_fields)
            count = validation_results.data_summary.excitation_distribution.(exc_fields{i});
            fprintf(fid, '  %s: %d 个数据集\n', exc_fields{i}, count);
        end
    end
    fprintf(fid, '\n');
    
    %% ====== 方法性能对比 ======
    fprintf(fid, ' 方法性能对比\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    
    id_perf = validation_results.identification_performance;
    theory_perf = validation_results.theory_performance;
    hybrid_perf = validation_results.hybrid_performance;
    
    fprintf(fid, '1. 数据驱动辨识方法:\n');
    fprintf(fid, '   相关系数: %.3f\n', id_perf.mean_correlation);
    fprintf(fid, '   验证状态: %s\n', bool_to_string(id_perf.validation_success));
    fprintf(fid, '   整体性能: %s\n', id_perf.overall_performance);
    fprintf(fid, '   物理评分: %.2f (%s)\n', id_perf.physics_check.overall_score, id_perf.physics_check.rating);
    fprintf(fid, '\n');
    
    fprintf(fid, '2. MuJoCo理论模型:\n');
    fprintf(fid, '   理论置信度: %.3f\n', theory_perf.confidence);
    fprintf(fid, '   估计相关性: %.3f\n', theory_perf.estimated_correlation);
    fprintf(fid, '   物理评分: %.2f (%s)\n', theory_perf.physics_check.overall_score, theory_perf.physics_check.rating);
    fprintf(fid, '\n');
    
    fprintf(fid, '3. 智能混合方法:\n');
    fprintf(fid, '   预测相关性: %.3f\n', hybrid_perf.predicted_correlation);
    fprintf(fid, '   融合质量: %.2f\n', hybrid_perf.fusion_quality);
    fprintf(fid, '   整体性能: %s\n', hybrid_perf.overall_performance);
    fprintf(fid, '   物理评分: %.2f (%s)\n', hybrid_perf.physics_check.overall_score, hybrid_perf.physics_check.rating);
    fprintf(fid, '\n');
    
    %% ====== 性能改进分析 ======
    fprintf(fid, ' 性能改进分析\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    comp = validation_results.comparison_analysis;
    
    fprintf(fid, '相对于纯辨识方法: %+.1f%%\n', comp.improvement_vs_identification * 100);
    fprintf(fid, '相对于纯理论方法: %+.1f%%\n', comp.improvement_vs_theory * 100);
    fprintf(fid, '相对于最佳单一方法: %+.1f%%\n', comp.improvement_vs_best_single * 100);
    fprintf(fid, '\n推荐等级: %s\n', comp.recommendation_level);
    
    if ~isempty(comp.advantages)
        fprintf(fid, '主要优势:\n');
        for i = 1:length(comp.advantages)
            fprintf(fid, '  • %s\n', comp.advantages{i});
        end
    end
    fprintf(fid, '\n');
    
    %% ====== 融合策略详情 ======
    fprintf(fid, ' 融合策略详情\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    fprintf(fid, '融合策略: %s\n', fusion_results.strategy.description);
    fprintf(fid, '辨识质量: %.3f\n', fusion_results.strategy.id_quality);
    fprintf(fid, '理论置信度: %.3f\n', fusion_results.strategy.theory_confidence);
    fprintf(fid, '\n');
    
    fprintf(fid, 'DOF级别权重分配:\n');
    dof_names = {'X位置', 'Y位置', 'Z位置', 'Roll', 'Pitch', 'Yaw', 'Armz关节', 'Army关节'};
    for i = 1:min(length(fusion_results.weights.dof_weights), length(dof_names))
        w_id = fusion_results.weights.dof_weights(i);
        fprintf(fid, '  %s: %.0f%%辨识 + %.0f%%理论\n', dof_names{i}, w_id*100, (1-w_id)*100);
    end
    fprintf(fid, '\n');
    
    %% ====== 物理一致性验证 ======
    fprintf(fid, ' 物理一致性验证\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    phys = validation_results.physics_validation;
    
    fprintf(fid, '动能正定性: %s\n', bool_to_string(phys.kinetic_energy_positive));
    fprintf(fid, '科氏矩阵反对称性: %s (误差: %.2e)\n', bool_to_string(phys.coriolis_antisymmetric), phys.coriolis_symmetry_error);
    fprintf(fid, '质量尺度合理性: %s\n', bool_to_string(phys.mass_scale_check));
    
    if isfield(phys, 'inertia_scale_check')
        fprintf(fid, '惯量尺度合理性: %s\n', bool_to_string(phys.inertia_scale_check));
    end
    if isfield(phys, 'moment_balance_reasonable')
        fprintf(fid, '力矩平衡合理性: %s (误差: %.1f%%)\n', bool_to_string(phys.moment_balance_reasonable), phys.moment_balance_error*100);
    end
    
    fprintf(fid, '\n物理一致性评分: %.2f (%s)\n\n', phys.consistency_score, phys.consistency_rating);
    
    %% ====== 数值稳定性分析 ======
    fprintf(fid, '数值稳定性分析\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    num = validation_results.numerical_analysis;
    
    fprintf(fid, '条件数:\n');
    fprintf(fid, '  辨识方法: %.2e\n', num.condition_numbers.identification);
    fprintf(fid, '  理论方法: %.2e\n', num.condition_numbers.theory);
    fprintf(fid, '  混合方法: %.2e\n', num.condition_numbers.hybrid);
    fprintf(fid, '  条件数改善: %s\n', bool_to_string(num.condition_improvement));
    fprintf(fid, '\n');
    
    fprintf(fid, '特征值稳定性: %s\n', bool_to_string(num.eigenvalue_stability));
    fprintf(fid, '数值稳定性评分: %.2f (%s)\n\n', num.stability_score, num.stability_rating);
    
    %% ====== 综合评价 ======
    fprintf(fid, ' 综合评价\n');
    fprintf(fid, '────────────────────────────────────────────────────────────────────\n');
    fprintf(fid, '综合评分: %.3f / 1.000\n', validation_results.overall_score);
    fprintf(fid, '推荐等级: %s\n', validation_results.recommendation.level);
    fprintf(fid, '\n推荐理由:\n');
    
    for i = 1:length(validation_results.recommendation.reasons)
        fprintf(fid, '  • %s\n', validation_results.recommendation.reasons{i});
    end
    
    if isfield(validation_results.recommendation, 'usage_suggestions') && ~isempty(validation_results.recommendation.usage_suggestions)
        fprintf(fid, '\n使用建议:\n');
        for i = 1:length(validation_results.recommendation.usage_suggestions)
            fprintf(fid, '  • %s\n', validation_results.recommendation.usage_suggestions{i});
        end
    end
    fprintf(fid, '\n');
    
    %% ====== 报告尾部 ======
    fprintf(fid, '======================================================================\n');
    fprintf(fid, '报告生成完成 - MCG系统辨识分析\n');
    fprintf(fid, '======================================================================\n');
    
    fclose(fid);
    fprintf('       综合分析报告: %s\n', report_file);
    
catch ME
    fclose(fid);
    fprintf('       报告生成失败: %s\n', ME.message);
end

end

function export_mcg_matrices(results_folder, M_ID, G_ID, C_ID, M_THEORY, G_THEORY, C_THEORY, M_HYBRID, G_HYBRID, C_HYBRID)
% 导出MCG矩阵

%% ====== 创建矩阵子文件夹 ======
matrices_folder = fullfile(results_folder, 'MCG_Matrices');
if ~exist(matrices_folder, 'dir')
    mkdir(matrices_folder);
end

%% ====== 导出为CSV文件 ======
try
    % 辨识结果
    writematrix(M_ID, fullfile(matrices_folder, 'M_Identified.csv'));
    writematrix(G_ID, fullfile(matrices_folder, 'G_Identified.csv'));
    writematrix(C_ID, fullfile(matrices_folder, 'C_Identified.csv'));
    
    % 理论结果
    writematrix(M_THEORY, fullfile(matrices_folder, 'M_Theory.csv'));
    writematrix(G_THEORY, fullfile(matrices_folder, 'G_Theory.csv'));
    writematrix(C_THEORY, fullfile(matrices_folder, 'C_Theory.csv'));
    
    % 混合结果
    writematrix(M_HYBRID, fullfile(matrices_folder, 'M_Hybrid.csv'));
    writematrix(G_HYBRID, fullfile(matrices_folder, 'G_Hybrid.csv'));
    writematrix(C_HYBRID, fullfile(matrices_folder, 'C_Hybrid.csv'));
    
    fprintf('       MCG矩阵已导出为CSV格式\n');
    
catch ME
    fprintf('      ️ CSV导出部分失败: %s\n', ME.message);
end

%% ====== 导出为MAT文件 ======
try
    matrices_mat_file = fullfile(matrices_folder, 'mcg_matrices.mat');
    save(matrices_mat_file, 'M_ID', 'G_ID', 'C_ID', 'M_THEORY', 'G_THEORY', 'C_THEORY', ...
         'M_HYBRID', 'G_HYBRID', 'C_HYBRID');
    
    fprintf('       MCG矩阵已保存为MAT文件\n');
    
catch ME
    fprintf('       MAT文件保存失败: %s\n', ME.message);
end

%% ====== 生成矩阵说明文件 ======
readme_file = fullfile(matrices_folder, 'README_Matrices.txt');
fid = fopen(readme_file, 'w');

if fid ~= -1
    fprintf(fid, 'MCG矩阵文件说明\n');
    fprintf(fid, '================\n\n');
    fprintf(fid, '本文件夹包含三种方法得到的MCG矩阵：\n\n');
    fprintf(fid, '1. 数据驱动辨识结果：\n');
    fprintf(fid, '   - M_Identified.csv/mat: 质量/惯量矩阵\n');
    fprintf(fid, '   - G_Identified.csv/mat: 重力向量\n');
    fprintf(fid, '   - C_Identified.csv/mat: 科氏矩阵\n\n');
    fprintf(fid, '2. MuJoCo理论模型：\n');
    fprintf(fid, '   - M_Theory.csv/mat: 理论质量/惯量矩阵\n');
    fprintf(fid, '   - G_Theory.csv/mat: 理论重力向量\n');
    fprintf(fid, '   - C_Theory.csv/mat: 理论科氏矩阵\n\n');
    fprintf(fid, '3. 智能混合结果（推荐使用）：\n');
    fprintf(fid, '   - M_Hybrid.csv/mat: 混合质量/惯量矩阵\n');
    fprintf(fid, '   - G_Hybrid.csv/mat: 混合重力向量\n');
    fprintf(fid, '   - C_Hybrid.csv/mat: 混合科氏矩阵\n\n');
    fprintf(fid, 'DOF顺序：\n');
    fprintf(fid, '1-3: 位置 (X前, Y左, Z上)\n');
    fprintf(fid, '4-6: 姿态 (Roll, Pitch, Yaw)\n');
    fprintf(fid, '7-8: 关节 (Armz, Army)\n');
    fprintf(fid, '9-10: 绳索摆角 (θx, θy)\n\n');
    fprintf(fid, '坐标系：FLU (Forward-Left-Up)\n');
    
    fclose(fid);
    fprintf('       矩阵说明文件已生成\n');
end

end

function generate_performance_visualization(results_folder, validation_results)
% 生成性能对比可视化

try
    %% ====== 方法性能对比雷达图 ======
    figure('Position', [100, 100, 1200, 800]);
    
    % 准备数据
    methods = {'数据驱动辨识', 'MuJoCo理论', '智能混合'};
    
    id_perf = validation_results.identification_performance;
    theory_perf = validation_results.theory_performance;
    hybrid_perf = validation_results.hybrid_performance;
    
    % 性能指标
    correlations = [id_perf.mean_correlation, theory_perf.estimated_correlation, hybrid_perf.predicted_correlation];
    physics_scores = [id_perf.physics_check.overall_score, theory_perf.physics_check.overall_score, hybrid_perf.physics_check.overall_score];
    
    % 子图1：相关性对比
    subplot(2, 2, 1);
    bar(correlations, 'FaceColor', [0.2, 0.6, 0.8]);
    set(gca, 'XTickLabel', methods);
    ylabel('相关系数');
    title('预测性能对比');
    ylim([0, 1]);
    grid on;
    
    % 添加数值标签
    for i = 1:length(correlations)
        text(i, correlations(i) + 0.02, sprintf('%.3f', correlations(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % 子图2：物理一致性对比
    subplot(2, 2, 2);
    bar(physics_scores, 'FaceColor', [0.8, 0.4, 0.2]);
    set(gca, 'XTickLabel', methods);
    ylabel('物理一致性评分');
    title('物理一致性对比');
    ylim([0, 1]);
    grid on;
    
    for i = 1:length(physics_scores)
        text(i, physics_scores(i) + 0.02, sprintf('%.2f', physics_scores(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % 子图3：改进幅度分析
    subplot(2, 2, 3);
    comp = validation_results.comparison_analysis;
    improvements = [comp.improvement_vs_identification, comp.improvement_vs_theory, comp.improvement_vs_best_single] * 100;
    improvement_labels = {'vs 辨识', 'vs 理论', 'vs 最佳单一'};
    
    colors = [0.2, 0.8, 0.2; 0.8, 0.6, 0.2; 0.6, 0.2, 0.8];
    bar_handle = bar(improvements);
    bar_handle.FaceColor = 'flat';
    bar_handle.CData = colors;
    
    set(gca, 'XTickLabel', improvement_labels);
    ylabel('性能提升 (%)');
    title('混合方法性能提升');
    grid on;
    
    for i = 1:length(improvements)
        if improvements(i) >= 0
            text(i, improvements(i) + 0.5, sprintf('+%.1f%%', improvements(i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        else
            text(i, improvements(i) - 0.5, sprintf('%.1f%%', improvements(i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    
    % 子图4：综合评价
    subplot(2, 2, 4);
    
    % 综合评分的各个组成部分
    overall_score = validation_results.overall_score;
    hybrid_contribution = hybrid_perf.predicted_correlation * 0.4;
    physics_contribution = validation_results.physics_validation.consistency_score * 0.25;
    numerical_contribution = validation_results.numerical_analysis.stability_score * 0.2;
    improvement_contribution = min(comp.improvement_vs_best_single, 0.1) * 1.5; % 改进显著性贡献
    
    contributions = [hybrid_contribution, physics_contribution, numerical_contribution, improvement_contribution];
    contribution_labels = {'性能(40%)', '物理(25%)', '数值(20%)', '改进(15%)'};
    
    pie_colors = [0.3, 0.7, 0.9; 0.9, 0.6, 0.3; 0.6, 0.9, 0.4; 0.9, 0.4, 0.6];
    pie(contributions, contribution_labels);
    colormap(pie_colors);
    title(sprintf('综合评分构成 (总分: %.3f)', overall_score));
    
    % 保存图像
    sgtitle('MCG系统辨识方法性能对比分析', 'FontSize', 16, 'FontWeight', 'bold');
    
    performance_fig_file = fullfile(results_folder, 'Performance_Comparison.png');
    saveas(gcf, performance_fig_file);
    
    performance_fig_file_fig = fullfile(results_folder, 'Performance_Comparison.fig');
    saveas(gcf, performance_fig_file_fig);
    
    close(gcf);
    
    fprintf('       性能对比图表已生成\n');
    
catch ME
    fprintf('       图表生成失败: %s\n', ME.message);
    if exist('gcf', 'var')
        close(gcf);
    end
end

end

function generate_usage_instructions(results_folder, validation_results)
% 生成使用说明

instructions_file = fullfile(results_folder, 'Usage_Instructions.txt');
fid = fopen(instructions_file, 'w');

if fid == -1
    fprintf('       无法创建使用说明文件\n');
    return;
end

try
    fprintf(fid, 'MCG系统辨识结果使用说明\n');
    fprintf(fid, '==========================\n\n');
    
    fprintf(fid, ' 结果文件结构\n');
    fprintf(fid, '────────────────────────────────────────\n');
    fprintf(fid, '本文件夹包含以下内容：\n');
    fprintf(fid, '• mcg_main_results.mat - 完整结果数据\n');
    fprintf(fid, '• MCG_Matrices/ - MCG矩阵文件\n');
    fprintf(fid, '• Performance_Comparison.png/.fig - 性能对比图表\n');
    fprintf(fid, '• MCG_Comprehensive_Analysis_Report.txt - 详细分析报告\n');
    fprintf(fid, '• Usage_Instructions.txt - 本使用说明\n\n');
    
    fprintf(fid, ' 推荐使用方法\n');
    fprintf(fid, '────────────────────────────────────────\n');
    
    % 根据验证结果给出具体建议
    overall_score = validation_results.overall_score;
    recommendation = validation_results.recommendation;
    
    if overall_score > 0.8
        fprintf(fid, ' 推荐等级：%s\n', recommendation.level);
        fprintf(fid, '推荐优先使用混合方法结果（M_HYBRID, G_HYBRID, C_HYBRID）\n\n');
        
        fprintf(fid, '原因：\n');
        for i = 1:length(recommendation.reasons)
            fprintf(fid, '• %s\n', recommendation.reasons{i});
        end
        
    elseif overall_score > 0.6
        fprintf(fid, ' 推荐等级：%s\n', recommendation.level);
        fprintf(fid, '可以使用混合方法结果，但建议进行进一步验证\n\n');
        
    else
        fprintf(fid, ' 推荐等级：%s\n', recommendation.level);
        fprintf(fid, '不推荐直接使用，建议改进数据质量后重新辨识\n\n');
    end
    
    fprintf(fid, ' 在MATLAB中加载和使用\n');
    fprintf(fid, '────────────────────────────────────────\n');
    fprintf(fid, '%% 加载完整结果\n');
    fprintf(fid, 'load(''mcg_main_results.mat'');\n\n');
    
    fprintf(fid, '%% 使用混合MCG模型进行动力学计算\n');
    fprintf(fid, '%% 假设当前状态为 q, q_dot, q_ddot, 控制输入为 u\n');
    fprintf(fid, 'tau_dynamics = M_HYBRID * q_ddot + C_HYBRID * q_dot + G_HYBRID;\n');
    fprintf(fid, 'residual = u - tau_dynamics; %% 控制误差\n\n');
    
    fprintf(fid, '%% 或者直接从工作空间使用\n');
    fprintf(fid, 'M = M_HYBRID_FINAL;\n');
    fprintf(fid, 'G = G_HYBRID_FINAL;\n');
    fprintf(fid, 'C = C_HYBRID_FINAL;\n\n');
    
    fprintf(fid, '🔧 控制器设计建议\n');
    fprintf(fid, '────────────────────────────────────────\n');
    
    if validation_results.numerical_analysis.stability_score > 0.7
        fprintf(fid, ' 数值稳定性良好，可直接用于控制器设计\n');
        fprintf(fid, '建议的控制律形式：\n');
        fprintf(fid, 'u = M_HYBRID * (q_ddot_d + Kd*(q_dot_d - q_dot) + Kp*(q_d - q)) + C_HYBRID * q_dot + G_HYBRID\n\n');
    else
        fprintf(fid, '数值稳定性一般，建议添加正则化项\n');
        fprintf(fid, 'M_reg = M_HYBRID + lambda * eye(size(M_HYBRID)); %% lambda = 1e-3 到 1e-6\n\n');
    end
    
    fprintf(fid, ' 坐标系说明\n');
    fprintf(fid, '────────────────────────────────────────\n');
    fprintf(fid, '本结果基于FLU坐标系：\n');
    fprintf(fid, '• X轴：前向 (Forward)\n');
    fprintf(fid, '• Y轴：左向 (Left)\n');
    fprintf(fid, '• Z轴：上向 (Up)\n\n');
    
    fprintf(fid, 'DOF编号对应：\n');
    fprintf(fid, '1-3: 位置 [X, Y, Z]\n');
    fprintf(fid, '4-6: 姿态 [Roll, Pitch, Yaw]\n');
    fprintf(fid, '7-8: 关节角 [Armz, Army]\n');
    fprintf(fid, '9-10: 绳索摆角 [θx, θy]\n\n');
    
    fprintf(fid, ' 注意事项\n');
    fprintf(fid, '────────────────────────────────────────\n');
    fprintf(fid, '• 重力项G_HYBRID(3)为负值（向上为正方向）\n');
    fprintf(fid, '• 质量矩阵M_HYBRID必须保持正定性\n');
    fprintf(fid, '• 科氏矩阵C_HYBRID具有反对称性质\n');
    fprintf(fid, '• 使用前请检查矩阵的条件数和特征值\n\n');
    
    if isfield(validation_results.recommendation, 'usage_suggestions')
        fprintf(fid, ' 应用建议\n');
        fprintf(fid, '────────────────────────────────────────\n');
        for i = 1:length(validation_results.recommendation.usage_suggestions)
            fprintf(fid, '• %s\n', validation_results.recommendation.usage_suggestions{i});
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '技术支持\n');
    fprintf(fid, '────────────────────────────────────────\n');
    fprintf(fid, '如有问题请参考：\n');
    fprintf(fid, '• MCG_Comprehensive_Analysis_Report.txt - 详细分析\n');
    fprintf(fid, '• Performance_Comparison.png - 性能对比图\n');
    fprintf(fid, '• 原始数据和中间结果保存在 mcg_main_results.mat 中\n\n');
    
    fprintf(fid, '生成时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    
    fclose(fid);
    fprintf('       使用说明文件已生成\n');
    
catch ME
    fclose(fid);
    fprintf('       使用说明生成失败: %s\n', ME.message);
end

end

function display_comprehensive_mcg_summary(validation_results, data_info)
% 显示综合MCG摘要

fprintf('\n=================== MCG系统辨识综合摘要 ===================\n');

%% ====== 数据概况 ======
fprintf('\ 数据概况:\n');
fprintf('   数据源: %s\n', data_info.data_folder);
fprintf('   数据集: %d 个 | 选择策略: %s\n', ...
        validation_results.data_summary.total_datasets, data_info.selection_strategy);

if isfield(validation_results.data_summary, 'quality_stats')
    stats = validation_results.data_summary.quality_stats;
    fprintf('   数据质量: %.3f ± %.3f (范围: [%.3f, %.3f])\n', ...
            stats.mean, stats.std, stats.min, stats.max);
end

%% ====== 方法性能对比 ======
fprintf('\n 方法性能对比:\n');

id_perf = validation_results.identification_performance;
theory_perf = validation_results.theory_performance;
hybrid_perf = validation_results.hybrid_performance;

fprintf('   ┌─────────────────────┬─────────────┬─────────────┬─────────────┐\n');
fprintf('   │      方法           │   相关性    │   物理评分  │   综合评价  │\n');
fprintf('   ├─────────────────────┼─────────────┼─────────────┼─────────────┤\n');
fprintf('   │ 数据驱动辨识        │    %.3f    │    %.2f     │     %s     │\n', ...
        id_perf.mean_correlation, id_perf.physics_check.overall_score, ...
        pad_string(id_perf.overall_performance, 7));
fprintf('   │ MuJoCo理论模型      │    %.3f    │    %.2f     │   理论基础   │\n', ...
        theory_perf.estimated_correlation, theory_perf.physics_check.overall_score);
fprintf('   │ 智能混合方法       │    %.3f    │    %.2f     │     %s     │\n', ...
        hybrid_perf.predicted_correlation, hybrid_perf.physics_check.overall_score, ...
        pad_string(hybrid_perf.overall_performance, 7));
fprintf('   └─────────────────────┴─────────────┴─────────────┴─────────────┘\n');

%% ====== 性能改进分析 ======
fprintf('\n📈 性能改进分析:\n');
comp = validation_results.comparison_analysis;

fprintf('   相比数据驱动辨识: %+.1f%%\n', comp.improvement_vs_identification * 100);
fprintf('   相比理论模型: %+.1f%%\n', comp.improvement_vs_theory * 100);
fprintf('   相比最佳单一方法: %+.1f%%\n', comp.improvement_vs_best_single * 100);

if ~isempty(comp.advantages)
    fprintf('   主要优势: ');
    for i = 1:length(comp.advantages)
        if i > 1, fprintf(', '); end
        fprintf('%s', comp.advantages{i});
    end
    fprintf('\n');
end

%% ====== 物理验证结果 ======
fprintf('\n⚖️ 物理验证结果:\n');
phys = validation_results.physics_validation;

status_symbols = {'no', 'yes'};
fprintf('   动能正定性: %s | 科氏反对称: %s | 质量尺度: %s\n', ...
        status_symbols{phys.kinetic_energy_positive + 1}, ...
        status_symbols{phys.coriolis_antisymmetric + 1}, ...
        status_symbols{phys.mass_scale_check + 1});

fprintf('   物理一致性: %.2f (%s)\n', phys.consistency_score, phys.consistency_rating);

%% ====== 数值稳定性 ======
fprintf('\n 数值稳定性:\n');
num = validation_results.numerical_analysis;

fprintf('   条件数: 辨识(%.1e) | 理论(%.1e) | 混合(%.1e)\n', ...
        num.condition_numbers.identification, num.condition_numbers.theory, num.condition_numbers.hybrid);
fprintf('   数值稳定性: %.2f (%s)\n', num.stability_score, num.stability_rating);

%% ====== 综合评价和推荐 ======
fprintf('\n 综合评价:\n');
fprintf('   综合评分: %.3f / 1.000\n', validation_results.overall_score);
fprintf('   推荐等级: %s\n', validation_results.recommendation.level);

if ~isempty(validation_results.recommendation.reasons)
    fprintf('   推荐理由: ');
    for i = 1:length(validation_results.recommendation.reasons)
        if i > 1, fprintf(', '); end
        fprintf('%s', validation_results.recommendation.reasons{i});
    end
    fprintf('\n');
end

%% ====== 使用建议 ======
if isfield(validation_results.recommendation, 'usage_suggestions') && ~isempty(validation_results.recommendation.usage_suggestions)
    fprintf('\n 使用建议:\n');
    for i = 1:length(validation_results.recommendation.usage_suggestions)
        fprintf('   • %s\n', validation_results.recommendation.usage_suggestions{i});
    end
end

%% ====== 结论 ======
fprintf('\n 结论:\n');

if validation_results.overall_score > 0.8
    fprintf('    MCG混合辨识取得优秀效果，强烈推荐使用\n');
    fprintf('    混合方法在性能和物理一致性方面均表现出色\n');
    
elseif validation_results.overall_score > 0.6
    fprintf('    MCG混合辨识取得良好效果，推荐使用\n');
    fprintf('    建议在应用前进行进一步验证\n');
    
elseif validation_results.overall_score > 0.4
    fprintf('    MCG混合辨识取得中等效果，需要谨慎使用\n');
    fprintf('    建议改进数据质量或调整融合策略\n');
    
else
    fprintf('    MCG混合辨识效果不理想，不推荐使用\n');
    fprintf('    建议重新采集数据或改进辨识方法\n');
end

fprintf('\n  结果文件已保存，请查看详细报告和矩阵文件\n');
fprintf('==========================================================\n\n');

end

%% ====== 辅助函数 ======
function file_size_mb = get_file_size_mb(filename)
% 获取文件大小（MB）

try
    file_info = dir(filename);
    file_size_mb = file_info.bytes / (1024 * 1024);
catch
    file_size_mb = 0;
end

end

function str = bool_to_string(bool_val)
% 布尔值转字符串

if bool_val
    str = '是';
else
    str = '否';
end

end

function padded_str = pad_string(input_str, target_length)
% 字符串填充到指定长度

current_length = length(input_str);
if current_length >= target_length
    padded_str = input_str(1:target_length);
else
    padding = repmat(' ', 1, target_length - current_length);
    padded_str = [input_str, padding];
end

end