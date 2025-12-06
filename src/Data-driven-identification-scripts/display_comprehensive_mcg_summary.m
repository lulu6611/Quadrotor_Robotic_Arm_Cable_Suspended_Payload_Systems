function display_comprehensive_mcg_summary(validation_results, data_info)
% 显示综合MCG摘要

fprintf('\n=================== MCG系统辨识综合摘要 ===================\n');

%% ====== 数据概况 ======
fprintf(' 数据概况:\n');
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
fprintf('   │ 智能混合方法      │    %.3f    │    %.2f     │     %s     │\n', ...
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
fprintf('\n 物理验证结果:\n');
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
fprintf('\n结论:\n');

if validation_results.overall_score > 0.8
    fprintf('   MCG混合辨识取得优秀效果，强烈推荐使用\n');
    fprintf('   混合方法在性能和物理一致性方面均表现出色\n');
    
elseif validation_results.overall_score > 0.6
    fprintf('   MCG混合辨识取得良好效果，推荐使用\n');
    fprintf('   建议在应用前进行进一步验证\n');
    
elseif validation_results.overall_score > 0.4
    fprintf('   MCG混合辨识取得中等效果，需要谨慎使用\n');
    fprintf('   建议改进数据质量或调整融合策略\n');
    
else
    fprintf('  MCG混合辨识效果不理想，不推荐使用\n');
    fprintf('   建议重新采集数据或改进辨识方法\n');
end

fprintf('\n 结果文件已保存，请查看详细报告和矩阵文件\n');
fprintf('==========================================================\n\n');

end

%% ====== 辅助函数 ======
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