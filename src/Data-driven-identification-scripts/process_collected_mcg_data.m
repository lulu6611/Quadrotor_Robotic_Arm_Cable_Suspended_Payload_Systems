function process_collected_mcg_data()
% 处理采集的MCG数据，执行完整的辨识流程

fprintf('========== MCG数据处理与辨识流程 ==========\n');
fprintf('开始时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ====== 第一步：数据加载和选择 ======
fprintf('\n 第一步：智能数据加载与选择...\n');
[selected_data, data_info] = smart_data_loader_and_selector();

if isempty(selected_data)
    error('未找到有效数据，请检查数据文件夹');
end

%% ====== 第二步：基于物理正确性的MCG辨识 ======
fprintf('\n 第二步：基于物理正确性的MCG辨识...\n');
[M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results] = perform_physics_based_mcg_identification(selected_data, data_info);

%% ====== 第三步：MuJoCo理论模型构建 ======
fprintf('\n 第三步：MuJoCo理论模型构建...\n');
[M_THEORY, G_THEORY, C_THEORY, theory_info] = build_mujoco_theoretical_model();

%% ====== 第四步：智能混合辨识 ======
fprintf('\n 第四步：智能混合MCG辨识...\n');
[M_HYBRID, G_HYBRID, C_HYBRID, fusion_results] = perform_intelligent_mcg_fusion(...
    M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results, ...
    M_THEORY, G_THEORY, C_THEORY, theory_info);

%% ====== 第五步：结果验证与分析 ======
fprintf('\n 第五步：结果验证与性能分析...\n');
validation_results = validate_and_analyze_mcg_results(...
    M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results, ...
    M_THEORY, G_THEORY, C_THEORY, theory_info, ...
    M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, ...
    selected_data, data_info);

%% ====== 第六步：保存和显示结果 ======
fprintf('\n 第六步：保存结果与生成报告...\n');
save_mcg_processing_results(M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, id_results, ...
                           M_THEORY, G_THEORY, C_THEORY, theory_info, ...
                           M_HYBRID, G_HYBRID, C_HYBRID, fusion_results, ...
                           validation_results, selected_data, data_info);

display_comprehensive_mcg_summary(validation_results, data_info);

fprintf('\n========== MCG处理流程完成 ==========\n');
fprintf(' 完成时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end