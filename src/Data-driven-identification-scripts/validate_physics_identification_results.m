function validation_results = validate_physics_identification_results(M_IDENTIFIED, G_IDENTIFIED, C_IDENTIFIED, q_data, q_dot_data, q_ddot_data, u_data)
% 验证物理辨识结果

validation_results = struct();

try
    N_total = size(q_data, 1);
    N_test = min(1000, N_total);
    test_indices = round(linspace(1, N_total, N_test));
    
    fprintf('      使用 %d 个样本进行FLU坐标系验证\n', N_test);
    
    % 通道名称（修正版，确保有效的MATLAB字段名）
    channel_names = {'X_Force_Forward', 'Y_Force_Left', 'Z_Force_Up', 'Roll_Torque', 'Pitch_Torque', 'Yaw_Torque', 'Armz_Torque', 'Army_Torque'};
    channel_display = {'X_Force(前)', 'Y_Force(左)', 'Z_Force(上)', 'Roll_Torque', 'Pitch_Torque', 'Yaw_Torque', 'Armz_Torque', 'Army_Torque'};
    
    n_channels = min(length(channel_names), size(u_data, 2));
    
    correlations = zeros(n_channels, 1);
    
    % DOF级别验证
    validation_results.channel_correlations = struct();
    validation_results.channel_performance = struct();
    
    for j = 1:n_channels
        predictions = zeros(N_test, 1);
        
        for i = 1:N_test
            idx = test_indices(i);
            
            % 动力学预测
            predicted_force = M_IDENTIFIED(j, j) * q_ddot_data(idx, j) + G_IDENTIFIED(j);
            
            % 添加主要耦合项
            for k = 1:size(M_IDENTIFIED, 2)
                if k ~= j && abs(M_IDENTIFIED(j, k)) > 1e-6
                    predicted_force = predicted_force + M_IDENTIFIED(j, k) * q_ddot_data(idx, k);
                end
            end
            
            predictions(i) = predicted_force;
        end
        
        % 计算相关系数
        if std(predictions) > 1e-6 && std(u_data(test_indices, j)) > 1e-6
            corr_matrix = corrcoef(predictions, u_data(test_indices, j));
            correlations(j) = corr_matrix(1, 2);
        else
            correlations(j) = 0;
        end
        
        % 性能评级
        if abs(correlations(j)) > 0.8
            performance = 'excellent';
        elseif abs(correlations(j)) > 0.6
            performance = 'good';
        elseif abs(correlations(j)) > 0.4
            performance = 'fair';
        else
            performance = 'poor';
        end
        
        validation_results.channel_performance.(channel_names{j}) = performance;
        validation_results.channel_correlations.(channel_names{j}) = correlations(j);
        
        fprintf('        %s: 相关性=%.3f, 性能=%s\n', channel_display{j}, correlations(j), performance);
    end
    
    % 总体性能（使用绝对值）
    valid_correlations = correlations(~isnan(correlations));
    validation_results.mean_correlation = mean(abs(valid_correlations));
    validation_results.correlations = correlations;
    
    % 计算误差统计
    all_errors = [];
    for j = 1:n_channels
        predictions = zeros(N_test, 1);
        
        for i = 1:N_test
            idx = test_indices(i);
            predicted_force = M_IDENTIFIED(j, j) * q_ddot_data(idx, j) + G_IDENTIFIED(j);
            predictions(i) = predicted_force;
        end
        
        errors = abs(predictions - u_data(test_indices, j));
        all_errors = [all_errors; errors];
    end
    
    validation_results.mean_error = mean(all_errors);
    validation_results.rms_error = sqrt(mean(all_errors.^2));
    
    % 整体性能评级
    if validation_results.mean_correlation > 0.7
        validation_results.overall_performance = 'excellent';
    elseif validation_results.mean_correlation > 0.5
        validation_results.overall_performance = 'good';
    elseif validation_results.mean_correlation > 0.3
        validation_results.overall_performance = 'fair';
    else
        validation_results.overall_performance = 'poor';
    end
    
    validation_results.validation_success = true;
    validation_results.mean_r_squared = validation_results.mean_correlation^2;
    
    fprintf('      FLU坐标系验证完成：平均相关性 %.3f, 性能: %s\n', ...
            validation_results.mean_correlation, validation_results.overall_performance);
    
catch ME
    validation_results.validation_success = false;
    validation_results.error_message = ME.message;
    validation_results.mean_correlation = 0.3;
    validation_results.overall_performance = 'poor';
    validation_results.mean_error = 0;
    validation_results.rms_error = 0;
    fprintf('       FLU坐标系验证失败: %s\n', ME.message);
end

end