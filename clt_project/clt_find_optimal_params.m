% clt_find_optimal_params.m
% このスクリプトは、CLT政策の最適パラメータを探索し、
% 厚生の表面グラフを出力します。計算の核心部は src/clt_solve_core.m を呼び出します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_;

fprintf('====================================================\n');
fprintf('Starting optimization for CLT Policy...\n');
fprintf('====================================================\n\n');

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- 設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼

% --- 1. 探索モードを選択 ---
% 'coarse' : 広範囲を粗く探索します
% 'fine'   : 指定した座標の周辺を詳細に探索します
search_mode = 'coarse'; 

% --- 2. ショック・シナリオを選択 ---
% 'beta' : 需要ショック (eps_beta_H = 0.02)
% 'a'    : 供給ショック (eps_a_H = -0.1)
shock_scenario = 'beta'; 

% --- フォルダ用ショック名の決定 ---
if strcmp(shock_scenario, 'beta')
    shock_folder = 'beta_shock';
else
    shock_folder = 'a_shock';
end

% --- 3. ショックパラメータの設定 ---
if strcmp(shock_scenario, 'beta')
    target_shock = 'eps_beta_H';
    shock_value  = 0.02; 
elseif strcmp(shock_scenario, 'a')
    target_shock = 'eps_a_H';
    shock_value  = -0.1; 
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% --- 4. 粗い探索 (coarse mode) の設定 ---
%phi_gap_H_range_coarse  = 0:0.001:0.01;
%phi_level_H_range_coarse = 0:0.000001:0.000003;

% --- 5. 詳細な探索 (fine mode) の設定 --- beta shock final
%center_phi_gap_H  = 0.17;
%center_phi_level_H = 5.06;
%fine_search_radius_gap_H   = 0.01;
%fine_search_step_gap_H     = 0.01;
%fine_search_radius_level_H = 0.01;
%fine_search_step_level_H   = 0.01;

% --- 5. 詳細な探索 (fine mode) の設定 --- a shock final
%center_phi_gap_H  = 0.008;
%center_phi_level_H = 0.0000001;
%fine_search_radius_gap_H   = 0.002;
%fine_search_step_gap_H     = 0.0005;
%fine_search_radius_level_H = 0.0000001;
%fine_search_step_level_H   = 0.0000001;

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

% --- 出力フォルダの定義と作成 ---
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- Simulation Settings Structure ---
% 命名規則に基づき先頭に clt を配置したモデルファイルを参照
settings.modnam = 'clt_model';
settings.modnamstar = 'clt_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.irfshock = target_shock;
settings.shockssequence = [ shock_value; zeros(399, 1) ];
settings.nperiods = 400;
settings.maxiter = 100;

% =========================================================================
%% --- Grid Setup based on Mode ---
% =========================================================================
if strcmp(search_mode, 'coarse')
    fprintf('--- Starting COARSE grid search ---\n');
    phi_gap_H_range_final   = phi_gap_H_range_coarse;
    phi_level_H_range_final = phi_level_H_range_coarse;
    plot_title_prefix = 'Coarse Search';
elseif strcmp(search_mode, 'fine')
    fprintf('--- Starting FINE grid search around (gap_H=%.4f, level_H=%.8f) ---\n', center_phi_gap_H, center_phi_level_H);
    phi_gap_H_range_final   = (center_phi_gap_H - fine_search_radius_gap_H):fine_search_step_gap_H:(center_phi_gap_H + fine_search_radius_gap_H);
    phi_level_H_range_final = (center_phi_level_H - fine_search_radius_level_H):fine_search_step_level_H:(center_phi_level_H + fine_search_radius_level_H);
    phi_gap_H_range_final(phi_gap_H_range_final <= 0) = [];
    phi_level_H_range_final(phi_level_H_range_final <= 0) = [];
    plot_title_prefix = 'Fine Search';
else
    error("Invalid search_mode. Choose 'coarse' or 'fine'.");
end

if isempty(phi_gap_H_range_final) || isempty(phi_level_H_range_final)
    error("The parameter range for the search is empty. Check your settings.");
end

% =========================================================================
%% --- Main Calculation Loop ---
% =========================================================================
results_welfare_H_final   = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));
results_welfare_F_final   = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));
results_welfare_SUM_final = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));

total_iterations = length(phi_level_H_range_final) * length(phi_gap_H_range_final);
current_iteration = 0;

root_dir = pwd;

for i = 1:length(phi_level_H_range_final)
    for j = 1:length(phi_gap_H_range_final)
        current_iteration = current_iteration + 1;
        phi_gap_H_value = phi_gap_H_range_final(j);
        phi_level_H_value = phi_level_H_range_final(i);

        fprintf('\nTrial (%d/%d): phi_gap_H = %.4f, phi_level_H = %.8f ... ', ...
            current_iteration, total_iterations, phi_gap_H_value, phi_level_H_value);

        % --- [計算核心 solve_core を関数として呼び出す] ---
        try
            % 命名規則を統一した clt_solve_core を呼び出し
            res = clt_solve_core(phi_gap_H_value, phi_level_H_value, settings);

            % 抽出
            welfare_H   = res.welfare_H_with_delta_pie;
            welfare_F   = res.welfare_F_with_delta_pie;
            welfare_SUM = res.welfare_SUM_with_delta_pie;
        catch ME
            fprintf('Error during trial: %s\n', ME.message);
            cd(root_dir);
            welfare_H   = -999; welfare_F   = -999; welfare_SUM = -999;
        end
        % --- [呼び出し終了] ---

        results_welfare_H_final(i, j)   = welfare_H;
        results_welfare_F_final(i, j)   = welfare_F;
        results_welfare_SUM_final(i, j) = welfare_SUM;

        if welfare_SUM > -999
            fprintf('Welfare H: %.4f, F: %.4f, Sum: %.4f\n', welfare_H, welfare_F, welfare_SUM);
        else
            fprintf('Calculation Error.\n');
        end
    end
end

if max(results_welfare_SUM_final(:)) <= -999
    fprintf('--- SEARCH FAILED: All parameter combinations resulted in an error. ---\n');
    return;
end

% =========================================================================
%% --- Find and Package Optimal Parameters ---
% =========================================================================
[max_welfare_H, idx_H] = max(results_welfare_H_final(:));
[row_H, col_H] = ind2sub(size(results_welfare_H_final), idx_H);
final_results.opt_phi_gap_H_H   = phi_gap_H_range_final(col_H);
final_results.opt_phi_level_H_H = phi_level_H_range_final(row_H);
final_results.max_welfare_H = max_welfare_H;
final_results.welfare_F_at_H_opt = results_welfare_F_final(row_H, col_H);
final_results.welfare_SUM_at_H_opt = results_welfare_SUM_final(row_H, col_H);

[max_welfare_F, idx_F] = max(results_welfare_F_final(:));
[row_F, col_F] = ind2sub(size(results_welfare_F_final), idx_F);
final_results.opt_phi_gap_H_F   = phi_gap_H_range_final(col_F);
final_results.opt_phi_level_H_F = phi_level_H_range_final(row_F);
final_results.max_welfare_F = max_welfare_F;
final_results.welfare_H_at_F_opt = results_welfare_H_final(row_F, col_F);
final_results.welfare_SUM_at_F_opt = results_welfare_SUM_final(row_F, col_F);

[max_welfare_SUM, idx_SUM] = max(results_welfare_SUM_final(:));
[row_SUM, col_SUM] = ind2sub(size(results_welfare_SUM_final), idx_SUM);
final_results.opt_phi_gap_H_SUM   = phi_gap_H_range_final(col_SUM);
final_results.opt_phi_level_H_SUM = phi_level_H_range_final(row_SUM);
final_results.max_welfare_SUM = max_welfare_SUM;
final_results.welfare_H_at_SUM_opt = results_welfare_H_final(row_SUM, col_SUM);
final_results.welfare_F_at_SUM_opt = results_welfare_F_final(row_SUM, col_SUM);

% 全データの保存 (all_data.mat)
save(fullfile(base_res_dir, 'all_data.mat'), ...
    'final_results', 'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'phi_gap_H_range_final', 'phi_level_H_range_final');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n--- CLT Final Optimal Policy Results ---\n');
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_gap_H = %.4f, phi_level_H = %.8f\n', final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameters: phi_gap_H = %.4f, phi_level_H = %.8f\n', final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);
fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_gap_H = %.4f, phi_level_H = %.8f\n', final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n', final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('--------------------------------------------------------------\n\n');

% =========================================================================
%% --- Plot and Save Graphs & Raw Data ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

title_pos_line1 = [0, 0.96, 1, 0.05];
title_pos_line2 = [0, 0.91, 1, 0.05];

% --- Graph 1: Home Welfare ---
h_home = figure('Name', 'CLT: Home Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_H_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_H = zlim;
plot3([final_results.opt_phi_gap_H_H, final_results.opt_phi_gap_H_H], ...
      [final_results.opt_phi_level_H_H, final_results.opt_phi_level_H_H], ...
      [z_limits_H(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
line1_H = sprintf('%s: Optimal for Home (%s)', plot_title_prefix, shock_scenario);
line2_H = sprintf('\\phi_{gap}^H=%.4f, \\phi_{level}^H=%.8f, Welfare=%.4f', ...
    final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_home, fullfile(base_res_dir, 'graph_home.pdf'), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_gap_H_range_final; data_y = phi_level_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H];
save(fullfile(base_res_dir, 'data_for_plotting_home.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'CLT: Foreign Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_F_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_F = zlim;
plot3([final_results.opt_phi_gap_H_F, final_results.opt_phi_gap_H_F], ...
      [final_results.opt_phi_level_H_F, final_results.opt_phi_level_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
line1_F = sprintf('%s: Optimal for Foreign (%s)', plot_title_prefix, shock_scenario);
line2_F = sprintf('\\phi_{gap}^H=%.4f, \\phi_{level}^H=%.8f, Welfare=%.4f', ...
    final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_foreign, fullfile(base_res_dir, 'graph_foreign.pdf'), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F];
save(fullfile(base_res_dir, 'data_for_plotting_foreign.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global (Sum) Welfare ---
h_global = figure('Name', 'CLT: Global Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_SUM_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_G = zlim;
plot3([final_results.opt_phi_gap_H_SUM, final_results.opt_phi_gap_H_SUM], ...
      [final_results.opt_phi_level_H_SUM, final_results.opt_phi_level_H_SUM], ...
      [z_limits_G(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
line1_G = sprintf('%s: Optimal for Global (Sum) (%s)', plot_title_prefix, shock_scenario);
line2_G = sprintf('\\phi_{gap}^H=%.4f, \\phi_{level}^H=%.8f, Welfare=%.4f', ...
    final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_global, fullfile(base_res_dir, 'graph_global.pdf'), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM];
save(fullfile(base_res_dir, 'data_for_plotting_global.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('CLT Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');