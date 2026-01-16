% it_ppi_find_optimal_params.m
% このスクリプトは、PPIベースのIT政策の最適パラメータを探索し、
% 厚生の2次元線グラフを出力します。計算の核心部は src/it_ppi_solve_core.m を呼び出します。

clear;
close all;
clc;

% --- [1] パス設定 ---
mpath = fileparts(mfilename('fullpath'));
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_;

fprintf('====================================================\n');
fprintf('Starting optimization for IT_PPI Policy...\n');
fprintf('====================================================\n\n');

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- 設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼

% --- 1. ショック・シナリオを選択 ---
% 'beta' : 需要ショック (eps_beta_H = 0.02)
% 'a'    : 供給ショック (eps_a_H = -0.1)
shock_scenario = 'beta'; 

% --- フォルダ用ショック名の決定 ---
if strcmp(shock_scenario, 'beta')
    shock_folder = 'beta_shock';
else
    shock_folder = 'a_shock';
end

% --- 2. ショックパラメータの設定 ---
if strcmp(shock_scenario, 'beta')
    target_shock = 'eps_beta_H';
    shock_value  = 0.02; 
elseif strcmp(shock_scenario, 'a')
    target_shock = 'eps_a_H';
    shock_value  = -0.1; 
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% --- 3. 探索範囲の設定 (1次元ラインサーチ用) --- beta final
% 開始する数字 : 間隔 : 終わりの数字
% IT-PPIではインフレ反応係数 phi_pi_H のみを最適化します
phi_pi_H_range = 2 : 10 : 300;

% --- 3. 探索範囲の設定 (1次元ラインサーチ用) --- a final
% 開始する数字 : 間隔 : 終わりの数字
% IT-PPIではインフレ反応係数 phi_pi_H のみを最適化します
%phi_pi_H_range = 1 : 0.005 : 1.02; 

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

% --- [修正] 出力フォルダの定義と作成 ---
% 構成: results/optimal_params/[shock_folder]/ 直下にすべて出力
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- Simulation Settings Structure ---
% 命名規則に基づき先頭に it_ppi を配置したモデルファイルを参照
settings.modnam = 'it_ppi_model';
settings.modnamstar = 'it_ppi_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.irfshock = target_shock; 
settings.shockssequence = [ shock_value; zeros(399, 1) ];
settings.nperiods = 400;
settings.maxiter = 100;

% =========================================================================
%% --- Main Calculation Loop (1D Line Search) ---
% =========================================================================
fprintf('--- Starting 1D Line search for phi_pi_H (PPI) for %s shock ---\n', shock_scenario);

% 結果を格納するベクトルを初期化
results_welfare_H   = zeros(length(phi_pi_H_range), 1);
results_welfare_F   = zeros(length(phi_pi_H_range), 1);
results_welfare_SUM = zeros(length(phi_pi_H_range), 1);

total_iterations = length(phi_pi_H_range);
root_dir = pwd;

% 1次元の探索ループ
for i = 1:total_iterations
    phi_pi_H_value = phi_pi_H_range(i);

    fprintf('\nTrial (%d/%d): phi_pi_H = %.4f ... ', i, total_iterations, phi_pi_H_value);

    % --- [計算核心 it_ppi_solve_core を関数として呼び出す] ---
    try
        % 命名規則を統一したコアエンジンを呼び出し
        res = it_ppi_solve_core(phi_pi_H_value, settings);

        % 命名規則 (with_delta_pie) の整合性をとって結果を抽出
        welfare_H   = res.welfare_H_with_delta_pie;
        welfare_F   = res.welfare_F_with_delta_pie;
        welfare_SUM = res.welfare_SUM_with_delta_pie;
    catch ME
        fprintf('Error during trial: %s\n', ME.message);
        cd(root_dir);
        welfare_H = -999; welfare_F = -999; welfare_SUM = -999;
    end
    % --- [呼び出し終了] ---

    results_welfare_H(i)   = welfare_H;
    results_welfare_F(i)   = welfare_F;
    results_welfare_SUM(i) = welfare_SUM;

    if welfare_SUM > -999
        fprintf('Welfare H: %.4f, F: %.4f, Sum: %.4f\n', welfare_H, welfare_F, welfare_SUM);
    else
        fprintf('Calculation Error.\n');
    end
end

if max(results_welfare_SUM(:)) <= -999
    fprintf('--- SEARCH FAILED: All parameter combinations resulted in an error. ---\n');
    return;
end

% =========================================================================
%% --- Find and Package Optimal Parameters ---
% =========================================================================
% 1. 自国（Home）の厚生を最大化するパラメータ
[max_welfare_H, idx_H] = max(results_welfare_H);
final_results.opt_phi_pi_H_H = phi_pi_H_range(idx_H);
final_results.max_welfare_H = max_welfare_H;
final_results.welfare_F_at_H_opt = results_welfare_F(idx_H);
final_results.welfare_SUM_at_H_opt = results_welfare_SUM(idx_H);

% 2. 外国（Foreign）の厚生を最大化するパラメータ
[max_welfare_F, idx_F] = max(results_welfare_F);
final_results.opt_phi_pi_H_F = phi_pi_H_range(idx_F);
final_results.max_welfare_F = max_welfare_F;
final_results.welfare_H_at_F_opt = results_welfare_H(idx_F);
final_results.welfare_SUM_at_F_opt = results_welfare_SUM(idx_F);

% 3. 全世界（Global Sum）の厚生を最大化するパラメータ
[max_welfare_SUM, idx_SUM] = max(results_welfare_SUM);
final_results.opt_phi_pi_H_SUM = phi_pi_H_range(idx_SUM);
final_results.max_welfare_SUM = max_welfare_SUM;
final_results.welfare_H_at_SUM_opt = results_welfare_H(idx_SUM);
final_results.welfare_F_at_SUM_opt = results_welfare_F(idx_SUM);

% 全データの保存 (all_data.mat)
save(fullfile(base_res_dir, 'all_data.mat'), ...
    'final_results', 'phi_pi_H_range', 'results_welfare_H', 'results_welfare_F', 'results_welfare_SUM');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n--- IT_PPI Final Optimal Policy Results (Shock: %s) ---\n', shock_scenario);
fprintf('1. Optimal for Home:\n');
fprintf('   Parameter: phi_pi_H = %.4f\n', final_results.opt_phi_pi_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameter: phi_pi_H = %.4f\n', final_results.opt_phi_pi_H_F);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);
fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameter: phi_pi_H = %.4f\n', final_results.opt_phi_pi_H_SUM);
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
h_home = figure('Name', 'IT_PPI: Home Welfare', 'Visible', 'off');
plot(phi_pi_H_range, results_welfare_H, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
plot(final_results.opt_phi_pi_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
grid on;
xlabel('\phi_{\pi}^{H}'); ylabel('Welfare');
line1_H = sprintf('IT_{PPI}: Optimal for Home (%s shock)', shock_scenario);
line2_H = sprintf('\\phi_{\\pi}^{H}=%.4f, Welfare=%.4f', final_results.opt_phi_pi_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_home, fullfile(base_res_dir, 'graph_home.pdf'), '-dpdf');
close(h_home);

% Home個別データ
data_x = phi_pi_H_range; data_y = results_welfare_H;
opt_pt = [final_results.opt_phi_pi_H_H, final_results.max_welfare_H];
save(fullfile(base_res_dir, 'data_for_plotting_home.mat'), 'data_x', 'data_y', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'IT_PPI: Foreign Welfare', 'Visible', 'off');
plot(phi_pi_H_range, results_welfare_F, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
plot(final_results.opt_phi_pi_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
grid on;
xlabel('\phi_{\pi}^{H}'); ylabel('Welfare');
line1_F = sprintf('IT_{PPI}: Optimal for Foreign (%s shock)', shock_scenario);
line2_F = sprintf('\\phi_{\\pi}^{H}=%.4f, Welfare=%.4f', final_results.opt_phi_pi_H_F, final_results.max_welfare_F);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_foreign, fullfile(base_res_dir, 'graph_foreign.pdf'), '-dpdf');
close(h_foreign);

% Foreign個別データ
data_x = phi_pi_H_range; data_y = results_welfare_F;
opt_pt = [final_results.opt_phi_pi_H_F, final_results.max_welfare_F];
save(fullfile(base_res_dir, 'data_for_plotting_foreign.mat'), 'data_x', 'data_y', 'opt_pt');

% --- Graph 3: Global (Sum) Welfare ---
h_global = figure('Name', 'IT_PPI: Global Welfare', 'Visible', 'off');
plot(phi_pi_H_range, results_welfare_SUM, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
plot(final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
grid on;
xlabel('\phi_{\pi}^{H}'); ylabel('Welfare');
line1_G = sprintf('IT_{PPI}: Optimal for Global (Sum) (%s shock)', shock_scenario);
line2_G = sprintf('\\phi_{\\pi}^{H}=%.4f, Welfare=%.4f', final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_global, fullfile(base_res_dir, 'graph_global.pdf'), '-dpdf');
close(h_global);

% Global個別データ
data_x = phi_pi_H_range; data_y = results_welfare_SUM;
opt_pt = [final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM];
save(fullfile(base_res_dir, 'data_for_plotting_global.mat'), 'data_x', 'data_y', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('IT_PPI Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');