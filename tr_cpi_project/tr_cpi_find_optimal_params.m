% tr_cpi_find_optimal_params.m
% このスクリプトは、TR CPI政策の最適パラメータを探索し、
% 厚生の表面グラフを出力します。計算の核心部は stage2/tr_cpi_stage2_solve_core.m を呼び出し、
% インフレ反応係数(phi_pi_H)と産出反応係数(phi_y_H)を最適化します。
% 構成を POLT の find スクリプトと完全に統一し、以下の手順を全自動で行います：
% 1. Stage 1 で伸縮価格下の潜在産出量パスを自動算出（目標値の固定）
% 2. Stage 2 (OccBin) エンジンを呼び出し、phi_pi_H と phi_y_H を最適化
% 3. 自国・外国・全世界のそれぞれの最適値を特定し、グラフとデータを詳細に保存

clear;
close all;
clc;

% --- [1] パス設定 (POLTと同一のルート基準構成) ---
mpath = fileparts(mfilename('fullpath'));
cd(mpath); % スクリプトの場所（Root）をカレントディレクトリに固定

% プロジェクト内の各フォルダへのパスを追加
addpath(fullfile(mpath, 'src', 'stage2'));
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_ options_;

fprintf('====================================================\n');
fprintf('Starting optimization for TR CPI Policy (ZLB Aware)...\n');
fprintf('====================================================\n\n');

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- 設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼

% --- 1. 探索モードを選択 ---
% 'coarse' : 広範囲を粗く探索
% 'fine'   : 特定座標の周辺を詳細に探索
search_mode = 'coarse'; 

% --- 2. ショック・シナリオを選択 ---
% 'beta' : 需要ショック (eps_beta_H = 0.02)
% 'a'    : 供給ショック (eps_a_H = -0.1)
if ~exist('shock_scenario', 'var')
    shock_scenario = 'beta'; 
end

% --- フォルダ用ショック名の決定 ---
if strcmp(shock_scenario, 'beta')
    shock_folder = 'beta_shock';
    target_shock = 'eps_beta_H';
    primary_shock_magnitude = 0.02; 
elseif strcmp(shock_scenario, 'a')
    shock_folder = 'a_shock';
    target_shock = 'eps_a_H';
    primary_shock_magnitude = -0.1; 
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% --- 3. 粗い探索 (coarse mode) の設定 --- beta final
%phi_pi_H_range_coarse   = 82:1:84;
%phi_y_H_range_coarse    = 0.4:0.1:0.6;

% --- 3. 粗い探索 (coarse mode) の設定 --- a final
%phi_pi_H_range_coarse   = 14000000000000:1000000000000:16000000000000;
%phi_y_H_range_coarse    = 0.0000000006:0.0000000001:0.0000000008;

% --- 4. 詳細な探索 (fine mode) の設定 ---
%center_phi_pi_H   = 1.5; 
%center_phi_y_H    = 0.5; 
%fine_search_radius_pi_H   = 0.4; 
%fine_search_step_pi_H     = 0.1; 
%fine_search_radius_y_H    = 0.4; 
%fine_search_step_y_H      = 0.1;

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

% --- 出力フォルダの定義と作成 ---
% 構成: results/optimal_params/[shock_folder]/ 直下にすべて出力
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- [1] Stage 1 (潜在産出量) パスの自動計算 (POLT構成と統一) ---
fprintf('\n--- Stage 1: Calculating Target Potential Path for %s shock ---\n', shock_scenario);

root_dir = pwd;
cd(fullfile(mpath, 'src', 'stage1'));

% 自然産出量の算出 (政策中立な伸縮価格モデルを実行)
% 【修正箇所】変数名を phi_pi_H_value, phi_y_H_value に統一
phi_pi_H_dummy = 1.5; phi_y_H_dummy = 0.5;
assignin('base', 'phi_pi_H_value', phi_pi_H_dummy);
assignin('base', 'phi_y_H_value', phi_y_H_dummy);

% Stage 1 モデルの実行
eval(['dynare tr_cpi_stage1_model.mod noclearall nolog']);

% ショックパスの計算 (simult_ を使用)
ex_flex = zeros(400, M_.exo_nbr);
s_idx = strmatch(target_shock, M_.exo_names, 'exact');
ex_flex(1, s_idx) = primary_shock_magnitude;
y_flex_sim = simult_(M_, options_, oo_.dr.ys, oo_.dr, ex_flex, 1);
y_H_idx = strmatch('y_H', M_.endo_names, 'exact');
eps_y_H_potential_vec = y_flex_sim(y_H_idx, 2:end)' - oo_.dr.ys(y_H_idx);

% 元のディレクトリ（Root）に戻る
cd(root_dir);

% 計算したパスを保存
save(fullfile(base_res_dir, 'stage1_potential_path.mat'), 'eps_y_H_potential_vec');
fprintf('Stage 1 完了。Stage 2 最適化を開始します...\n');

%% --- [2] Simulation Settings Structure ---
% Stage 2 用の非線形・ZLB考慮設定
settings.modnam = 'tr_cpi_stage2_model';
settings.modnamstar = 'tr_cpi_stage2_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.nperiods = 400;
settings.maxiter = 100;

% ショック定義 (主ショック + Stage1から引き継いだ潜在産出量パス)
settings.irfshock = char(target_shock, 'eps_y_H_potential');

% ショックシーケンスの構築 (2列行列: 第1列は主ショック、第2列は潜在産出量パス)
shockssequence = zeros(settings.nperiods, 2);
shockssequence(1, 1) = primary_shock_magnitude;
len_pot = min(settings.nperiods, length(eps_y_H_potential_vec));
shockssequence(1:len_pot, 2) = eps_y_H_potential_vec(1:len_pot);
settings.shockssequence = shockssequence;

% =========================================================================
%% --- Grid Setup based on Mode ---
% =========================================================================
if strcmp(search_mode, 'coarse')
    fprintf('--- Starting COARSE grid search ---\n');
    phi_pi_H_range_final   = phi_pi_H_range_coarse;
    phi_y_H_range_final = phi_y_H_range_coarse;
    plot_title_prefix = 'Coarse Search';
elseif strcmp(search_mode, 'fine')
    fprintf('--- Starting FINE grid search around (pi_H=%.2f, y_H=%.2f) ---\n', center_phi_pi_H, center_phi_y_H);
    phi_pi_H_range_final   = (center_phi_pi_H - fine_search_radius_pi_H):fine_search_step_pi_H:(center_phi_pi_H + fine_search_radius_pi_H);
    phi_y_H_range_final = (center_phi_y_H - fine_search_radius_y_H):fine_search_step_y_H:(center_phi_y_H + fine_search_radius_y_H);
    
    % 制約の適用
    phi_pi_H_range_final(phi_pi_H_range_final <= 1) = []; % テイラー原則の維持
    phi_y_H_range_final(phi_y_H_range_final < 0) = []; 
    plot_title_prefix = 'Fine Search';
else
    error("Invalid search_mode. Choose 'coarse' or 'fine'.");
end

if isempty(phi_pi_H_range_final) || isempty(phi_y_H_range_final)
    error("The parameter range for the search is empty. Check your settings.");
end

% =========================================================================
%% --- Main Calculation Loop ---
% =========================================================================
results_welfare_H_final   = zeros(length(phi_y_H_range_final), length(phi_pi_H_range_final));
results_welfare_F_final   = zeros(length(phi_y_H_range_final), length(phi_pi_H_range_final));
results_welfare_SUM_final = zeros(length(phi_y_H_range_final), length(phi_pi_H_range_final));

total_iterations = length(phi_y_H_range_final) * length(phi_pi_H_range_final);
current_iteration = 0;
root_dir = pwd;

for i = 1:length(phi_y_H_range_final)
    for j = 1:length(phi_pi_H_range_final)
        current_iteration = current_iteration + 1;
        phi_pi_H_value = phi_pi_H_range_final(j);
        phi_y_H_value = phi_y_H_range_final(i);

        fprintf('\nTrial (%d/%d): phi_pi_H = %.2f, phi_y_H = %.2f ... ', ...
            current_iteration, total_iterations, phi_pi_H_value, phi_y_H_value);

        % --- [計算核心 tr_cpi_stage2_solve_core を呼び出す] ---
        try
            % stage2/tr_cpi_stage2_solve_core.m を呼び出し
            res = tr_cpi_stage2_solve_core(phi_pi_H_value, phi_y_H_value, settings);

            % エンジンの命名規則 (with_delta_pie) に基づいて抽出
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
% 1. Optimal for Home
[max_w_H, idx_H] = max(results_welfare_H_final(:));
[row_H, col_H] = ind2sub(size(results_welfare_H_final), idx_H);
final_results.opt_phi_pi_H_H = phi_pi_H_range_final(col_H);
final_results.opt_phi_y_H_H  = phi_y_H_range_final(row_H);
final_results.max_welfare_H  = max_w_H;
final_results.welfare_F_at_H_opt = results_welfare_F_final(row_H, col_H);
final_results.welfare_SUM_at_H_opt = results_welfare_SUM_final(row_H, col_H);

% 2. Optimal for Foreign
[max_w_F, idx_F] = max(results_welfare_F_final(:));
[row_F, col_F] = ind2sub(size(results_welfare_F_final), idx_F);
final_results.opt_phi_pi_H_F = phi_pi_H_range_final(col_F);
final_results.opt_phi_y_H_F  = phi_y_H_range_final(row_F);
final_results.max_welfare_F  = max_w_F;
final_results.welfare_H_at_F_opt = results_welfare_H_final(row_F, col_F);
final_results.welfare_SUM_at_F_opt = results_welfare_SUM_final(row_F, col_F);

% 3. Optimal for Global (Sum)
[max_w_SUM, idx_SUM] = max(results_welfare_SUM_final(:));
[row_SUM, col_SUM] = ind2sub(size(results_welfare_SUM_final), idx_SUM);
final_results.opt_phi_pi_H_SUM = phi_pi_H_range_final(col_SUM);
final_results.opt_phi_y_H_SUM  = phi_y_H_range_final(row_SUM);
final_results.max_welfare_SUM  = max_w_SUM;
final_results.welfare_H_at_SUM_opt = results_welfare_H_final(row_SUM, col_SUM);
final_results.welfare_F_at_SUM_opt = results_welfare_F_final(row_SUM, col_SUM);

% 全データの保存 (all_optimization_data.mat)
save(fullfile(base_res_dir, 'all_optimization_data.mat'), ...
    'final_results', 'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'phi_pi_H_range_final', 'phi_y_H_range_final', 'shock_scenario');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n--- TR CPI Final Optimal Policy Results (Shock: %s) ---\n', shock_scenario);
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_pi_H = %.2f, phi_y_H = %.2f\n', final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameters: phi_pi_H = %.2f, phi_y_H = %.2f\n', final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);
fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_pi_H = %.2f, phi_y_H = %.2f\n', final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM);
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
h_home = figure('Name', 'TR CPI: Home Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, results_welfare_H_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_H = zlim;
plot3([final_results.opt_phi_pi_H_H, final_results.opt_phi_pi_H_H], ...
      [final_results.opt_phi_y_H_H, final_results.opt_phi_y_H_H], ...
      [z_limits_H(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{\pi}^H'); ylabel('\phi_{y}^H'); zlabel('Welfare');
line1_H = sprintf('%s: Optimal for Home (%s shock)', plot_title_prefix, shock_scenario);
line2_H = sprintf('\\phi_{\\pi}^H=%.12f, \\phi_{y}^H=%.12f, Welfare=%.4f', ...
                  final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_home, fullfile(base_res_dir, 'graph_home.pdf'), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_pi_H_range_final; data_y = phi_y_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H];
save(fullfile(base_res_dir, 'data_for_plotting_home.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'TR CPI: Foreign Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, results_welfare_F_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_F = zlim;
plot3([final_results.opt_phi_pi_H_F, final_results.opt_phi_pi_H_F], ...
      [final_results.opt_phi_y_H_F, final_results.opt_phi_y_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{\pi}^H'); ylabel('\phi_{y}^H'); zlabel('Welfare');
line1_F = sprintf('%s: Optimal for Foreign (%s shock)', plot_title_prefix, shock_scenario);
line2_F = sprintf('\\phi_{\\pi}^H=%.12f, \\phi_{y}^H=%.12f, Welfare=%.4f', ...
                  final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_foreign, fullfile(base_res_dir, 'graph_foreign.pdf'), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F];
save(fullfile(base_res_dir, 'data_for_plotting_foreign.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global (Sum) Welfare ---
h_global = figure('Name', 'TR CPI: Global Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, results_welfare_SUM_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_G = zlim;
plot3([final_results.opt_phi_pi_H_SUM, final_results.opt_phi_pi_H_SUM], ...
      [final_results.opt_phi_y_H_SUM, final_results.opt_phi_y_H_SUM], ...
      [z_limits_G(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{\pi}^H'); ylabel('\phi_{y}^H'); zlabel('Welfare');
line1_G = sprintf('%s: Optimal for Global (Sum) (%s shock)', plot_title_prefix, shock_scenario);
line2_G = sprintf('\\phi_{\\pi}^H=%.12f, \\phi_{y}^H=%.12f, Welfare=%.4f', ...
                  final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_global, fullfile(base_res_dir, 'graph_global.pdf'), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM];
save(fullfile(base_res_dir, 'data_for_plotting_global.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('TR CPI Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');