% polt_find_optimal_params.m
% このスクリプトは、POLT政策の最適パラメータを探索します。
% Stage 1 で潜在産出量パスを算出し、Stage 2 の計算エンジンとして src/stage2/polt_stage2_solve_core.m を呼び出します。

clear;
close all;
clc;

% スクリプトがある場所を特定してカレントフォルダを移動
mpath = fileparts(mfilename('fullpath'));
cd(mpath);

% --- パスの設定 ---
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_ options_;

fprintf('====================================================\n');
fprintf('Starting optimization for POLT Policy...\n');
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
shock_scenario = 'beta'; 

% --- 3. 粗い探索 (coarse mode) の設定 --- a final
%phi_gap_H_range_coarse   = 0.005:0.0001:0.0052;
%phi_level_H_range_coarse = 0:1:1;

% --- 4. 詳細な探索 (fine mode) の設定 --- beta final
%center_phi_gap_H   = 0.3; 
%center_phi_level_H = 0.06; 
%fine_search_radius_gap_H   = 0.02; 
%fine_search_step_gap_H     = 0.01; 
%fine_search_radius_level_H = 0.01; 
%fine_search_step_level_H   = 0.005;

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

%% --- [1] ショック設定およびフォルダ作成 ---
if strcmp(shock_scenario, 'beta')
    target_shock = 'eps_beta_H';
    shock_value  = 0.02;
    shock_folder = 'beta_shock';
elseif strcmp(shock_scenario, 'a')
    target_shock = 'eps_a_H';
    shock_value  = -0.1;
    shock_folder = 'a_shock';
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% 出力フォルダの定義: results/optimal_params/[shock_folder]/ 直下にすべて出力
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- [2] POLT Stage 1: 潜在産出量パスの自動計算 ---
fprintf('\n--- POLT Stage 1: Calculating target path for %s shock ---\n', shock_scenario);

% Stage 1 フォルダへ移動して実行
root_dir = pwd;
cd(fullfile(mpath, 'src', 'stage1'));

% パラメータダミー
phi_gap_H_value = 0; phi_level_H_value = 0;
assignin('base', 'phi_gap_H_value', phi_gap_H_value);
assignin('base', 'phi_level_H_value', phi_level_H_value);

% Stage 1 モデルの実行 (以前定義した polt_stage1_model を参照)
eval(['dynare polt_stage1_model.mod noclearall nolog']);

% ショックパスの計算 (simult_ を使用)
ex_flex = zeros(400, M_.exo_nbr);
s_idx = strmatch(target_shock, M_.exo_names, 'exact');
ex_flex(1, s_idx) = shock_value;
y_flex_sim = simult_(M_, options_, oo_.dr.ys, oo_.dr, ex_flex, 1);
y_H_idx = strmatch('y_H', M_.endo_names, 'exact');
eps_y_H_potential_vec = y_flex_sim(y_H_idx, 2:end)' - oo_.dr.ys(y_H_idx);

% ルートに戻る
cd(root_dir);

% Stage 1 のパスを保存
save(fullfile(base_res_dir, 'stage1_potential_path.mat'), 'eps_y_H_potential_vec');

fprintf('Stage 1 完了。Stage 2 最適化を開始します...\n');

%% --- [3] Stage 2 (OccBin) 用の共通設定 ---
% Stage 2 のエンジン（src/stage2/polt_stage2_solve_core.m）をパスに追加
addpath(fullfile(mpath, 'src', 'stage2'));

settings.modnam = 'polt_stage2_model';
settings.modnamstar = 'polt_stage2_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';

% 構造ショックと潜在産出量ショックの2つをセット
settings.irfshock = char(target_shock, 'eps_y_H_potential');
settings.shockssequence = zeros(400, 2);
settings.shockssequence(1, 1) = shock_value;
settings.shockssequence(:, 2) = eps_y_H_potential_vec;
settings.nperiods = 400;
settings.maxiter = 100;

% =========================================================================
%% --- [4] グリッド設定とループ実行 ---
% =========================================================================
if strcmp(search_mode, 'coarse')
    fprintf('--- Starting COARSE grid search ---\n');
    phi_gap_H_range_final   = phi_gap_H_range_coarse;
    phi_level_H_range_final = phi_level_H_range_coarse;
    plot_title_prefix = 'POLT Coarse Search';
elseif strcmp(search_mode, 'fine')
    fprintf('--- Starting FINE grid search around (gap_H=%.2f, level_H=%.2f) ---\n', center_phi_gap_H, center_phi_level_H);
    phi_gap_H_range_final   = (center_phi_gap_H - fine_search_radius_gap_H):fine_search_step_gap_H:(center_phi_gap_H + fine_search_radius_gap_H);
    phi_level_H_range_final = (center_phi_level_H - fine_search_radius_level_H):fine_search_step_level_H:(center_phi_level_H + fine_search_radius_level_H);
    phi_gap_H_range_final(phi_gap_H_range_final <= 0) = [];
    phi_level_H_range_final(phi_level_H_range_final <= 0) = [];
    plot_title_prefix = 'POLT Fine Search';
else
    error("Invalid search_mode. Choose 'coarse' or 'fine'.");
end

if isempty(phi_gap_H_range_final) || isempty(phi_level_H_range_final)
    error("The parameter range for the search is empty. Check your settings.");
end

results_welfare_H_final   = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));
results_welfare_F_final   = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));
results_welfare_SUM_final = zeros(length(phi_level_H_range_final), length(phi_gap_H_range_final));

total_iterations = length(phi_level_H_range_final) * length(phi_gap_H_range_final);
current_iter = 0;

for i = 1:length(phi_level_H_range_final)
    for j = 1:length(phi_gap_H_range_final)
        current_iter = current_iter + 1;
        phi_gap_H_val = phi_gap_H_range_final(j);
        phi_level_H_val = phi_level_H_range_final(i);
        
        fprintf('\nTrial (%d/%d): phi_gap_H = %.2f, phi_level_H = %.2f ... ', ...
            current_iter, total_iterations, phi_gap_H_val, phi_level_H_val);
        
        try
            % src/stage2/polt_stage2_solve_core.m を呼び出す
            res = polt_stage2_solve_core(phi_gap_H_val, phi_level_H_val, settings);
            
            welfare_H   = res.welfare_H_with_delta_pie;
            welfare_F   = res.welfare_F_with_delta_pie;
            welfare_SUM = res.welfare_SUM_with_delta_pie;
        catch ME
            fprintf('Error during trial: %s\n', ME.message);
            cd(mpath);
            welfare_H = -999; welfare_F = -999; welfare_SUM = -999;
        end
        
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
%% --- [5] 最適パラメータの特定とパッケージング ---
% =========================================================================
% Home厚生最大
[max_w_H, idx_H] = max(results_welfare_H_final(:));
[row_H, col_H] = ind2sub(size(results_welfare_H_final), idx_H);
final_results.opt_phi_gap_H_H   = phi_gap_H_range_final(col_H);
final_results.opt_phi_level_H_H = phi_level_H_range_final(row_H);
final_results.max_welfare_H = max_w_H;
final_results.welfare_F_at_H_opt = results_welfare_F_final(row_H, col_H);
final_results.welfare_SUM_at_H_opt = results_welfare_SUM_final(row_H, col_H);

% Foreign厚生最大
[max_w_F, idx_F] = max(results_welfare_F_final(:));
[row_F, col_F] = ind2sub(size(results_welfare_F_final), idx_F);
final_results.opt_phi_gap_H_F   = phi_gap_H_range_final(col_F);
final_results.opt_phi_level_H_F = phi_level_H_range_final(row_F);
final_results.max_welfare_F = max_w_F;

% Global厚生最大
[max_w_SUM, idx_SUM] = max(results_welfare_SUM_final(:));
[row_SUM, col_SUM] = ind2sub(size(results_welfare_SUM_final), idx_SUM);
final_results.opt_phi_gap_H_SUM   = phi_gap_H_range_final(col_SUM);
final_results.opt_phi_level_H_SUM = phi_level_H_range_final(row_SUM);
final_results.max_welfare_SUM = max_w_SUM;
final_results.welfare_H_at_SUM_opt = results_welfare_H_final(row_SUM, col_SUM);
final_results.welfare_F_at_SUM_opt = results_welfare_F_final(row_SUM, col_SUM);

% 全データの保存 (all_data.mat)
save(fullfile(base_res_dir, 'all_data.mat'), ...
    'final_results', 'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'phi_gap_H_range_final', 'phi_level_H_range_final');

% =========================================================================
%% --- [6] ログ出力 ---
% =========================================================================
fprintf('\n--- POLT Final Optimal Policy Results (Shock: %s) ---\n', shock_scenario);
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_gap_H = %.2f, phi_level_H = %.2f\n', final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_gap_H = %.2f, phi_level_H = %.2f\n', final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n', final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('--------------------------------------------------------------\n\n');

% =========================================================================
%% --- [7] グラフ描画と保存 ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

title_pos_line1 = [0, 0.96, 1, 0.05];
title_pos_line2 = [0, 0.91, 1, 0.05];

% --- Graph 1: Home Welfare ---
h_home = figure('Name', 'POLT: Home Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_H_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits = zlim;
plot3([final_results.opt_phi_gap_H_H, final_results.opt_phi_gap_H_H], ...
      [final_results.opt_phi_level_H_H, final_results.opt_phi_level_H_H], ...
      [z_limits(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
line1 = sprintf('%s: Optimal for Home (%s)', plot_title_prefix, shock_scenario);
line2 = sprintf('\\phi_{gap}^H=%.10f, \\phi_{level}^H=%.10f, Welfare=%.4f', final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_home, fullfile(base_res_dir, 'graph_home.pdf'), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_gap_H_range_final; data_y = phi_level_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H];
save(fullfile(base_res_dir, 'data_for_plotting_home.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'POLT: Foreign Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_F_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits_F = zlim;
plot3([final_results.opt_phi_gap_H_F, final_results.opt_phi_gap_H_F], ...
      [final_results.opt_phi_level_H_F, final_results.opt_phi_level_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
print(h_foreign, fullfile(base_res_dir, 'graph_foreign.pdf'), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F];
save(fullfile(base_res_dir, 'data_for_plotting_foreign.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global Welfare ---
h_global = figure('Name', 'POLT: Global Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, results_welfare_SUM_final);
colorbar; view(3); hold on;
plot3(final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
z_limits = zlim;
plot3([final_results.opt_phi_gap_H_SUM, final_results.opt_phi_gap_H_SUM], ...
      [final_results.opt_phi_level_H_SUM, final_results.opt_phi_level_H_SUM], ...
      [z_limits(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
hold off;
xlabel('\phi_{gap}^H'); ylabel('\phi_{level}^H'); zlabel('Welfare');
line1 = sprintf('%s: Optimal for Global (%s)', plot_title_prefix, shock_scenario);
line2 = sprintf('\\phi_{gap}^H=%.10f, \\phi_{level}^H=%.10f, Welfare=%.4f', final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
print(h_global, fullfile(base_res_dir, 'graph_global.pdf'), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM];
save(fullfile(base_res_dir, 'data_for_plotting_global.mat'), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('POLT Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');