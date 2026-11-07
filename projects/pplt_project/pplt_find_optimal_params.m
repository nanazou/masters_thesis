% pplt_find_optimal_params.m
% このスクリプトは、PPLT政策の最適パラメータを探索し、
% 厚生の表面グラフを出力します。計算の核心部は src/pplt_solve_core.m を呼び出します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
cd(mpath); % スクリプトの場所（Root）をカレントディレクトリに固定
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: pplt_project -> pplt)
[~, proj_folder_name] = fileparts(mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_;

fprintf('====================================================\n');
fprintf('Starting optimization for PPLT Policy...\n');
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
shock_scenario = 'a'; % <--- ここを 'beta' または 'a' に書き換えて切り替えます

% --- フォルダ用ショック名の決定 ---
if strcmp(shock_scenario, 'beta')
    shock_folder = 'beta_shock';
else
    shock_folder = 'a_shock';
end

% --- 3. ショックパラメータの設定 (run_pplt.m の設定と同期) ---
if strcmp(shock_scenario, 'beta')
    target_shock = 'eps_beta_H';
    shock_value  = 0.02; 
elseif strcmp(shock_scenario, 'a')
    target_shock = 'eps_a_H';
    shock_value  = -0.1; 
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% --- 4. 粗い探索 (coarse mode) の設定 --- a shock (0,0)近傍の探索 (5 × 4 = 20点)
phi_gap_H_range_coarse   = [0, 5e-6, 1e-5, 1.5e-5, 2e-5];
phi_level_H_range_coarse = [0, 2e-8, 5e-8, 1e-7];

% --- 4. 粗い探索 (coarse mode) の設定 --- 旧 a final
%phi_gap_H_range_coarse   = 0.000000000001:0.000000000001:0.000000000005;
%phi_level_H_range_coarse = 0.000000000001:0.000000000001:0.000000000004;

% --- 4. 粗い探索 (coarse mode) の設定 --- beta (75〜95のピーク確認: 3×3 = 全9点)
%phi_gap_H_range_coarse   = [75, 85, 95];
%phi_level_H_range_coarse = [1300, 1500, 1700];

% --- 4. 粗い探索 (coarse mode) の設定 --- beta final
%phi_gap_H_range_coarse   = 90:10:120;
%phi_level_H_range_coarse = 490:10:510;

% --- 5. 詳細な探索 (fine mode) の設定 ---
center_phi_gap_H   = 80.5; 
center_phi_level_H = 33.1; 
fine_search_radius_gap_H   = 2; 
fine_search_step_gap_H     = 0.1; 
fine_search_radius_level_H = 2; 
fine_search_step_level_H   = 0.1;

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

% --- 出力フォルダの定義と作成 ---
% 構成: results/optimal_params/[shock_folder]/ 直下にすべて出力
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- Simulation Settings Structure ---
% 独自名に変更した統一ファイル名 ('pplt_model', 'pplt_model_zlb') を参照
settings.modnam = 'pplt_model';
settings.modnamstar = 'pplt_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.irfshock = target_shock; 
settings.shockssequence = shock_value; 
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
    fprintf('--- Starting FINE grid search around (gap_H=%.4g, level_H=%.4g) ---\n', center_phi_gap_H, center_phi_level_H);
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

        fprintf('\nTrial (%d/%d): phi_gap_H = %.4g, phi_level_H = %.4g ... ', ...
            current_iteration, total_iterations, phi_gap_H_value, phi_level_H_value);

        % --- [計算核心 pplt_solve_core を呼び出す] ---
        try
            % 共通エンジン pplt_solve_core.m を呼び出し (4つの戻り値を取得)
            [res, irf_res, Mbase_, oobase_] = pplt_solve_core(phi_gap_H_value, phi_level_H_value, settings);

            % ▼▼▼【固有値判定: lambda > 0 (正の安定固有値のみ許可)】▼▼▼
            % 安定固有値（絶対値 < 1）のなかに負の実部を持つもの（lambda < 0）があるか検査
            stable_eigs = oobase_.dr.eigval(abs(oobase_.dr.eigval) < 1 - 1e-5);
            if any(real(stable_eigs) < -1e-5)
                error('Negative stable eigenvalue (lambda < 0) detected.');
            end

            % ▼▼▼【OccBinの未収束（サイクリング打ち切り）の厳密な検出】▼▼▼
            % 収束していれば、ZLBバインド期(i_H_pie <= -i_H_ss + 1e-6)と
            % i_H_notional_pie < -i_H_ss の期間が完全に一致する。
            % maxiterに達して未収束で終わった点はこれが不一致になるため除外する。
            i_H_ss_val = Mbase_.params(strcmp('i_H_ss', Mbase_.param_names));
            zlb_actual   = (irf_res.i_H_pie <= -i_H_ss_val + 1e-6);
            zlb_notional = (irf_res.i_H_notional_pie < -i_H_ss_val - 1e-6);
            if any(zlb_actual ~= zlb_notional)
                error('OccBin did not converge within maxiter (regime cycling).');
            end

            % フェイルセーフ：Non-ZLB金利の2周期チャタリング（振動）の検出
            di_lin = diff(irf_res.i_H_lin);
            zigzag_count = sum((di_lin(1:end-1) .* di_lin(2:end) < 0) & (abs(di_lin(1:end-1)) > 1e-6));
            if zigzag_count > 10
                error('Oscillating solution (zigzag chattering) detected.');
            end
            % ▲▲▲【固有値判定ここまで】▲▲▲

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
[max_welfare_H, idx_H] = max(results_welfare_H_final(:));
[row_H, col_H] = ind2sub(size(results_welfare_H_final), idx_H);
final_results.opt_phi_gap_H_H   = phi_gap_H_range_final(col_H);
final_results.opt_phi_level_H_H = phi_level_H_range_final(row_H);
final_results.max_welfare_H = max_welfare_H;
final_results.welfare_F_at_H_opt = results_welfare_F_final(row_H, col_H);
final_results.welfare_SUM_at_H_opt = results_welfare_SUM_final(row_H, col_H);

% 2. Optimal for Foreign
[max_welfare_F, idx_F] = max(results_welfare_F_final(:));
[row_F, col_F] = ind2sub(size(results_welfare_F_final), idx_F);
final_results.opt_phi_gap_H_F   = phi_gap_H_range_final(col_F);
final_results.opt_phi_level_H_F = phi_level_H_range_final(row_F);
final_results.max_welfare_F = max_welfare_F;
final_results.welfare_H_at_F_opt = results_welfare_H_final(row_F, col_F);
final_results.welfare_SUM_at_F_opt = results_welfare_SUM_final(row_F, col_F);

% 3. Optimal for Global (Sum)
[max_welfare_SUM, idx_SUM] = max(results_welfare_SUM_final(:));
[row_SUM, col_SUM] = ind2sub(size(results_welfare_SUM_final), idx_SUM);
final_results.opt_phi_gap_H_SUM   = phi_gap_H_range_final(col_SUM);
final_results.opt_phi_level_H_SUM = phi_level_H_range_final(row_SUM);
final_results.max_welfare_SUM = max_welfare_SUM;
final_results.welfare_H_at_SUM_opt = results_welfare_H_final(row_SUM, col_SUM);
final_results.welfare_F_at_SUM_opt = results_welfare_F_final(row_SUM, col_SUM);

% 全データの保存 (例: pplt_opt_all_data_a_shock.mat)
all_data_filename = sprintf('%s_opt_all_data_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, all_data_filename), ...
    'final_results', 'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'phi_gap_H_range_final', 'phi_level_H_range_final');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n--- PPLT Final Optimal Policy Results (Shock: %s) ---\n', shock_scenario);
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);
fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n', final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('--------------------------------------------------------------\n\n');

% =========================================================================
%% --- Plot and Save Graphs & Raw Data ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

title_pos_line1 = [0, 0.94, 1, 0.05];
title_pos_line2 = [0, 0.89, 1, 0.05];

% 政策名（大文字）と探索モードの簡潔な表記
policy_name_upper = upper(policy_key);
mode_label = [upper(search_mode(1)), lower(search_mode(2:end))]; % 'Coarse' or 'Fine'

% 角での数字の衝突を防ぐため、X軸・Y軸の表示範囲の両端に8%の余白（パディング）を計算
x_min = min(phi_gap_H_range_final);   x_max = max(phi_gap_H_range_final);
y_min = min(phi_level_H_range_final); y_max = max(phi_level_H_range_final);
dx = max(x_max - x_min, 1e-15);
dy = max(y_max - y_min, 1e-15);
x_lim_pad = [x_min - 0.08 * dx, x_max + 0.08 * dx];
y_lim_pad = [y_min - 0.08 * dy, y_max + 0.08 * dy];

% プロット用にエラー値 (-999) を NaN に置換した配列を用意 (Z軸の潰れ防止)
plot_welfare_H   = results_welfare_H_final;   plot_welfare_H(plot_welfare_H <= -999)     = NaN;
plot_welfare_F   = results_welfare_F_final;   plot_welfare_F(plot_welfare_F <= -999)     = NaN;
plot_welfare_SUM = results_welfare_SUM_final; plot_welfare_SUM(plot_welfare_SUM <= -999) = NaN;

% --- Graph 1: Home Welfare ---
h_home = figure('Name', 'PPLT: Home Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_H);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_H = zlim;
dz_H = max(z_limits_H(2) - z_limits_H(1), 1e-6);
zlim([z_limits_H(1) - 0.08 * dz_H, z_limits_H(2) + 0.08 * dz_H]);
z_limits_H = zlim;
plot3(final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_H, final_results.opt_phi_gap_H_H], ...
      [final_results.opt_phi_level_H_H, final_results.opt_phi_level_H_H], ...
      [z_limits_H(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
% Z軸ラベル(Welfare)は数字との重なりを防ぐためZ軸の頂点上に横書きで配置
text(x_lim_pad(1), y_lim_pad(2), z_limits_H(2) + 0.08 * dz_H, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{gap}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{level}^H'), 'FontSize', 10);
line1_H = sprintf('%s: Home Welfare (%s shock, %s)', policy_name_upper, shock_scenario, mode_label);
line2_H = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_home_filename = sprintf('%s_opt_graph_home_%s.pdf', policy_key, shock_folder);
print(h_home, fullfile(base_res_dir, graph_home_filename), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_gap_H_range_final; data_y = phi_level_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H];
mat_home_filename = sprintf('%s_opt_data_for_plotting_home_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_home_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'PPLT: Foreign Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_F);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_F = zlim;
dz_F = max(z_limits_F(2) - z_limits_F(1), 1e-6);
zlim([z_limits_F(1) - 0.08 * dz_F, z_limits_F(2) + 0.08 * dz_F]);
z_limits_F = zlim;
plot3(final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_F, final_results.opt_phi_gap_H_F], ...
      [final_results.opt_phi_level_H_F, final_results.opt_phi_level_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
text(x_lim_pad(1), y_lim_pad(2), z_limits_F(2) + 0.08 * dz_F, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{gap}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{level}^H'), 'FontSize', 10);
line1_F = sprintf('%s: Foreign Welfare (%s shock, %s)', policy_name_upper, shock_scenario, mode_label);
line2_F = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_foreign_filename = sprintf('%s_opt_graph_foreign_%s.pdf', policy_key, shock_folder);
print(h_foreign, fullfile(base_res_dir, graph_foreign_filename), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F];
mat_foreign_filename = sprintf('%s_opt_data_for_plotting_foreign_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_foreign_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global (Sum) Welfare ---
h_global = figure('Name', 'PPLT: Global Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_SUM);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_G = zlim;
dz_G = max(z_limits_G(2) - z_limits_G(1), 1e-6);
zlim([z_limits_G(1) - 0.08 * dz_G, z_limits_G(2) + 0.08 * dz_G]);
z_limits_G = zlim;
plot3(final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_SUM, final_results.opt_phi_gap_H_SUM], ...
      [final_results.opt_phi_level_H_SUM, final_results.opt_phi_level_H_SUM], ...
      [z_limits_G(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
text(x_lim_pad(1), y_lim_pad(2), z_limits_G(2) + 0.08 * dz_G, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{gap}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{level}^H'), 'FontSize', 10);
line1_G = sprintf('%s: Global Welfare (%s shock, %s)', policy_name_upper, shock_scenario, mode_label);
line2_G = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_global_filename = sprintf('%s_opt_graph_global_%s.pdf', policy_key, shock_folder);
print(h_global, fullfile(base_res_dir, graph_global_filename), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM];
mat_global_filename = sprintf('%s_opt_data_for_plotting_global_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_global_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('PPLT Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');