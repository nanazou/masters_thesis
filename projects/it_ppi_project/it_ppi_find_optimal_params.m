% it_ppi_find_optimal_params.m
% このスクリプトは、IT_PPI政策において需要ショック(beta)と供給ショック(a)の
% 両方を考慮した「期待厚生 (Expected Welfare)」を最大化する共通パラメータ(phi_pi_H)を探索し、
% 厚生の2次元線グラフを出力します。計算の核心部は src/it_ppi_solve_core.m を呼び出します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
cd(mpath); % スクリプトの場所（Root）をカレントディレクトリに固定
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: it_ppi_project -> it_ppi)
[~, proj_folder_name] = fileparts(mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_;

fprintf('====================================================\n');
fprintf('Starting Expected Welfare optimization for IT_PPI Policy...\n');
fprintf('====================================================\n\n');

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- 設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼

% --- 1. 両ショックのパラメータと重み設定 ---
beta_shock_name  = 'eps_beta_H';
beta_shock_value = 0.02;   % 需要ショック (割引因子上昇)
a_shock_name     = 'eps_a_H';
a_shock_value    = -0.1;   % 供給ショック (生産性低下)

weight_beta = 0.5;         % 期待厚生の重み (betaショック)
weight_a    = 0.5;         % 期待厚生の重み (aショック)

% --- 2. 探索範囲の設定 (1次元ラインサーチ: 10点) ---
phi_pi_H_range = [3000, 5000, 7000, 9000, 11000, 13000, 14000, 16000, 18000, 20000];

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

if isempty(phi_pi_H_range)
    error("The parameter range for the search is empty. Check your settings.");
end

% --- 出力フォルダの定義と作成 (サブフォルダを作らず optimal_params 直下に保存) ---
base_res_dir = fullfile(mpath, 'results', 'optimal_params');
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- Simulation Settings Structures (Beta & A Shocks) ---
settings_base.modnam = 'it_ppi_model';
settings_base.modnamstar = 'it_ppi_model_zlb';
settings_base.constraint = 'i_H_notional < -i_H_ss';
settings_base.constraint_relax = 'i_H_notional > -i_H_ss';
settings_base.nperiods = 400;
settings_base.maxiter = 100;

% Betaショック用設定
settings_beta = settings_base;
settings_beta.irfshock = beta_shock_name;
settings_beta.shockssequence = beta_shock_value;

% Aショック用設定
settings_a = settings_base;
settings_a.irfshock = a_shock_name;
settings_a.shockssequence = a_shock_value;

% =========================================================================
%% --- Main Calculation Loop (1D Line Search for Both Beta and A Shocks) ---
% =========================================================================
fprintf('--- Starting 1D Line search for phi_pi_H (Expected Welfare) ---\n');

total_iterations = length(phi_pi_H_range);

% 期待厚生 (Expected Welfare = 0.5 * W_beta + 0.5 * W_a)
results_welfare_H   = zeros(total_iterations, 1);
results_welfare_F   = zeros(total_iterations, 1);
results_welfare_SUM = zeros(total_iterations, 1);

% 内訳保存用 (Betaショック単独・Aショック単独・単純和)
results_welfare_H_beta   = zeros(total_iterations, 1);
results_welfare_F_beta   = zeros(total_iterations, 1);
results_welfare_SUM_beta = zeros(total_iterations, 1);

results_welfare_H_a      = zeros(total_iterations, 1);
results_welfare_F_a      = zeros(total_iterations, 1);
results_welfare_SUM_a    = zeros(total_iterations, 1);

results_welfare_H_total  = zeros(total_iterations, 1);

root_dir = pwd;

for i = 1:total_iterations
    phi_pi_H_value = phi_pi_H_range(i);

    fprintf('\nTrial (%d/%d): phi_pi_H = %.4g\n', i, total_iterations, phi_pi_H_value);

    try
        % -------------------------------------------------------------
        % [Step 1] Beta Shock の計算と3段階フィルター検査
        % -------------------------------------------------------------
        [res_beta, irf_beta, Mbase_b, oobase_b] = it_ppi_solve_core(phi_pi_H_value, settings_beta);

        % 1. 固有値判定: lambda > 0 (正の安定固有値のみ許可)
        stable_eigs_b = oobase_b.dr.eigval(abs(oobase_b.dr.eigval) < 1 - 1e-5);
        if any(real(stable_eigs_b) < -1e-5)
            error('[Beta Shock] Negative stable eigenvalue (lambda < 0) detected.');
        end

        % 2. OccBinの未収束（サイクリング打ち切り）の厳密な検出
        i_H_ss_b = Mbase_b.params(strcmp('i_H_ss', Mbase_b.param_names));
        zlb_actual_b   = (irf_beta.i_H_pie <= -i_H_ss_b + 1e-6);
        zlb_notional_b = (irf_beta.i_H_notional_pie < -i_H_ss_b - 1e-6);
        if any(zlb_actual_b ~= zlb_notional_b)
            error('[Beta Shock] OccBin did not converge within maxiter (regime cycling).');
        end

        % 3. Non-ZLB金利の2周期チャタリング（振動）の検出
        di_lin_b = diff(irf_beta.i_H_lin);
        zigzag_b = sum((di_lin_b(1:end-1) .* di_lin_b(2:end) < 0) & (abs(di_lin_b(1:end-1)) > 1e-6));
        if zigzag_b > 10
            error('[Beta Shock] Oscillating solution (zigzag chattering) detected.');
        end

        % -------------------------------------------------------------
        % [Step 2] A Shock の計算と3段階フィルター検査
        % -------------------------------------------------------------
        [res_a, irf_a, Mbase_a, oobase_a] = it_ppi_solve_core(phi_pi_H_value, settings_a);

        % 1. 固有値判定: lambda > 0
        stable_eigs_a = oobase_a.dr.eigval(abs(oobase_a.dr.eigval) < 1 - 1e-5);
        if any(real(stable_eigs_a) < -1e-5)
            error('[A Shock] Negative stable eigenvalue (lambda < 0) detected.');
        end

        % 2. OccBinの未収束の厳密な検出
        i_H_ss_a = Mbase_a.params(strcmp('i_H_ss', Mbase_a.param_names));
        zlb_actual_a   = (irf_a.i_H_pie <= -i_H_ss_a + 1e-6);
        zlb_notional_a = (irf_a.i_H_notional_pie < -i_H_ss_a - 1e-6);
        if any(zlb_actual_a ~= zlb_notional_a)
            error('[A Shock] OccBin did not converge within maxiter (regime cycling).');
        end

        % 3. Non-ZLB金利の2周期チャタリング（振動）の検出
        di_lin_a = diff(irf_a.i_H_lin);
        zigzag_a = sum((di_lin_a(1:end-1) .* di_lin_a(2:end) < 0) & (abs(di_lin_a(1:end-1)) > 1e-6));
        if zigzag_a > 10
            error('[A Shock] Oscillating solution (zigzag chattering) detected.');
        end

        % -------------------------------------------------------------
        % [Step 3] 両ショック通過時のみ各厚生と期待厚生を抽出・計算
        % -------------------------------------------------------------
        w_H_b   = res_beta.welfare_H_with_delta_pie;
        w_F_b   = res_beta.welfare_F_with_delta_pie;
        w_SUM_b = res_beta.welfare_SUM_with_delta_pie;

        w_H_a   = res_a.welfare_H_with_delta_pie;
        w_F_a   = res_a.welfare_F_with_delta_pie;
        w_SUM_a = res_a.welfare_SUM_with_delta_pie;

        welfare_H     = weight_beta * w_H_b   + weight_a * w_H_a;
        welfare_F     = weight_beta * w_F_b   + weight_a * w_F_a;
        welfare_SUM   = weight_beta * w_SUM_b + weight_a * w_SUM_a;
        welfare_H_tot = w_H_b + w_H_a;

    catch ME
        fprintf('  -> Excluded: %s\n', ME.message);
        cd(root_dir);
        w_H_b = -999; w_F_b = -999; w_SUM_b = -999;
        w_H_a = -999; w_F_a = -999; w_SUM_a = -999;
        welfare_H = -999; welfare_F = -999; welfare_SUM = -999;
        welfare_H_tot = -999;
    end

    % 配列への格納
    results_welfare_H_beta(i)   = w_H_b;
    results_welfare_F_beta(i)   = w_F_b;
    results_welfare_SUM_beta(i) = w_SUM_b;

    results_welfare_H_a(i)      = w_H_a;
    results_welfare_F_a(i)      = w_F_a;
    results_welfare_SUM_a(i)    = w_SUM_a;

    results_welfare_H(i)        = welfare_H;
    results_welfare_F(i)        = welfare_F;
    results_welfare_SUM(i)      = welfare_SUM;
    results_welfare_H_total(i)  = welfare_H_tot;

    if welfare_SUM > -999
        fprintf('  -> [OK] Home (Beta: %.6f, A: %.6f) => Expected H: %.6f (Total: %.6f)\n', ...
            w_H_b, w_H_a, welfare_H, welfare_H_tot);
    end
end

if max(results_welfare_SUM(:)) <= -999
    fprintf('--- SEARCH FAILED: All parameter combinations resulted in an error. ---\n');
    return;
end

% =========================================================================
%% --- Find and Package Optimal Parameters ---
% =========================================================================
% 1. 自国（Home）の期待厚生を最大化するパラメータ
[max_welfare_H, idx_H] = max(results_welfare_H);
final_results.opt_phi_pi_H_H          = phi_pi_H_range(idx_H);
final_results.max_welfare_H           = max_welfare_H;
final_results.max_welfare_H_total     = results_welfare_H_total(idx_H);
final_results.welfare_H_beta_at_H_opt = results_welfare_H_beta(idx_H);
final_results.welfare_H_a_at_H_opt    = results_welfare_H_a(idx_H);
final_results.welfare_F_at_H_opt      = results_welfare_F(idx_H);
final_results.welfare_SUM_at_H_opt    = results_welfare_SUM(idx_H);

% 2. 外国（Foreign）の期待厚生を最大化するパラメータ
[max_welfare_F, idx_F] = max(results_welfare_F);
final_results.opt_phi_pi_H_F          = phi_pi_H_range(idx_F);
final_results.max_welfare_F           = max_welfare_F;
final_results.welfare_F_beta_at_F_opt = results_welfare_F_beta(idx_F);
final_results.welfare_F_a_at_F_opt    = results_welfare_F_a(idx_F);
final_results.welfare_H_at_F_opt      = results_welfare_H(idx_F);
final_results.welfare_SUM_at_F_opt    = results_welfare_SUM(idx_F);

% 3. 全世界（Global Sum）の期待厚生を最大化するパラメータ
[max_welfare_SUM, idx_SUM] = max(results_welfare_SUM);
final_results.opt_phi_pi_H_SUM          = phi_pi_H_range(idx_SUM);
final_results.max_welfare_SUM           = max_welfare_SUM;
final_results.welfare_SUM_beta_at_S_opt = results_welfare_SUM_beta(idx_SUM);
final_results.welfare_SUM_a_at_S_opt    = results_welfare_SUM_a(idx_SUM);
final_results.welfare_H_at_SUM_opt      = results_welfare_H(idx_SUM);
final_results.welfare_F_at_SUM_opt      = results_welfare_F(idx_SUM);

% 全データの保存 (例: it_ppi_opt_all_data.mat)
all_data_filename = sprintf('%s_opt_all_data.mat', policy_key);
save(fullfile(base_res_dir, all_data_filename), ...
    'final_results', 'phi_pi_H_range', ...
    'results_welfare_H', 'results_welfare_F', 'results_welfare_SUM', ...
    'results_welfare_H_beta', 'results_welfare_F_beta', 'results_welfare_SUM_beta', ...
    'results_welfare_H_a', 'results_welfare_F_a', 'results_welfare_SUM_a', ...
    'results_welfare_H_total');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n==============================================================\n');
fprintf('--- IT_PPI Final Optimal Policy Results (Expected Welfare) ---\n');
fprintf('==============================================================\n');
fprintf('1. Optimal for Home:\n');
fprintf('   Parameter: phi_pi_H = %.4g\n', final_results.opt_phi_pi_H_H);
fprintf('   Expected Welfare H (0.5*Beta + 0.5*A) = %.6f (Total = %.6f)\n', ...
    final_results.max_welfare_H, final_results.max_welfare_H_total);
fprintf('     - Breakdown: Beta Shock H = %.6f, A Shock H = %.6f\n', ...
    final_results.welfare_H_beta_at_H_opt, final_results.welfare_H_a_at_H_opt);
fprintf('   Expected Foreign = %.6f, Expected Global Sum = %.6f\n\n', ...
    final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);

fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameter: phi_pi_H = %.4g\n', final_results.opt_phi_pi_H_F);
fprintf('   Resulting Expected Welfare: Home = %.6f, Foreign = %.6f, Sum = %.6f\n\n', ...
    final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);

fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameter: phi_pi_H = %.4g\n', final_results.opt_phi_pi_H_SUM);
fprintf('   Resulting Expected Welfare: Home = %.6f, Foreign = %.6f, Sum = %.6f\n', ...
    final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('==============================================================\n\n');

% =========================================================================
%% --- Plot and Save Graphs & Raw Data ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

% タイトルとパラメータ行の配置座標 (グラフ左端 x = 0.08 に揃えて 'left' で配置)
title_pos_line1 = [0.08, 0.945, 0.88, 0.045];
title_pos_line2 = [0.08, 0.895, 0.88, 0.045];

% 軸の端での数字の衝突を防ぐため、X軸の表示範囲の両端に8%の余白（パディング）を計算
x_min = min(phi_pi_H_range); x_max = max(phi_pi_H_range);
dx = max(x_max - x_min, 1e-12);
x_lim_pad = [x_min - 0.08 * dx, x_max + 0.08 * dx];
x_mid = mean(x_lim_pad);

% プロット用にエラー値 (-999) を NaN に置換した配列を用意 (Y軸の潰れ防止)
plot_welfare_H   = results_welfare_H;   plot_welfare_H(plot_welfare_H <= -999)     = NaN;
plot_welfare_F   = results_welfare_F;   plot_welfare_F(plot_welfare_F <= -999)     = NaN;
plot_welfare_SUM = results_welfare_SUM; plot_welfare_SUM(plot_welfare_SUM <= -999) = NaN;

% --- Helper function for 2D axis styling ---
function style_2d_axes(ax, x_lim_pad, x_mid, y_limits, dy)
    set(ax, 'Position', [0.16, 0.18, 0.74, 0.62], 'FontSize', 8);
    % Y軸ラベル (Welfare) をY軸の頂点上に配置
    text(x_lim_pad(1), y_limits(2) + 0.06 * dy, 'Welfare', ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
    % X軸ラベル (\phi_{\pi}^H) を改行なしでX軸下中央に明示配置 (gnuplotで消えるのを防止)
    text(x_mid, y_limits(1) - 0.14 * dy, '\phi_{\pi}^{H}', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Interpreter', 'tex');
end

% --- Graph 1: Home Expected Welfare ---
h_home = figure('Name', 'IT_PPI: Home Expected Welfare', 'Visible', 'off');
plot(phi_pi_H_range, plot_welfare_H, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
xlim(x_lim_pad);
y_limits_H = ylim;
dy_H = max(y_limits_H(2) - y_limits_H(1), 1e-6);
ylim([y_limits_H(1) - 0.08 * dy_H, y_limits_H(2) + 0.08 * dy_H]);
y_limits_H = ylim;
plot(final_results.opt_phi_pi_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
style_2d_axes(gca, x_lim_pad, x_mid, y_limits_H, dy_H);
hold off;
grid on;
line1_H = 'IT_{PPI}: Home Expected Welfare';
line2_H = sprintf('\\phi_{\\pi}^{H} = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_pi_H_H, final_results.max_welfare_H, ...
    final_results.welfare_H_beta_at_H_opt, final_results.welfare_H_a_at_H_opt);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_home_filename = sprintf('%s_opt_graph_home.pdf', policy_key);
print(h_home, fullfile(base_res_dir, graph_home_filename), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_pi_H_range; data_y = results_welfare_H;
opt_pt = [final_results.opt_phi_pi_H_H, final_results.max_welfare_H];
mat_home_filename = sprintf('%s_opt_data_for_plotting_home.mat', policy_key);
save(fullfile(base_res_dir, mat_home_filename), 'data_x', 'data_y', 'opt_pt');

% --- Graph 2: Foreign Expected Welfare ---
h_foreign = figure('Name', 'IT_PPI: Foreign Expected Welfare', 'Visible', 'off');
plot(phi_pi_H_range, plot_welfare_F, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
xlim(x_lim_pad);
y_limits_F = ylim;
dy_F = max(y_limits_F(2) - y_limits_F(1), 1e-6);
ylim([y_limits_F(1) - 0.08 * dy_F, y_limits_F(2) + 0.08 * dy_F]);
y_limits_F = ylim;
plot(final_results.opt_phi_pi_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
style_2d_axes(gca, x_lim_pad, x_mid, y_limits_F, dy_F);
hold off;
grid on;
line1_F = 'IT_{PPI}: Foreign Expected Welfare';
line2_F = sprintf('\\phi_{\\pi}^{H} = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_pi_H_F, final_results.max_welfare_F, ...
    final_results.welfare_F_beta_at_F_opt, final_results.welfare_F_a_at_F_opt);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_foreign_filename = sprintf('%s_opt_graph_foreign.pdf', policy_key);
print(h_foreign, fullfile(base_res_dir, graph_foreign_filename), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_x = phi_pi_H_range; data_y = results_welfare_F;
opt_pt = [final_results.opt_phi_pi_H_F, final_results.max_welfare_F];
mat_foreign_filename = sprintf('%s_opt_data_for_plotting_foreign.mat', policy_key);
save(fullfile(base_res_dir, mat_foreign_filename), 'data_x', 'data_y', 'opt_pt');

% --- Graph 3: Global (Sum) Expected Welfare ---
h_global = figure('Name', 'IT_PPI: Global Expected Welfare', 'Visible', 'off');
plot(phi_pi_H_range, plot_welfare_SUM, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
xlim(x_lim_pad);
y_limits_G = ylim;
dy_G = max(y_limits_G(2) - y_limits_G(1), 1e-6);
ylim([y_limits_G(1) - 0.08 * dy_G, y_limits_G(2) + 0.08 * dy_G]);
y_limits_G = ylim;
plot(final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
style_2d_axes(gca, x_lim_pad, x_mid, y_limits_G, dy_G);
hold off;
grid on;
line1_G = 'IT_{PPI}: Global Expected Welfare';
line2_G = sprintf('\\phi_{\\pi}^{H} = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM, ...
    final_results.welfare_SUM_beta_at_S_opt, final_results.welfare_SUM_a_at_S_opt);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_global_filename = sprintf('%s_opt_graph_global.pdf', policy_key);
print(h_global, fullfile(base_res_dir, graph_global_filename), '-dpdf');
close(h_global);

% Global個別データ保存
data_x = phi_pi_H_range; data_y = results_welfare_SUM;
opt_pt = [final_results.opt_phi_pi_H_SUM, final_results.max_welfare_SUM];
mat_global_filename = sprintf('%s_opt_data_for_plotting_global.mat', policy_key);
save(fullfile(base_res_dir, mat_global_filename), 'data_x', 'data_y', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('IT_PPI Policy expected welfare optimization has finished.\n');
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');