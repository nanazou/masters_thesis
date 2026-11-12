% olt_find_optimal_params.m
% このスクリプトは、OLT政策において需要ショック(beta)と供給ショック(a)の
% 両方を考慮した「期待厚生 (Expected Welfare)」を最大化する共通パラメータを探索し、
% 厚生の表面グラフを出力します。計算の核心部は src/olt_solve_core.m を呼び出します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
cd(mpath); % スクリプトの場所（Root）をカレントディレクトリに固定
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: olt_project -> olt)
[~, proj_folder_name] = fileparts(mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_;

fprintf('====================================================\n');
fprintf('Starting Expected Welfare optimization for OLT Policy...\n');
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

% --- 2. 探索グリッドの設定 (4 x 4 = 16点) ---
phi_gap_H_range_final   = [1000000, 2000000, 3000000, 4000000];
phi_level_H_range_final = [1000000, 2000000, 3000000, 4000000];

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

if isempty(phi_gap_H_range_final) || isempty(phi_level_H_range_final)
    error("The parameter range for the search is empty. Check your settings.");
end

% --- 出力フォルダの定義と作成 (サブフォルダを作らず optimal_params 直下に保存) ---
base_res_dir = fullfile(mpath, 'results', 'optimal_params');
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- Simulation Settings Structures (Beta & A Shocks) ---
settings_base.modnam = 'olt_model';
settings_base.modnamstar = 'olt_model_zlb';
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
%% --- Main Calculation Loop (Evaluating Both Beta and A Shocks) ---
% =========================================================================
n_level = length(phi_level_H_range_final);
n_gap   = length(phi_gap_H_range_final);

% 期待厚生 (Expected Welfare = 0.5 * W_beta + 0.5 * W_a)
results_welfare_H_final   = zeros(n_level, n_gap);
results_welfare_F_final   = zeros(n_level, n_gap);
results_welfare_SUM_final = zeros(n_level, n_gap);

% 内訳保存用 (Betaショック単独・Aショック単独・単純和)
results_welfare_H_beta    = zeros(n_level, n_gap);
results_welfare_F_beta    = zeros(n_level, n_gap);
results_welfare_SUM_beta  = zeros(n_level, n_gap);

results_welfare_H_a       = zeros(n_level, n_gap);
results_welfare_F_a       = zeros(n_level, n_gap);
results_welfare_SUM_a     = zeros(n_level, n_gap);

results_welfare_H_total   = zeros(n_level, n_gap);

total_iterations = n_level * n_gap;
current_iteration = 0;
root_dir = pwd;

for i = 1:n_level
    for j = 1:n_gap
        current_iteration = current_iteration + 1;
        phi_gap_H_value   = phi_gap_H_range_final(j);
        phi_level_H_value = phi_level_H_range_final(i);

        fprintf('\nTrial (%d/%d): phi_gap_H = %.4g, phi_level_H = %.4g\n', ...
            current_iteration, total_iterations, phi_gap_H_value, phi_level_H_value);

        try
            % -------------------------------------------------------------
            % [Step 1] Beta Shock の計算と3段階フィルター検査
            % -------------------------------------------------------------
            [res_beta, irf_beta, Mbase_b, oobase_b] = olt_solve_core(phi_gap_H_value, phi_level_H_value, settings_beta);

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
            [res_a, irf_a, Mbase_a, oobase_a] = olt_solve_core(phi_gap_H_value, phi_level_H_value, settings_a);

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
        results_welfare_H_beta(i, j)    = w_H_b;
        results_welfare_F_beta(i, j)    = w_F_b;
        results_welfare_SUM_beta(i, j)  = w_SUM_b;

        results_welfare_H_a(i, j)       = w_H_a;
        results_welfare_F_a(i, j)       = w_F_a;
        results_welfare_SUM_a(i, j)     = w_SUM_a;

        results_welfare_H_final(i, j)   = welfare_H;
        results_welfare_F_final(i, j)   = welfare_F;
        results_welfare_SUM_final(i, j) = welfare_SUM;
        results_welfare_H_total(i, j)   = welfare_H_tot;

        if welfare_SUM > -999
            fprintf('  -> [OK] Home (Beta: %.6f, A: %.6f) => Expected H: %.6f (Total: %.6f)\n', ...
                w_H_b, w_H_a, welfare_H, welfare_H_tot);
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
% 1. Optimal for Home (Expected Welfare)
[max_welfare_H, idx_H] = max(results_welfare_H_final(:));
[row_H, col_H] = ind2sub(size(results_welfare_H_final), idx_H);
final_results.opt_phi_gap_H_H         = phi_gap_H_range_final(col_H);
final_results.opt_phi_level_H_H       = phi_level_H_range_final(row_H);
final_results.max_welfare_H           = max_welfare_H;
final_results.max_welfare_H_total     = results_welfare_H_total(row_H, col_H);
final_results.welfare_H_beta_at_H_opt = results_welfare_H_beta(row_H, col_H);
final_results.welfare_H_a_at_H_opt    = results_welfare_H_a(row_H, col_H);
final_results.welfare_F_at_H_opt      = results_welfare_F_final(row_H, col_H);
final_results.welfare_SUM_at_H_opt    = results_welfare_SUM_final(row_H, col_H);

% 2. Optimal for Foreign (Expected Welfare)
[max_welfare_F, idx_F] = max(results_welfare_F_final(:));
[row_F, col_F] = ind2sub(size(results_welfare_F_final), idx_F);
final_results.opt_phi_gap_H_F         = phi_gap_H_range_final(col_F);
final_results.opt_phi_level_H_F       = phi_level_H_range_final(row_F);
final_results.max_welfare_F           = max_welfare_F;
final_results.welfare_F_beta_at_F_opt = results_welfare_F_beta(row_F, col_F);
final_results.welfare_F_a_at_F_opt    = results_welfare_F_a(row_F, col_F);
final_results.welfare_H_at_F_opt      = results_welfare_H_final(row_F, col_F);
final_results.welfare_SUM_at_F_opt    = results_welfare_SUM_final(row_F, col_F);

% 3. Optimal for Global (Expected Sum)
[max_welfare_SUM, idx_SUM] = max(results_welfare_SUM_final(:));
[row_SUM, col_SUM] = ind2sub(size(results_welfare_SUM_final), idx_SUM);
final_results.opt_phi_gap_H_SUM         = phi_gap_H_range_final(col_SUM);
final_results.opt_phi_level_H_SUM       = phi_level_H_range_final(row_SUM);
final_results.max_welfare_SUM           = max_welfare_SUM;
final_results.welfare_SUM_beta_at_S_opt = results_welfare_SUM_beta(row_SUM, col_SUM);
final_results.welfare_SUM_a_at_S_opt    = results_welfare_SUM_a(row_SUM, col_SUM);
final_results.welfare_H_at_SUM_opt      = results_welfare_H_final(row_SUM, col_SUM);
final_results.welfare_F_at_SUM_opt      = results_welfare_F_final(row_SUM, col_SUM);

% 全データの保存 (例: olt_opt_all_data.mat)
all_data_filename = sprintf('%s_opt_all_data.mat', policy_key);
save(fullfile(base_res_dir, all_data_filename), ...
    'final_results', ...
    'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'results_welfare_H_beta', 'results_welfare_F_beta', 'results_welfare_SUM_beta', ...
    'results_welfare_H_a', 'results_welfare_F_a', 'results_welfare_SUM_a', ...
    'results_welfare_H_total', ...
    'phi_gap_H_range_final', 'phi_level_H_range_final');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n==============================================================\n');
fprintf('--- OLT Final Optimal Policy Results (Expected Welfare) ---\n');
fprintf('==============================================================\n');
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', ...
    final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H);
fprintf('   Expected Welfare H (0.5*Beta + 0.5*A) = %.6f (Total = %.6f)\n', ...
    final_results.max_welfare_H, final_results.max_welfare_H_total);
fprintf('     - Breakdown: Beta Shock H = %.6f, A Shock H = %.6f\n', ...
    final_results.welfare_H_beta_at_H_opt, final_results.welfare_H_a_at_H_opt);
fprintf('   Expected Foreign = %.6f, Expected Global Sum = %.6f\n\n', ...
    final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);

fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', ...
    final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F);
fprintf('   Resulting Expected Welfare: Home = %.6f, Foreign = %.6f, Sum = %.6f\n\n', ...
    final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);

fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_gap_H = %.4g, phi_level_H = %.4g\n', ...
    final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM);
fprintf('   Resulting Expected Welfare: Home = %.6f, Foreign = %.6f, Sum = %.6f\n', ...
    final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('==============================================================\n\n');

% =========================================================================
%% --- Plot and Save Graphs & Raw Data ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

% タイトルとパラメータ行の配置座標
% グラフ左端 (x = 0.08) に揃えて 'left' で配置し、カラーバーより十分高い位置に設定
title_pos_line1 = [0.08, 0.945, 0.88, 0.045];
title_pos_line2 = [0.08, 0.895, 0.88, 0.045];

policy_name_upper = upper(policy_key);

% 角での数字の衝突を防ぐため、X軸・Y軸の表示範囲の両端に8%の余白（パディング）を計算
x_min = min(phi_gap_H_range_final);   x_max = max(phi_gap_H_range_final);
y_min = min(phi_level_H_range_final); y_max = max(phi_level_H_range_final);
dx = max(x_max - x_min, 1e-12);
dy = max(y_max - y_min, 1e-12);
x_lim_pad = [x_min - 0.08 * dx, x_max + 0.08 * dx];
y_lim_pad = [y_min - 0.08 * dy, y_max + 0.08 * dy];

% 軸ラベル（\phi_{gap}^H と \phi_{level}^H）を配置する3D座標の計算
x_mid = mean(x_lim_pad);
y_mid = mean(y_lim_pad);
% X軸ラベル (\phi_{gap}^H) の座標 (カスタム目盛り数字のさらに右下外側)
xlab_pos_x = x_mid + 0.14 * dx;
xlab_pos_y = y_lim_pad(1) - 0.28 * dy;
% Y軸ラベル (\phi_{level}^H) の座標 (左下辺のさらに左下外側)
ylab_pos_x = x_lim_pad(1) - 0.22 * dx;
ylab_pos_y = y_mid - 0.10 * dy;

% プロット用にエラー値 (-999) を NaN に置換した配列を用意 (Z軸の潰れ防止)
plot_welfare_H   = results_welfare_H_final;   plot_welfare_H(plot_welfare_H <= -999)     = NaN;
plot_welfare_F   = results_welfare_F_final;   plot_welfare_F(plot_welfare_F <= -999)     = NaN;
plot_welfare_SUM = results_welfare_SUM_final; plot_welfare_SUM(plot_welfare_SUM <= -999) = NaN;

% --- Helper function for axis styling (Custom 3D X-Tick Labels & Axis Labels) ---
function style_3d_axes(ax, cb, x_lim_pad, y_lim_pad, dx, dy, z_limits, dz, xlab_pos_x, xlab_pos_y, ylab_pos_x, ylab_pos_y)
    % 1. 3Dグラフ本体とカラーバーの高さを抑え、上部タイトル領域との間に十分な余白を確保
    set(ax, 'Position', [0.16, 0.16, 0.58, 0.62], 'FontSize', 8);
    if ~isempty(cb) && ishandle(cb)
        set(cb, 'Position', [0.79, 0.12, 0.035, 0.68]);
    end

    % 2. phi_gap 軸 (X軸) の目盛り数字を軸からしっかり離すため、
    %    gnuplot標準の窮屈な xticklabel を消し、軸から右下外側へオフセットした位置に text() で直接描画
    xt = get(ax, 'xtick');
    xt = xt(xt >= x_lim_pad(1) & xt <= x_lim_pad(2));
    set(ax, 'xtick', xt, 'xticklabel', {});
    for k = 1:length(xt)
        text(xt(k) + 0.03 * dx, y_lim_pad(1) - 0.10 * dy, z_limits(1) - 0.03 * dz, ...
             sprintf('%g', xt(k)), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    % 3. Z軸ラベル (Welfare) をZ軸頂点の上に配置
    text(x_lim_pad(1), y_lim_pad(2), z_limits(2) + 0.08 * dz, 'Welfare', ...
         'HorizontalAlignment', 'center', 'FontSize', 10);

    % 4. X軸ラベル (\phi_{gap}^H) と Y軸ラベル (\phi_{level}^H) を3D空間の軸外側に明示配置
    text(xlab_pos_x, xlab_pos_y, z_limits(1) - 0.08 * dz, '\phi_{gap}^H', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Interpreter', 'tex');
    text(ylab_pos_x, ylab_pos_y, z_limits(1) - 0.06 * dz, '\phi_{level}^H', ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Interpreter', 'tex');
end

% --- Graph 1: Home Expected Welfare ---
h_home = figure('Name', 'OLT: Home Expected Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_H);
cb_H = colorbar; view(-37.5, 25); hold on;
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_H = zlim;
dz_H = max(z_limits_H(2) - z_limits_H(1), 1e-6);
zlim([z_limits_H(1) - 0.08 * dz_H, z_limits_H(2) + 0.08 * dz_H]);
z_limits_H = zlim;
plot3(final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H, ...
      'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_H, final_results.opt_phi_gap_H_H], ...
      [final_results.opt_phi_level_H_H, final_results.opt_phi_level_H_H], ...
      [z_limits_H(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
style_3d_axes(gca, cb_H, x_lim_pad, y_lim_pad, dx, dy, z_limits_H, dz_H, xlab_pos_x, xlab_pos_y, ylab_pos_x, ylab_pos_y);
hold off;
line1_H = sprintf('%s: Home Expected Welfare', policy_name_upper);
line2_H = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H, ...
    final_results.welfare_H_beta_at_H_opt, final_results.welfare_H_a_at_H_opt);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_home_filename = sprintf('%s_opt_graph_home.pdf', policy_key);
print(h_home, fullfile(base_res_dir, graph_home_filename), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_gap_H_range_final; data_y = phi_level_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_gap_H_H, final_results.opt_phi_level_H_H, final_results.max_welfare_H];
mat_home_filename = sprintf('%s_opt_data_for_plotting_home.mat', policy_key);
save(fullfile(base_res_dir, mat_home_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Expected Welfare ---
h_foreign = figure('Name', 'OLT: Foreign Expected Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_F);
cb_F = colorbar; view(-37.5, 25); hold on;
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_F = zlim;
dz_F = max(z_limits_F(2) - z_limits_F(1), 1e-6);
zlim([z_limits_F(1) - 0.08 * dz_F, z_limits_F(2) + 0.08 * dz_F]);
z_limits_F = zlim;
plot3(final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F, ...
      'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_F, final_results.opt_phi_gap_H_F], ...
      [final_results.opt_phi_level_H_F, final_results.opt_phi_level_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
style_3d_axes(gca, cb_F, x_lim_pad, y_lim_pad, dx, dy, z_limits_F, dz_F, xlab_pos_x, xlab_pos_y, ylab_pos_x, ylab_pos_y);
hold off;
line1_F = sprintf('%s: Foreign Expected Welfare', policy_name_upper);
line2_F = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F, ...
    final_results.welfare_F_beta_at_F_opt, final_results.welfare_F_a_at_F_opt);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_foreign_filename = sprintf('%s_opt_graph_foreign.pdf', policy_key);
print(h_foreign, fullfile(base_res_dir, graph_foreign_filename), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_gap_H_F, final_results.opt_phi_level_H_F, final_results.max_welfare_F];
mat_foreign_filename = sprintf('%s_opt_data_for_plotting_foreign.mat', policy_key);
save(fullfile(base_res_dir, mat_foreign_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global (Sum) Expected Welfare ---
h_global = figure('Name', 'OLT: Global Expected Welfare', 'Visible', 'off');
surf(phi_gap_H_range_final, phi_level_H_range_final, plot_welfare_SUM);
cb_G = colorbar; view(-37.5, 25); hold on;
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_G = zlim;
dz_G = max(z_limits_G(2) - z_limits_G(1), 1e-6);
zlim([z_limits_G(1) - 0.08 * dz_G, z_limits_G(2) + 0.08 * dz_G]);
z_limits_G = zlim;
plot3(final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM, ...
      'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_gap_H_SUM, final_results.opt_phi_gap_H_SUM], ...
      [final_results.opt_phi_level_H_SUM, final_results.opt_phi_level_H_SUM], ...
      [z_limits_G(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
style_3d_axes(gca, cb_G, x_lim_pad, y_lim_pad, dx, dy, z_limits_G, dz_G, xlab_pos_x, xlab_pos_y, ylab_pos_x, ylab_pos_y);
hold off;
line1_G = sprintf('%s: Global Expected Welfare', policy_name_upper);
line2_G = sprintf('\\phi_{gap}^H = %.4g, \\phi_{level}^H = %.4g, Exp. Welfare = %.4f (\\beta: %.4f, a: %.4f)', ...
    final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM, ...
    final_results.welfare_SUM_beta_at_S_opt, final_results.welfare_SUM_a_at_S_opt);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 10);
graph_global_filename = sprintf('%s_opt_graph_global.pdf', policy_key);
print(h_global, fullfile(base_res_dir, graph_global_filename), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_gap_H_SUM, final_results.opt_phi_level_H_SUM, final_results.max_welfare_SUM];
mat_global_filename = sprintf('%s_opt_data_for_plotting_global.mat', policy_key);
save(fullfile(base_res_dir, mat_global_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('OLT Policy expected welfare optimization has finished.\n');
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');