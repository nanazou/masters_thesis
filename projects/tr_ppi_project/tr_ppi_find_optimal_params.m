% tr_ppi_find_optimal_params.m

% このスクリプトは、TR PPI政策の最適パラメータを探索し、
% 厚生の表面グラフを出力します。計算の核心部は stage2/tr_ppi_stage2_solve_core.m を呼び出し、
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

% 各フォルダへのパスを追加
addpath(fullfile(mpath, 'src', 'stage2'));
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: tr_ppi_project -> tr_ppi)
[~, proj_folder_name] = fileparts(mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_ options_;

fprintf('====================================================\n');
fprintf('Starting optimization for TR PPI Policy (ZLB Aware)...\n');
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
shock_scenario = 'a'; 

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

% --- 3. 粗い探索 (coarse mode) の設定 --- a shock (TR_CPIと揃えた境界ズーム: 6 × 3 = 18点)
phi_pi_H_range_coarse = [1.0050, 1.0052, 1.0054, 1.0056, 1.0058, 1.0060];
phi_y_H_range_coarse  = [0, 1e-6, 2e-6];

% --- 3. 粗い探索 (coarse mode) の設定 --- beta
%phi_pi_H_range_coarse = [0, 1e-4];
%phi_y_H_range_coarse  = (0.5 : 0.5 : 8.0) * 1e16;

% --- 3. 粗い探索 (coarse mode) の設定 --- 旧 a final
%phi_pi_H_range_coarse   = 0:1:1;
%phi_y_H_range_coarse    = 0.005:0.0001:0.0053;

% --- 4. 詳細な探索 (fine mode) の設定 ---
center_phi_pi_H   = 1.5; 
center_phi_y_H    = 0.5; 
fine_search_radius_pi_H   = 0.4; 
fine_search_step_pi_H     = 0.1; 
fine_search_radius_y_H    = 0.4; 
fine_search_step_y_H      = 0.1;

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
% --- 設定エリアここまで ---

% --- 出力フォルダの定義と作成 ---
base_res_dir = fullfile(mpath, 'results', 'optimal_params', shock_folder);
if ~exist(base_res_dir, 'dir'), mkdir(base_res_dir); end

%% --- [1] Stage 1 (潜在産出量) パスの自動計算 (POLT構成と統一) ---
fprintf('\n--- Stage 1: Calculating Target Potential Path for %s shock ---\n', shock_scenario);

root_dir = pwd;
cd(fullfile(mpath, 'src', 'stage1'));

% 自然産出量の算出用ダミーパラメータ (TR用に統一)
phi_pi_H_dummy = 1.5; phi_y_H_dummy = 0.5;
assignin('base', 'phi_pi_H_value', phi_pi_H_dummy);
assignin('base', 'phi_y_H_value', phi_y_H_dummy);

% Stage 1 モデルの実行
eval(['dynare tr_ppi_stage1_model.mod noclearall nolog']);

% ショックパスの計算 (simult_ を使用)
ex_flex = zeros(400, M_.exo_nbr);
s_idx = strmatch(target_shock, M_.exo_names, 'exact');
ex_flex(1, s_idx) = primary_shock_magnitude;
y_flex_sim = simult_(M_, options_, oo_.dr.ys, oo_.dr, ex_flex, 1);
y_H_idx = strmatch('y_H', M_.endo_names, 'exact');
eps_y_H_potential_vec = y_flex_sim(y_H_idx, 2:end)' - oo_.dr.ys(y_H_idx);

% 元のディレクトリに戻る
cd(root_dir);

% 計算したパスを保存 (例: tr_ppi_opt_stage1_potential_path_a_shock.mat)
stage1_path_filename = sprintf('%s_opt_stage1_potential_path_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, stage1_path_filename), 'eps_y_H_potential_vec');
fprintf('Stage 1 完了。Stage 2 最適化を開始します...\n');

%% --- [2] Simulation Settings Structure ---
% 【修正箇所】実際のファイル名 tr_ppi_stage2_model に合わせて修正
settings.modnam = 'tr_ppi_stage2_model';
settings.modnamstar = 'tr_ppi_stage2_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.nperiods = 400;
settings.maxiter = 100;

% ショック定義 (1×2行ベクトルにして400回ループを回避し高速化)
settings.irfshock = char(target_shock, 'eps_y_H_potential');
settings.shockssequence = [primary_shock_magnitude, eps_y_H_potential_vec(1)];

% ★減衰パラメータの定義を追加
if strcmp(shock_scenario, 'beta')
    assignin('base', 'rho_potential_value', 0.977);
elseif strcmp(shock_scenario, 'a')
    assignin('base', 'rho_potential_value', 0.950);
end

% =========================================================================
%% --- Grid Setup based on Mode ---
% =========================================================================
if strcmp(search_mode, 'coarse')
    fprintf('--- Starting COARSE grid search ---\n');
    phi_pi_H_range_final   = phi_pi_H_range_coarse;
    phi_y_H_range_final = phi_y_H_range_coarse;
    plot_title_prefix = 'Coarse Search';
elseif strcmp(search_mode, 'fine')
    fprintf('--- Starting FINE grid search around (pi_H=%.4g, y_H=%.4g) ---\n', center_phi_pi_H, center_phi_y_H);
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

        fprintf('\nTrial (%d/%d): phi_pi_H = %.4g, phi_y_H = %.4g ... ', ...
            current_iteration, total_iterations, phi_pi_H_value, phi_y_H_value);

        try
            % 共通エンジン tr_ppi_stage2_solve_core.m を呼び出し (4つの戻り値を取得)
            [res, irf_res, Mbase_, oobase_] = tr_ppi_stage2_solve_core(phi_pi_H_value, phi_y_H_value, settings);

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

            % 厚生データの抽出
            welfare_H   = res.welfare_H_with_delta_pie;
            welfare_F   = res.welfare_F_with_delta_pie;
            welfare_SUM = res.welfare_SUM_with_delta_pie;
        catch ME
            fprintf('Error during trial: %s\n', ME.message);
            cd(root_dir);
            welfare_H   = -999; welfare_F   = -999; welfare_SUM = -999;
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

% 全データの保存 (例: tr_ppi_opt_all_data_a_shock.mat)
all_data_filename = sprintf('%s_opt_all_data_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, all_data_filename), ...
    'final_results', 'results_welfare_H_final', 'results_welfare_F_final', 'results_welfare_SUM_final', ...
    'phi_pi_H_range_final', 'phi_y_H_range_final');

% =========================================================================
%% --- Print Final Results to Log ---
% =========================================================================
fprintf('\n--- TR PPI Final Optimal Policy Results (Shock: %s) ---\n', shock_scenario);
fprintf('1. Optimal for Home:\n');
fprintf('   Parameters: phi_pi_H = %.4g, phi_y_H = %.4g\n', final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.max_welfare_H, final_results.welfare_F_at_H_opt, final_results.welfare_SUM_at_H_opt);
fprintf('2. Optimal for Foreign:\n');
fprintf('   Parameters: phi_pi_H = %.4g, phi_y_H = %.4g\n', final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n\n', final_results.welfare_H_at_F_opt, final_results.max_welfare_F, final_results.welfare_SUM_at_F_opt);
fprintf('3. Optimal for Global (Sum):\n');
fprintf('   Parameters: phi_pi_H = %.4g, phi_y_H = %.4g\n', final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM);
fprintf('   Resulting Welfare: Home = %.4f, Foreign = %.4f, Sum = %.4f\n', final_results.welfare_H_at_SUM_opt, final_results.welfare_F_at_SUM_opt, final_results.max_welfare_SUM);
fprintf('--------------------------------------------------------------\n\n');

% =========================================================================
%% --- Plot and Save Graphs & Raw Data ---
% =========================================================================
fprintf('Plotting and saving results...\n');
graphics_toolkit("gnuplot");

title_pos_line1 = [0, 0.94, 1, 0.05];
title_pos_line2 = [0, 0.89, 1, 0.05];

% 探索モードの簡潔な表記
mode_label = [upper(search_mode(1)), lower(search_mode(2:end))]; % 'Coarse' or 'Fine'

% 角での数字の衝突を防ぐため、X軸・Y軸の表示範囲の両端に8%の余白（パディング）を計算
x_min = min(phi_pi_H_range_final); x_max = max(phi_pi_H_range_final);
y_min = min(phi_y_H_range_final);  y_max = max(phi_y_H_range_final);
dx = max(x_max - x_min, 1e-15);
dy = max(y_max - y_min, 1e-15);
x_lim_pad = [x_min - 0.08 * dx, x_max + 0.08 * dx];
y_lim_pad = [y_min - 0.08 * dy, y_max + 0.08 * dy];

% プロット用にエラー値 (-999) を NaN に置換した配列を用意 (Z軸の潰れ防止)
plot_welfare_H   = results_welfare_H_final;   plot_welfare_H(plot_welfare_H <= -999)     = NaN;
plot_welfare_F   = results_welfare_F_final;   plot_welfare_F(plot_welfare_F <= -999)     = NaN;
plot_welfare_SUM = results_welfare_SUM_final; plot_welfare_SUM(plot_welfare_SUM <= -999) = NaN;

% --- Graph 1: Home Welfare ---
h_home = figure('Name', 'TR PPI: Home Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, plot_welfare_H);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_H = zlim;
dz_H = max(z_limits_H(2) - z_limits_H(1), 1e-6);
zlim([z_limits_H(1) - 0.08 * dz_H, z_limits_H(2) + 0.08 * dz_H]);
z_limits_H = zlim;
plot3(final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_pi_H_H, final_results.opt_phi_pi_H_H], ...
      [final_results.opt_phi_y_H_H, final_results.opt_phi_y_H_H], ...
      [z_limits_H(1), final_results.max_welfare_H], 'k--', 'LineWidth', 1);
% Z軸ラベル(Welfare)は数字との重なりを防ぐためZ軸の頂点上に横書きで配置
text(x_lim_pad(1), y_lim_pad(2), z_limits_H(2) + 0.08 * dz_H, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{\\pi}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{y}^H'), 'FontSize', 10);
line1_H = sprintf('TR_{PPI}: Home Welfare (%s shock, %s)', shock_scenario, mode_label);
line2_H = sprintf('\\phi_{\\pi}^H = %.4g, \\phi_{y}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H);
annotation('textbox', title_pos_line1, 'String', line1_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_H, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_home_filename = sprintf('%s_opt_graph_home_%s.pdf', policy_key, shock_folder);
print(h_home, fullfile(base_res_dir, graph_home_filename), '-dpdf');
close(h_home);

% Home個別データ保存
data_x = phi_pi_H_range_final; data_y = phi_y_H_range_final; data_z = results_welfare_H_final;
opt_pt = [final_results.opt_phi_pi_H_H, final_results.opt_phi_y_H_H, final_results.max_welfare_H];
mat_home_filename = sprintf('%s_opt_data_for_plotting_home_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_home_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 2: Foreign Welfare ---
h_foreign = figure('Name', 'TR PPI: Foreign Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, plot_welfare_F);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_F = zlim;
dz_F = max(z_limits_F(2) - z_limits_F(1), 1e-6);
zlim([z_limits_F(1) - 0.08 * dz_F, z_limits_F(2) + 0.08 * dz_F]);
z_limits_F = zlim;
plot3(final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_pi_H_F, final_results.opt_phi_pi_H_F], ...
      [final_results.opt_phi_y_H_F, final_results.opt_phi_y_H_F], ...
      [z_limits_F(1), final_results.max_welfare_F], 'k--', 'LineWidth', 1);
text(x_lim_pad(1), y_lim_pad(2), z_limits_F(2) + 0.08 * dz_F, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{\\pi}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{y}^H'), 'FontSize', 10);
line1_F = sprintf('TR_{PPI}: Foreign Welfare (%s shock, %s)', shock_scenario, mode_label);
line2_F = sprintf('\\phi_{\\pi}^H = %.4g, \\phi_{y}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F);
annotation('textbox', title_pos_line1, 'String', line1_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_F, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_foreign_filename = sprintf('%s_opt_graph_foreign_%s.pdf', policy_key, shock_folder);
print(h_foreign, fullfile(base_res_dir, graph_foreign_filename), '-dpdf');
close(h_foreign);

% Foreign個別データ保存
data_z = results_welfare_F_final;
opt_pt = [final_results.opt_phi_pi_H_F, final_results.opt_phi_y_H_F, final_results.max_welfare_F];
mat_foreign_filename = sprintf('%s_opt_data_for_plotting_foreign_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_foreign_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

% --- Graph 3: Global (Sum) Welfare ---
h_global = figure('Name', 'TR PPI: Global Welfare', 'Visible', 'off');
surf(phi_pi_H_range_final, phi_y_H_range_final, plot_welfare_SUM);
colorbar; view(-37.5, 25); hold on;
set(gca, 'Position', [0.16, 0.16, 0.62, 0.66], 'FontSize', 8);
xlim(x_lim_pad); ylim(y_lim_pad);
z_limits_G = zlim;
dz_G = max(z_limits_G(2) - z_limits_G(1), 1e-6);
zlim([z_limits_G(1) - 0.08 * dz_G, z_limits_G(2) + 0.08 * dz_G]);
z_limits_G = zlim;
plot3(final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
plot3([final_results.opt_phi_pi_H_SUM, final_results.opt_phi_pi_H_SUM], ...
      [final_results.opt_phi_y_H_SUM, final_results.opt_phi_y_H_SUM], ...
      [z_limits_G(1), final_results.max_welfare_SUM], 'k--', 'LineWidth', 1);
text(x_lim_pad(1), y_lim_pad(2), z_limits_G(2) + 0.08 * dz_G, 'Welfare', ...
     'HorizontalAlignment', 'center', 'FontSize', 10);
hold off;
xlabel(sprintf('\n\\phi_{\\pi}^H'), 'FontSize', 10);
ylabel(sprintf('\n\\phi_{y}^H'), 'FontSize', 10);
line1_G = sprintf('TR_{PPI}: Global Welfare (%s shock, %s)', shock_scenario, mode_label);
line2_G = sprintf('\\phi_{\\pi}^H = %.4g, \\phi_{y}^H = %.4g, Welfare = %.4f', ...
                  final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM);
annotation('textbox', title_pos_line1, 'String', line1_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
annotation('textbox', title_pos_line2, 'String', line2_G, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
graph_global_filename = sprintf('%s_opt_graph_global_%s.pdf', policy_key, shock_folder);
print(h_global, fullfile(base_res_dir, graph_global_filename), '-dpdf');
close(h_global);

% Global個別データ保存
data_z = results_welfare_SUM_final;
opt_pt = [final_results.opt_phi_pi_H_SUM, final_results.opt_phi_y_H_SUM, final_results.max_welfare_SUM];
mat_global_filename = sprintf('%s_opt_data_for_plotting_global_%s.mat', policy_key, shock_folder);
save(fullfile(base_res_dir, mat_global_filename), 'data_x', 'data_y', 'data_z', 'opt_pt');

fprintf('\n====================================================\n');
fprintf('TR PPI Policy optimization has finished for %s.\n', shock_scenario);
fprintf('All results (PDF and .mat) saved in: %s\n', base_res_dir);
fprintf('====================================================\n');