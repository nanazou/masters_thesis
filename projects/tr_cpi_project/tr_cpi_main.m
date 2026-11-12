% tr_cpi_main.m
%
% 目的:
%   1. results/optimal_params/tr_cpi_opt_all_data.mat から最適パラメータを自動読み込み
%      （または手動指定したパラメータを使用）
%   2. 需要ショック（beta）と供給ショック（a）のそれぞれについて、
%      src/tr_cpi_stage1_main.m (潜在産出量算出) と
%      src/tr_cpi_stage2_main.m (共通パラメータでの政策シミュレーション) を順番に実行
%   3. 両ショックの厚生と期待厚生 (Expected Welfare) のサマリーを出力

clear;
close all;
clc;

% --- [1] 実行環境の固定 ---
root_mpath = fileparts(mfilename('fullpath'));
cd(root_mpath);

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: tr_cpi_project -> tr_cpi)
[~, proj_folder_name] = fileparts(root_mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_ options_;

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- パラメータ読み込み設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% true  : tr_cpi_find_optimal_params.m で求めた最適パラメータ (.mat) を自動読み込み
% false : 下の manual_phi_* に書いた値を手動で使用
use_auto_opt_params = true;

% 手動実行時（use_auto_opt_params = false のとき）に使う値
manual_phi_pi_H = 150;
manual_phi_y_H  = 1.0;
% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

opt_mat_path = fullfile(root_mpath, 'results', 'optimal_params', sprintf('%s_opt_all_data.mat', policy_key));

if use_auto_opt_params
    if ~isfile(opt_mat_path)
        error(['最適パラメータファイルが見つかりません:\n  %s\n', ...
               '先に %s_find_optimal_params.m を実行するか、use_auto_opt_params = false に設定してください。'], ...
               opt_mat_path, policy_key);
    end
    opt_data = load(opt_mat_path, 'final_results');
    phi_pi_H_common  = opt_data.final_results.opt_phi_pi_H_H;
    phi_y_H_common   = opt_data.final_results.opt_phi_y_H_H;
    param_source_str = sprintf('Auto-loaded from %s_opt_all_data.mat', policy_key);
else
    phi_pi_H_common  = manual_phi_pi_H;
    phi_y_H_common   = manual_phi_y_H;
    param_source_str = 'Manually specified in tr_cpi_main.m';
end

fprintf('========================================================================\n');
fprintf('  %s MAIN SIMULATION START (Stage 1 & Stage 2)\n', upper(policy_key));
fprintf('  Parameter Source : %s\n', param_source_str);
fprintf('  Applied Params   : phi_pi_H = %.6g, phi_y_H = %.8g\n', phi_pi_H_common, phi_y_H_common);
fprintf('========================================================================\n');

% 実行したいショックシナリオのリスト ('beta' = 需要ショック, 'a' = 供給ショック)
scenarios = {'beta', 'a'};

% サマリー記録用構造体
summary_welfares = struct();

% シナリオ（Beta と A）のループを開始
for s = 1:length(scenarios)
    
    % 子スクリプト（tr_cpi_stage1_main, tr_cpi_stage2_main）に渡す共通変数をセット
    shock_scenario  = scenarios{s};
    phi_pi_H_custom = phi_pi_H_common;
    phi_y_H_custom  = phi_y_H_common;
    
    fprintf('\n************************************************************************\n');
    fprintf('  【実行中シナリオ (%d/%d)】: [%s]\n', s, length(scenarios), upper(shock_scenario));
    fprintf('  使用パラメータ: phi_pi_H = %.6g, phi_y_H = %.8g\n', phi_pi_H_custom, phi_y_H_custom);
    fprintf('************************************************************************\n');
    
    % --- Step 1: Stage 1 (Flexible Price / Potential Output) の実行 ---
    fprintf('\n>> Step 1: src/tr_cpi_stage1_main.m を実行しています (潜在産出量パスの計算)...\n');
    try
        run(fullfile(root_mpath, 'src', 'tr_cpi_stage1_main.m'));
        cd(root_mpath);
    catch ME
        cd(root_mpath);
        fprintf('\n[エラー発生] Stage 1 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue;
    end
    
    % --- Step 2: Stage 2 (Sticky Price + ZLB / Policy Simulation) の実行 ---
    fprintf('\n>> Step 2: src/tr_cpi_stage2_main.m を実行しています (政策シミュレーション)...\n');
    try
        run(fullfile(root_mpath, 'src', 'tr_cpi_stage2_main.m'));
        cd(root_mpath);
        summary_welfares.(shock_folder) = welfares;
    catch ME
        cd(root_mpath);
        fprintf('\n[エラー発生] Stage 2 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue;
    end
    
    fprintf('\n--- シナリオ [%s] の全計算と保存が完了しました ---\n', upper(shock_scenario));
end

%% --- 両ショック統合の実行サマリー表示 ---
if isfield(summary_welfares, 'beta_shock') && isfield(summary_welfares, 'a_shock')
    w_H_beta_pie = summary_welfares.beta_shock.welfare_H_with_delta_pie;
    w_H_a_pie    = summary_welfares.a_shock.welfare_H_with_delta_pie;
    w_H_exp_pie  = 0.5 * (w_H_beta_pie + w_H_a_pie);

    fprintf('\n========================================================================\n');
    fprintf('  全てのシナリオ（beta_shock, a_shock）の実行が完了しました。\n');
    fprintf('  [適用パラメータ] phi_pi_H = %.6g, phi_y_H = %.8g (%s)\n', ...
        phi_pi_H_common, phi_y_H_common, param_source_str);
    fprintf('  [Home厚生 (Delta有り)] Beta: %.6f, A: %.6f => 期待厚生: %.6f\n', ...
        w_H_beta_pie, w_H_a_pie, w_H_exp_pie);
    fprintf('========================================================================\n');
end

fprintf('結果の確認先 (project_root基準):\n');
fprintf('  - 最適パラメータ探索結果 : results/optimal_params/\n');
fprintf('  - Stage 1 結果           : results/stage1/[shock_folder]/\n');
fprintf('  - Stage 2 結果           : results/stage2/[shock_folder]/\n');
fprintf('========================================================================\n');