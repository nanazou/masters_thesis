% polt_main.m
% POLT (Potential Output Level Targeting) Comprehensive Execution Script
%
% このスクリプトは results/optimal_params/polt_opt_all_data.mat から
% 最適パラメータを自動で読み込み（または手動指定値を使用し）、
% 需要ショック（beta）と供給ショック（a）のそれぞれについて
% src/polt_stage1_main.m と src/polt_stage2_main.m を順番に実行します。

clear;
close all;
clc;

% --- [1] 実行環境の固定 ---
root_mpath = fileparts(mfilename('fullpath'));
cd(root_mpath);

[~, proj_folder_name] = fileparts(root_mpath);
policy_key = strrep(proj_folder_name, '_project', '');

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- パラメータ読み込み設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% true  : polt_find_optimal_params.m で求めた最適パラメータ (.mat) を自動読み込み
% false : 下の manual_phi_* に書いた値を手動で使用
use_auto_opt_params = true;

% 手動実行時（use_auto_opt_params = false のとき）に使う値
manual_phi_gap_H   = 1000000;
manual_phi_level_H = 1000000;
% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

opt_mat_path = fullfile(root_mpath, 'results', 'optimal_params', sprintf('%s_opt_all_data.mat', policy_key));

if use_auto_opt_params
    if ~isfile(opt_mat_path)
        error(['最適パラメータファイルが見つかりません:\n  %s\n', ...
               '先に %s_find_optimal_params.m を実行するか、use_auto_opt_params = false に設定してください。'], ...
               opt_mat_path, policy_key);
    end
    opt_data = load(opt_mat_path, 'final_results');
    phi_gap_H_common   = opt_data.final_results.opt_phi_gap_H_H;
    phi_level_H_common = opt_data.final_results.opt_phi_level_H_H;
    param_source_str   = sprintf('Auto-loaded from %s_opt_all_data.mat', policy_key);
else
    phi_gap_H_common   = manual_phi_gap_H;
    phi_level_H_common = manual_phi_level_H;
    param_source_str   = 'Manually specified in polt_main.m';
end

% 実行するショックシナリオ
scenarios = {'beta', 'a'};
summary_welfares = struct();

fprintf('========================================================================\n');
fprintf('  %s (Potential Output Level Targeting) 全自動シミュレーション開始\n', upper(policy_key));
fprintf('  Parameter Source : %s\n', param_source_str);
fprintf('  Applied Params   : phi_gap_H = %.6g, phi_level_H = %.8g\n', phi_gap_H_common, phi_level_H_common);
fprintf('========================================================================\n');

% シナリオ（Beta と A）のループを開始
for s = 1:length(scenarios)
    
    % 子スクリプト（polt_stage1_main, polt_stage2_main）に渡す共通変数をセット
    shock_scenario = scenarios{s};
    
    fprintf('\n************************************************************************\n');
    fprintf('  【実行中シナリオ (%d/%d)】: [%s]\n', s, length(scenarios), upper(shock_scenario));
    fprintf('************************************************************************\n');
    
    % --- Step 1: Stage 1 (Flexible Price / Potential Output) の実行 ---
    fprintf('\n>> Step 1: src/polt_stage1_main.m を実行しています (潜在産出量パスの計算)...\n');
    try
        run(fullfile(root_mpath, 'src', 'polt_stage1_main.m'));
    catch ME
        fprintf('\n[エラー発生] Stage 1 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue;
    end
    
    % --- Step 2: Stage 2 (Sticky Price + ZLB / Policy Simulation) の実行 ---
    fprintf('\n>> Step 2: src/polt_stage2_main.m を実行しています (政策シミュレーション)...\n');
    try
        run(fullfile(root_mpath, 'src', 'polt_stage2_main.m'));
        summary_welfares.([shock_scenario, '_shock']) = welfares;
    catch ME
        fprintf('\n[エラー発生] Stage 2 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue;
    end
    
    fprintf('\n--- シナリオ [%s] の全計算と保存が完了しました ---\n', upper(shock_scenario));
end

%% --- 両ショック統合の実行サマリー表示 ---
fprintf('\n========================================================================\n');
fprintf('  全てのシミュレーション工程が正常に終了しました。\n');
fprintf('  [適用パラメータ] phi_gap_H = %.6g, phi_level_H = %.8g (%s)\n', ...
    phi_gap_H_common, phi_level_H_common, param_source_str);
if isfield(summary_welfares, 'beta_shock') && isfield(summary_welfares, 'a_shock')
    w_H_beta_pie = summary_welfares.beta_shock.welfare_H_with_delta_pie;
    w_H_a_pie    = summary_welfares.a_shock.welfare_H_with_delta_pie;
    w_H_exp_pie  = 0.5 * (w_H_beta_pie + w_H_a_pie);
    fprintf('  [Home厚生 (Delta有り)] Beta: %.6f, A: %.6f => 期待厚生: %.6f\n', ...
        w_H_beta_pie, w_H_a_pie, w_H_exp_pie);
end
fprintf('------------------------------------------------------------------------\n');
fprintf('  - Stage 1 結果: results/stage1/[shock_folder]/\n');
fprintf('  - Stage 2 結果: results/stage2/[shock_folder]/\n');
fprintf('========================================================================\n');