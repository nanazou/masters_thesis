% tr_ppi_main.m
%
% 目的: 
%   1. 需要ショック（Betaショック）に対する Stage 1 と Stage 2 を実行
%   2. 供給ショック（Aショック）に対する Stage 1 と Stage 2 を実行
%
% 動作の仕組み:
%   このスクリプトは 'shock_scenario' 変数を定義した上で、
%   src フォルダ内の 'tr_ppi_stage1_main.m' と 'tr_ppi_stage2_main.m' を順番に呼び出します。
%   各ステージの結果は、results フォルダ内の各ショック専用フォルダに集約されます。

clear;
close all;
clc;

% --- [1] 実行環境の固定 ---
% スクリプトが存在する場所（Rootフォルダ）をカレントディレクトリに設定
root_mpath = fileparts(mfilename('fullpath'));
cd(root_mpath);

% 実行したいショックシナリオのリストを定義
% 'beta' = 需要ショック, 'a' = 供給ショック
scenarios = {'beta', 'a'};

fprintf('========================================================================\n');
fprintf('         Taylor Rule (TR PPI) プロジェクト 全自動シミュレーション開始\n');
fprintf('========================================================================\n');

% シナリオ（Beta と A）のループを開始
for s = 1:length(scenarios)
    
    % 子スクリプト（tr_ppi_stage1_main, tr_ppi_stage2_main）に渡す共通変数をセット
    shock_scenario = scenarios{s};
    
    fprintf('\n************************************************************************\n');
    fprintf('  【実行中シナリオ】: [%s]\n', upper(shock_scenario));
    fprintf('************************************************************************\n');
    
    % --- Step 1: Stage 1 (Flexible Price / Potential Output) の実行 ---
    fprintf('\n>> Step 1: src/tr_ppi_stage1_main.m を実行しています (潜在産出量パスの計算)...\n');
    try
        % 'run' コマンドにより、現在のワークスペースの変数(shock_scenario)を維持したまま実行
        run(fullfile(root_mpath, 'src', 'tr_ppi_stage1_main.m'));
    catch ME
        fprintf('\n[エラー発生] Stage 1 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue; % 次のショックシナリオへ進む
    end
    
    % --- Step 2: Stage 2 (Sticky Price + ZLB / Policy Simulation) の実行 ---
    fprintf('\n>> Step 2: src/tr_ppi_stage2_main.m を実行しています (政策シミュレーション)...\n');
    try
        % Stage 1 で保存された 'potential_output.mat' が自動的にロードされます
        run(fullfile(root_mpath, 'src', 'tr_ppi_stage2_main.m'));
    catch ME
        fprintf('\n[エラー発生] Stage 2 (%s) でエラーが発生しました:\n', shock_scenario);
        disp(ME.message);
        continue; % 次のショックシナリオへ進む
    end
    
    fprintf('\n--- シナリオ [%s] の全計算と保存が完了しました ---\n', upper(shock_scenario));
    
end

fprintf('\n========================================================================\n');
fprintf('         全てのシミュレーション工程が正常に終了しました。\n');
fprintf('========================================================================\n');
fprintf('結果の確認先 (project_root基準):\n');
fprintf('  - Stage 1 結果 (潜在産出量データ/グラフ):\n');
fprintf('    results/stage1/[shock_scenario]/\n');
fprintf('  - Stage 2 結果 (最終IRFグラフ/CSV/厚生データ):\n');
fprintf('    results/stage2/[shock_scenario]/\n');
fprintf('========================================================================\n');