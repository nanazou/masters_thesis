% src/tr_cpi_stage1_main.m
% このスクリプトは、選択されたショックに基づいて潜在産出量（Potential Output）を計算します。
% 成果物を results > stage1 > [shock] 直下に集約します。

% 注: このスクリプトは src フォルダ内に配置され、実行されます。

% clear; % 統合メインスクリプトから呼ばれる際、変数を維持するためコメントアウト
close all;
% clc;

% --- [1] パス設定 ---
% スクリプトの場所（src フォルダ）を基準にパスを設定
mpath = fileparts(mfilename('fullpath'));
addpath(fullfile(mpath, 'stage1'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_ options_;

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- 設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼

% 外部（統合メインスクリプト）から shock_scenario が指定されていない場合のデフォルト値
if ~exist('shock_scenario', 'var')
    % 'beta' : 需要ショック (eps_beta_H = 0.02)
    % 'a'    : 供給ショック (eps_a_H = -0.1)
    shock_scenario = 'beta'; 
end

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

% --- [2] ショックパラメータおよびフォルダ名の自動設定 ---
if strcmp(shock_scenario, 'beta')
    target_shock = 'eps_beta_H';
    shock_value  = 0.02; 
    shock_folder = 'beta_shock'; % フォルダ構成用
elseif strcmp(shock_scenario, 'a')
    target_shock = 'eps_a_H';
    shock_value  = -0.1; 
    shock_folder = 'a_shock';    % フォルダ構成用
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% 出力パスを project_root/results/stage1/[shock] 直下に設定
output_root_folder = fullfile(mpath, '..', 'results', 'stage1', shock_folder);

% フォルダの自動生成 (存在しない場合)
if ~exist(output_root_folder, 'dir'), mkdir(output_root_folder); end

% シミュレーション設定の構築
settings.stage = 1;
settings.modnam = 'tr_cpi_stage1_model'; % src/stage1/tr_cpi_stage1_model.mod を参照
settings.nperiods = 400;
settings.irfshock = target_shock;
settings.shockssequence = shock_value;

% Stage 1 用ダミーパラメータ (形式統一のために保持)
phi_gap_H_value = 0;
phi_level_H_value = 0;

% --- [3] 【核心】stage1/tr_cpi_stage1_solve_core.m の呼び出し ---
fprintf('Stage 1 (Flexible Price) シミュレーションを開始します。シナリオ: %s (%s)\n', shock_scenario, shock_folder);
[welfares, irf_results, Mbase_, oobase_] = tr_cpi_stage1_solve_core(phi_gap_H_value, phi_level_H_value, settings);

% 厚生変数を個別に展開
welfare_H_with_delta_pie = welfares.welfare_H_with_delta_pie;
welfare_H_with_delta_lin = welfares.welfare_H_with_delta_lin;
welfare_H_without_delta_pie = welfares.welfare_H_without_delta_pie;
welfare_H_without_delta_lin = welfares.welfare_H_without_delta_lin;
welfare_F_with_delta_pie = welfares.welfare_F_with_delta_pie;
welfare_F_with_delta_lin = welfares.welfare_F_with_delta_lin;
welfare_F_without_delta_pie = welfares.welfare_F_without_delta_pie;
welfare_F_without_delta_lin = welfares.welfare_F_without_delta_lin;

% IRF結果構造体の全フィールドを現在のワークスペースに展開
fnames = fieldnames(irf_results);
for k = 1:length(fnames)
    eval([fnames{k}, ' = irf_results.', fnames{k}, ';']);
end

% --- [4] Stage 2 連携用データの保存 ---
% Stage 2 のショックとして使用する潜在産出量データ
% 伸縮的価格モデルなので y_H_pie がそのまま潜在産出量 (偏差) となる
eps_y_H_potential = y_H_pie; 

% Stage 2 が直接参照するファイル名で保存
save_path_potential = fullfile(output_root_folder, 'potential_output.mat');
save(save_path_potential, 'eps_y_H_potential');

fprintf('Stage 1 潜在産出量を保存しました: %s\n', save_path_potential);

% --- [5] グラフの描画 ---
fprintf('グラフを作成し、PDFに保存中...\n');

% 表示変数のリスト (全ショック変数および統計ベース物価指数を網羅)
    var_list = { ...
        'e_slash_star', ... 
        'lambda_H', ...
        'lambda_H_slash_star', ...
        'lambda_F_slash_star', ...
        'lambda_F', ...
        'utility_H_with_delta', 'utility_F_with_delta', ...
        'utility_H_without_delta', 'utility_F_without_delta', ... 
        'gamma_H', ...
        'p_H', 'p_F_star', ...
        'i_H_notional', ...
        'i_H', 'i_F', ...
        'c_H_W', 'c_F_W', ...
        'c_H_H', 'c_F_H', ...
        'c_H_F', 'c_F_F', ...
        'y_H', 'y_F', ...
        'l_H', 'l_F', ...
        'p_H_W', 'p_F_W_star', ...
        'p_H_W_bar', 'p_F_W_star_bar', ... % 統計ベースの消費者物価指数
        'p_H_bar', 'p_F_star_bar', ...
        'p_H_bar_y_H', 'p_H_W_c_H_W', ... 
        'pi_H', 'pi_F_star', ...
        'pi_H_W', 'pi_F_W_star', ...
        'p_H_tilde', 'p_F_star_tilde', ...
        'v_H', 'v_F', ...
        'w_H', 'w_F', ...
        'p_F_bar', 'p_H_star_bar', ...
        't_H', 't_F', ...
        'b_H', 'b_F_star', ... 
        'beta_H', 'beta_F', ...
        'a_H', 'a_F', ...
        'tau_H', 'tau_F', ...
        'eps_beta_H', 'eps_beta_F', ...
        'eps_a_H', 'eps_a_F', ...
        'eps_tau_H', 'eps_tau_F', ...
        'eps_i_H', 'eps_i_F', ...
        'Delta_H', 'Delta_F', ...
    };

graphics_toolkit("gnuplot");

% 出力ファイル名の設定 (results > stage1 > [shock] > graphs.pdf)
output_filename = fullfile(output_root_folder, 'graphs.pdf'); 

if isfile(output_filename), delete(output_filename); end

plots_per_page = 6;
rows = 2;
cols = 3;
num_vars = length(var_list);
nperiods = settings.nperiods;

for i = 1:plots_per_page:num_vars
    h = figure('Visible', 'off', 'PaperType', 'a4', 'PaperOrientation', 'landscape', 'PaperUnits', 'centimeters');
    
    graph_width = 25;
    graph_height = 18;
    left_margin = (29.7 - graph_width) / 2;
    bottom_margin = (21.0 - graph_height) / 2;
    set(h, 'PaperPosition', [left_margin, bottom_margin, graph_width, graph_height]);
    
    for j = 1:plots_per_page
        plot_index = i + j - 1;
        if plot_index > num_vars, break; end
        subplot(rows, cols, j);
        
        var_name = var_list{plot_index};
        
        if strncmp(var_name, 'eps_', 4)
            data_shock = zeros(nperiods, 1);
            if strcmp(var_name, settings.irfshock)
                data_shock(1) = shock_value;
            end
            plot(1:nperiods, data_shock, 'b-', 'LineWidth', 1.5);
        else
            plot(1:nperiods, eval([var_name, '_pie']), 'b-', 'LineWidth', 1.5);
            hold on;
            plot(1:nperiods, eval([var_name, '_lin']), 'r--', 'LineWidth', 1.5);
        end
        grid on;
        xlim([0, nperiods]);
        
        % --- タイトル生成ロジック ---
        if strncmp(var_name, 'eps_', 4)
            parts = strsplit(var_name, '_');
            middle = parts{2};
            tex_middle = middle; 
            switch middle
                case 'beta', tex_middle = '\beta';
                case 'tau',  tex_middle = '\tau';
                case 'chi',  tex_middle = '\chi';
            end
            if length(parts) > 2
                endpoint = parts{3};
                latex_title = ['\epsilon^{', tex_middle, '^{', endpoint, '}}'];
            else
                latex_title = ['\epsilon_{', tex_middle, '}'];
            end
        elseif strcmp(var_name, 'p_H_bar_y_H') 
            latex_title = 'Nominal GDP (\bar{p}^H y^H)';
        elseif strcmp(var_name, 'p_H_W_c_H_W') 
            latex_title = 'Nominal Total Cons (p^{H \to W} c^{H \to W})';
        elseif strcmp(var_name, 'p_H_W_bar')
            latex_title = 'Normalized CPI (\bar{p}^{H \rightarrow W})';
        elseif strcmp(var_name, 'p_F_W_star_bar') % 修正
            latex_title = 'Normalized CPI (\bar{p}^{F \rightarrow W*})';
        elseif strcmp(var_name, 'pi_H_W')
            latex_title = 'CPI Inflation (\pi^{H \rightarrow W})';
        elseif strcmp(var_name, 'pi_F_W_star')
            latex_title = 'CPI Inflation (\pi^{F \rightarrow W*})';
        elseif strcmp(var_name, 'p_H_tilde')
            latex_title = '\tilde{p}^H';
        elseif strcmp(var_name, 'p_F_star_tilde') % 修正
            latex_title = '\tilde{p}^{F*}';
        elseif strcmp(var_name, 'i_H_notional')
            latex_title = 'i^{H,notional}';
        elseif endsWith(var_name, '_bar') || contains(var_name, '_bar_') || contains(var_name, '_star_bar') % 修正
            temp_name = var_name;
            parts = strsplit(temp_name, '_');
            base = parts{1};
            superscript = '';
            if contains(var_name, '_W')
                superscript = [superscript, 'H \rightarrow W'];
            end
            if contains(var_name, 'star')
                if isempty(superscript), superscript = '*'; else superscript = [superscript, '*']; end
            end
            
            if isempty(superscript)
                latex_title = [base, '\_bar'];
            else
                latex_title = [base, '\_bar^{', superscript, '}'];
            end
        elseif strcmp(var_name, 'lambda_H_slash_star')
            latex_title = '\lambda^{H/*}';
        elseif strcmp(var_name, 'lambda_F')
            latex_title = '\lambda^F';
        elseif strcmp(var_name, 'lambda_F_slash_star')
            latex_title = '\lambda^{F/*}';
        elseif strcmp(var_name, 'pi_F_star')
            latex_title = '\pi^{F*}';
        elseif strcmp(var_name, 'b_F_star')
            latex_title = 'b^{F*}';
        elseif strcmp(var_name, 'e_slash_star')
            latex_title = 'e^{/*}';
        elseif strcmp(var_name, 'utility_H_without_delta')
            latex_title = 'utility^H (without \Delta)';
        elseif strcmp(var_name, 'utility_F_without_delta')
            latex_title = 'utility^F (without \Delta)';
        elseif strcmp(var_name, 'utility_H_with_delta')
            latex_title = 'utility^H (with \Delta)';
        elseif strcmp(var_name, 'utility_F_with_delta')
            latex_title = 'utility^F (with \Delta)';
        else
            parts = strsplit(var_name, '_');
            base = parts{1};
            switch base
                case 'pi', base = '\pi';
                case 'lambda', base = '\lambda';
                case 'gamma', base = '\gamma';
                case 'beta', base = '\beta';
                case 'mu', base = '\mu';
                case 'chi', base = '\chi';
                case 'tau', base = '\tau';
                case 'utility', base = 'utility';
                case 'Delta', base = '\Delta'; 
            end
            if length(parts) > 1
                superscript = strjoin(parts(2:end), '_');
                superscript = strrep(superscript, '_', '\rightarrow'); 
                superscript = strrep(superscript, '\rightarrowstar', '*'); 
                latex_title = [base '^{' superscript '}'];
            else
                latex_title = base;
            end
        end
        
        set(gca, 'FontSize', 8);
        title(latex_title, 'Interpreter', 'tex');

        if plot_index == 1
            legend('Stage 1 (Flex)', 'Location', 'northeast');
        end
        
        % 厚生値の表示
        if strcmp(var_name, 'utility_H_with_delta')
            xlim_vals = get(gca, 'XLim'); ylim_vals = get(gca, 'YLim');
            x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
            y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
            text_H_pie = sprintf('Welfare (with Delta): %.4f', welfare_H_with_delta_pie);
            text(x_pos, y_pos_top, text_H_pie, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
        elseif strcmp(var_name, 'utility_H_without_delta')
            xlim_vals = get(gca, 'XLim'); ylim_vals = get(gca, 'YLim');
            x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
            y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
            text_H_pie = sprintf('Welfare (without Delta): %.4f', welfare_H_without_delta_pie);
            text(x_pos, y_pos_top, text_H_pie, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
    
    page_num = floor((i-1) / plots_per_page) + 1;
    if page_num == 1
        print(h, output_filename, '-dpdf');
    else
        print(h, output_filename, '-dpdf', '-append');
    end
    close(h);
end

% --- [6] 全データの保存 (.mat) ---
fprintf('全てのシミュレーション結果を.matファイルに保存中...\n');

vars_to_save = {};
for i = 1:length(var_list)
    v_name = var_list{i};
    if ~strncmp(v_name, 'eps_', 4) 
        if exist([v_name, '_pie'], 'var')
            vars_to_save{end+1} = [v_name, '_pie']; 
            vars_to_save{end+1} = [v_name, '_lin']; 
        end
    end
end

% 厚生変数のリスト作成
w_list = {'welfare_H_with_delta_pie','welfare_H_with_delta_lin','welfare_F_with_delta_pie','welfare_F_with_delta_lin',...
          'welfare_H_without_delta_pie','welfare_H_without_delta_lin','welfare_F_without_delta_pie','welfare_F_without_delta_lin'};
vars_to_save = [vars_to_save, w_list];

% 簡潔な名前 (data_for_plotting.mat) で保存
full_save_path = fullfile(output_root_folder, 'data_for_plotting.mat');
save(full_save_path, vars_to_save{:});
fprintf('全てのデータが "%s" に保存されました。\n', full_save_path);

% --- [7] CSVへの書き出し ---
fprintf('CSVファイルへの書き出しを開始します...\n');
try
    S = load(full_save_path);
    header = {}; data_matrix = []; 
    for i = 1:length(vars_to_save)
        v_name = vars_to_save{i};
        if ~strncmp(v_name, 'welfare_', 8)
            if isfield(S, v_name) && isvector(S.(v_name)) && length(S.(v_name)) == nperiods
                header{end+1} = v_name;
                data_matrix = [data_matrix, S.(v_name)(:)]; 
            end
        end
    end
    
    csv_save_path = fullfile(output_root_folder, 'data_for_plotting.csv');
    fid = fopen(csv_save_path, 'w');
    if fid == -1, error('CSVファイルを開けませんでした: %s', csv_save_path); end
    fprintf(fid, '%s', header{1}); 
    for i = 2:length(header), fprintf(fid, ',%s', header{i}); end
    fprintf(fid, '\n'); 
    fclose(fid);
    dlmwrite(csv_save_path, data_matrix, '-append', 'delimiter', ',', 'precision', '%.8g');
    
    % --- [8] 厚生データのテキスト保存 ---
    scalar_save_path = fullfile(output_root_folder, 'data_for_plotting_welfare.txt');
    fid_scalar = fopen(scalar_save_path, 'w');
    if fid_scalar ~= -1
        fprintf(fid_scalar, 'Welfare (Utility) Calculations (Stage 1 Flex - Scenario: %s):\n\n', shock_scenario);
        fprintf(fid_scalar, '--- with Price Dispersion Cost (Delta included) ---\n');
        fprintf(fid_scalar, 'welfare_H_with_delta_pie: %.8f\n', S.welfare_H_with_delta_pie);
        fprintf(fid_scalar, 'welfare_H_with_delta_lin: %.8f\n', S.welfare_H_with_delta_lin);
        fprintf(fid_scalar, 'welfare_F_with_delta_pie: %.8f\n', S.welfare_F_with_delta_pie);
        fprintf(fid_scalar, 'welfare_F_with_delta_lin: %.8f\n', S.welfare_F_with_delta_lin);
        fprintf(fid_scalar, '\n');
        fprintf(fid_scalar, '--- without Price Dispersion Cost (Delta = 0) ---\n');
        fprintf(fid_scalar, 'welfare_H_without_delta_pie: %.8f\n', S.welfare_H_without_delta_pie);
        fprintf(fid_scalar, 'welfare_H_without_delta_lin: %.8f\n', S.welfare_H_without_delta_lin);
        fprintf(fid_scalar, 'welfare_F_without_delta_pie: %.8f\n', S.welfare_F_without_delta_pie);
        fprintf(fid_scalar, 'welfare_F_without_delta_lin: %.8f\n', S.welfare_F_without_delta_lin);
        fclose(fid_scalar);
    end
catch ME
    fprintf('CSV書き出しエラー: %s\n', ME.message);
end

fprintf('DONE. Stage 1 [%s] の処理が完了しました。フォルダ: %s\n', shock_scenario, output_root_folder);