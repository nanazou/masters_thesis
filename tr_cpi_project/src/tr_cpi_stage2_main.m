% src/tr_cpi_stage2_main.m
% このスクリプトは Stage 1 の潜在産出量データを読み込み、tr_cpi_stage2_solve_core.m を実行します。
% 成果物を results > stage2 > [shock] 直下に集約します。

% 注: このスクリプトは src フォルダ内に配置され、実行されます。

% clear;
close all;
% clc;

% --- [1] パス設定 ---
% スクリプトの場所（src フォルダ）を基準にパスを設定
mpath = fileparts(mfilename('fullpath'));
addpath(fullfile(mpath, 'stage2'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_ options_;

if ~exist('shock_scenario', 'var')
    shock_scenario = 'beta'; 
end

if strcmp(shock_scenario, 'beta')
    shock_folder = 'beta_shock';
elseif strcmp(shock_scenario, 'a')
    shock_folder = 'a_shock';
else
    error("Invalid shock_scenario. Choose 'beta' or 'a'.");
end

% ▼▼▼【計算に使用するパラメータ設定 (TR CPI STAGE2用)】▼▼▼
phi_pi_H_value = 83;
phi_y_H_value  = 0.5;
% ▲▲▲【設定完了】▲▲▲

% 出力パスを project_root/results/stage2/[shock] 直下に設定
output_root_folder = fullfile(mpath, '..', 'results', 'stage2', shock_folder);

if ~exist(output_root_folder, 'dir'), mkdir(output_root_folder); end

% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

%% --- [2] Stage 1 データの読み込み ---
% Stage 1 の成果物 (results/stage1/[shock]/potential_output.mat) を参照
stage1_file = fullfile(mpath, '..', 'results', 'stage1', shock_folder, 'potential_output.mat');

if ~exist(stage1_file, 'file')
    error('Stage 1 の結果ファイルが見つかりません: %s\n先に tr_cpi_stage1_main.m を実行してください。', stage1_file);
end

Stage1Data = load(stage1_file); 

if isfield(Stage1Data, 'eps_y_H_potential')
    eps_y_H_potential_vec = Stage1Data.eps_y_H_potential(:);
else
    error('読み込んだファイル内に eps_y_H_potential が見つかりません。');
end
fprintf('Stage 1 潜在産出量データをロードしました。シナリオ: %s\n', shock_scenario);

%% --- [3] シミュレーションの設定 ---
% 独自名に変更した統一ファイル名 ('tr_cpi_stage2_model', 'tr_cpi_stage2_model_zlb') を参照
settings.modnam = 'tr_cpi_stage2_model';
settings.modnamstar = 'tr_cpi_stage2_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';
settings.nperiods = 400;
settings.maxiter = 100;

% ショックの定義 (主ショック + Stage1から引き継いだ潜在産出量パス)
if strcmp(shock_scenario, 'beta')
    settings.irfshock = char('eps_beta_H', 'eps_y_H_potential');
    primary_shock_magnitude = 0.02;
elseif strcmp(shock_scenario, 'a')
    settings.irfshock = char('eps_a_H', 'eps_y_H_potential');
    primary_shock_magnitude = -0.1;
end

% ショックシーケンスの構築 (2つの外生変数を同時に投入)
% 第1列: 主ショック (t=1のみ), 第2列: 自然産出量パス (全期間)
shockssequence = zeros(settings.nperiods, 2);
shockssequence(1, 1) = primary_shock_magnitude; 
len_data = min(settings.nperiods, length(eps_y_H_potential_vec));
shockssequence(1:len_data, 2) = eps_y_H_potential_vec(1:len_data);

settings.shockssequence = shockssequence;

%% --- [4] 【核心】stage2/tr_cpi_stage2_solve_core.m の呼び出し ---
fprintf('Stage 2 (Sticky Price + ZLB) シミュレーションを開始します...\n');
[welfares, irf_results, Mbase_, oobase_] = tr_cpi_stage2_solve_core(phi_pi_H_value, phi_y_H_value, settings);

% 取得した厚生変数を個別に展開
welfare_H_with_delta_pie = welfares.welfare_H_with_delta_pie;
welfare_H_with_delta_lin = welfares.welfare_H_with_delta_lin;
welfare_H_without_delta_pie = welfares.welfare_H_without_delta_pie;
welfare_H_without_delta_lin = welfares.welfare_H_without_delta_lin;
welfare_F_with_delta_pie = welfares.welfare_F_with_delta_pie;
welfare_F_with_delta_lin = welfares.welfare_F_with_delta_lin;
welfare_F_without_delta_pie = welfares.welfare_F_without_delta_pie;
welfare_F_without_delta_lin = welfares.welfare_F_without_delta_lin;

% irf_results構造体の中身を現在のワークスペースに展開
fnames = fieldnames(irf_results);
for k = 1:length(fnames)
    eval([fnames{k}, ' = irf_results.', fnames{k}, ';']);
end

%% --- [5] グラフの描画 ---
fprintf('グラフを作成し、PDFに保存中...\n');

% 表示変数のリスト (統計ベース指数および全ショック変数を網羅)
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
        'p_H_W_bar', 'p_F_W_star_bar', ... % 統計ベースのCPI
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
        'y_H_potential', ...
        'eps_beta_H', 'eps_beta_F', ...
        'eps_a_H', 'eps_a_F', ...
        'eps_tau_H', 'eps_tau_F', ...
        'eps_i_H', 'eps_i_F', ...
        'eps_y_H_potential', ...
        'Delta_H', 'Delta_F', ...
    };

graphics_toolkit("gnuplot");

% 出力ファイル名の設定 (results > stage2 > [shock] > graphs.pdf)
output_filename = fullfile(output_root_folder, 'graphs.pdf'); 

if isfile(output_filename), delete(output_filename); end

plots_per_page = 6;
rows = 2;
cols = 3;
num_vars = length(var_list);
nperiods = settings.nperiods;
irfshock = settings.irfshock;

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
            shock_row_idx = 0;
            for r = 1:size(irfshock, 1)
                if strcmp(strtrim(irfshock(r,:)), var_name)
                    shock_row_idx = r; break;
                end
            end
            if shock_row_idx > 0
                plot(1:nperiods, shockssequence(:, shock_row_idx), 'b-', 'LineWidth', 1.5);
            else
                plot(1:nperiods, zeros(nperiods, 1), 'b-', 'LineWidth', 1.5);
            end
        else
            plot(1:nperiods, eval([var_name, '_pie']), 'b-', 'LineWidth', 1.5);
            hold on;
            plot(1:nperiods, eval([var_name, '_lin']), 'r--', 'LineWidth', 1.5);
        end
        grid on;
        xlim([0, nperiods]);
        
        % --- タイトル生成ロジック (一切の要約・省略なし) ---
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
        elseif strcmp(var_name, 'p_F_W_star_bar')
            latex_title = 'Normalized CPI (\bar{p}^{F \rightarrow W*})';
        elseif strcmp(var_name, 'pi_H_W')
            latex_title = 'CPI Inflation (\pi^{H \rightarrow W})';
        elseif strcmp(var_name, 'pi_F_W_star')
            latex_title = 'CPI Inflation (\pi^{F \rightarrow W*})';
        elseif strcmp(var_name, 'p_H_tilde')
            latex_title = '\tilde{p}^H';
        elseif strcmp(var_name, 'p_F_star_tilde')
            latex_title = '\tilde{p}^{F*}';
        elseif strcmp(var_name, 'i_H_notional')
            latex_title = 'i^{H,notional}';
        elseif endsWith(var_name, '_bar') || contains(var_name, '_bar_') || endsWith(var_name, '_star_bar') % 修正
            % 命名規則変更に対応した LaTeX ラベル生成
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
        elseif strcmp(var_name, 'y_H_potential')
            latex_title = 'y^{H,pot}';
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
                superscript = strrep(superscript, 'star', '*'); 
                latex_title = [base '^{' superscript '}'];
            else
                latex_title = base;
            end
        end
        
        set(gca, 'FontSize', 8);
        title(latex_title, 'Interpreter', 'tex');

        if plot_index == 1
            legend('ZLB (OccBin)', 'Non-ZLB', 'Location', 'northeast');
        end
        
        % 厚生値のグラフ内表示
        if strcmp(var_name, 'utility_H_with_delta')
            xlim_vals = get(gca, 'XLim'); ylim_vals = get(gca, 'YLim');
            x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
            y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
            y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
            text_H_pie_str = sprintf('Welfare (with Delta, ZLB): %.4f', welfare_H_with_delta_pie);
            text_H_lin_str = sprintf('Welfare (with Delta, Non-ZLB): %.4f', welfare_H_with_delta_lin);
            text(x_pos, y_pos_bottom, text_H_lin_str, 'Color', 'red', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            text(x_pos, y_pos_top, text_H_pie_str, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
        elseif strcmp(var_name, 'utility_H_without_delta')
            xlim_vals = get(gca, 'XLim'); ylim_vals = get(gca, 'YLim');
            x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
            y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
            y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
            text_H_pie_str = sprintf('Welfare (without Delta, ZLB): %.4f', welfare_H_without_delta_pie);
            text_H_lin_str = sprintf('Welfare (without Delta, Non-ZLB): %.4f', welfare_H_without_delta_lin);
            text(x_pos, y_pos_bottom, text_H_lin_str, 'Color', 'red', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            text(x_pos, y_pos_top, text_H_pie_str, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
    
    if i == 1
        print(h, output_filename, '-dpdf');
    else
        print(h, output_filename, '-dpdf', '-append');
    end
    close(h);
end

%% --- [6] データの保存 (.mat) ---
fprintf('全てのシミュレーション結果を.matファイルに保存中...\n');

vars_to_save = {};
for i = 1:length(var_list)
    vn = var_list{i};
    if ~strncmp(vn, 'eps_', 4) 
        if exist([vn, '_pie'], 'var')
            vars_to_save{end+1} = [vn, '_pie']; 
            vars_to_save{end+1} = [vn, '_lin']; 
        end
    end
end

% 厚生変数の追加
w_list = {'welfare_H_with_delta_pie','welfare_H_with_delta_lin','welfare_F_with_delta_pie','welfare_F_with_delta_lin',...
          'welfare_H_without_delta_pie','welfare_H_without_delta_lin','welfare_F_without_delta_pie','welfare_F_without_delta_lin'};
vars_to_save = [vars_to_save, w_list];

% 成果物フォルダ直下に保存
full_save_path = fullfile(output_root_folder, 'data_for_plotting.mat');
save(full_save_path, vars_to_save{:});
fprintf('全てのデータが "%s" に保存されました。\n', full_save_path);

%% --- [7] CSVへの書き出し ---
fprintf('CSVファイルへの書き出しを開始します...\n');
try
    S = load(full_save_path);
    header = {}; data_matrix = []; 
    for i = 1:length(vars_to_save)
        vn = vars_to_save{i};
        if ~strncmp(vn, 'welfare_', 8)
            if isfield(S, vn) && isvector(S.(vn)) && length(S.(vn)) == nperiods
                header{end+1} = vn;
                data_matrix = [data_matrix, S.(vn)(:)]; 
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
        fprintf(fid_scalar, 'Welfare (Utility) Calculations (Stage 2 Sticky - Scenario: %s):\n\n', shock_scenario);
        fprintf(fid_scalar, '--- with Price Dispersion Cost (Delta included) ---\n');
        fprintf(fid_scalar, 'welfare_H_with_delta_pie (ZLB): %.8f\n', S.welfare_H_with_delta_pie);
        fprintf(fid_scalar, 'welfare_H_with_delta_lin (Non-ZLB): %.8f\n', S.welfare_H_with_delta_lin);
        fprintf(fid_scalar, 'welfare_F_with_delta_pie (ZLB): %.8f\n', S.welfare_F_with_delta_pie);
        fprintf(fid_scalar, 'welfare_F_with_delta_lin (Non-ZLB): %.8f\n', S.welfare_F_with_delta_lin);
        fprintf(fid_scalar, '\n');
        fprintf(fid_scalar, '--- without Price Dispersion Cost (Delta = 0) ---\n');
        fprintf(fid_scalar, 'welfare_H_without_delta_pie (ZLB): %.8f\n', S.welfare_H_without_delta_pie);
        fprintf(fid_scalar, 'welfare_H_without_delta_lin (Non-ZLB): %.8f\n', S.welfare_H_without_delta_lin);
        fprintf(fid_scalar, 'welfare_F_without_delta_pie (ZLB): %.8f\n', S.welfare_F_without_delta_pie);
        fprintf(fid_scalar, 'welfare_F_without_delta_lin (Non-ZLB): %.8f\n', S.welfare_F_without_delta_lin);
        fclose(fid_scalar);
    end
catch ME
    fprintf('CSV書き出しエラー: %s\n', ME.message);
end

fprintf('DONE. Stage 2 [%s] の処理が完了しました。フォルダ: %s\n', shock_scenario, output_root_folder);