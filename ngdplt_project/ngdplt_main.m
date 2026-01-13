% ngdplt_main.m
% このスクリプトは ngdplt_solve_core.m を呼び出し、特定のパラメータでの結果を詳細に出力します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

global M_ oo_;

%% --- 共通シミュレーション設定 ---
% 独自名に変更した統一ファイル名 ('ngdplt_model', 'ngdplt_model_zlb') を参照
settings.modnam = 'ngdplt_model';
settings.modnamstar = 'ngdplt_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';

settings.nperiods = 400;
settings.maxiter = 100;

% ▼▼▼【計算に使用するパラメータ設定 (NGDPLT用)】▼▼▼
phi_gap_H_value = 6.5;
phi_level_H_value = 13.8;
% ▲▲▲【設定完了】▲▲▲

%% --- シナリオ定義 ---
% 1: beta_shock (需要ショック), 2: a_shock (供給ショック)
scenarios = struct();

% シナリオ1: Beta Shock
scenarios(1).name = 'beta_shock';
scenarios(1).irfshock = 'eps_beta_H';
scenarios(1).shock_val = 0.02; % 正のショック(割引因子上昇=需要減退)

% シナリオ2: A Shock
scenarios(2).name = 'a_shock';
scenarios(2).irfshock = 'eps_a_H';
scenarios(2).shock_val = -0.1; % 負のショック(生産性低下)

%% --- ループ実行開始 ---
for s_idx = 1:length(scenarios)
    
    % --- 現在のシナリオ設定の読み込み ---
    current_scen = scenarios(s_idx);
    fprintf('\n====================================================\n');
    fprintf('NGDPLT シナリオ実行中: %s (%s = %.2f)\n', current_scen.name, current_scen.irfshock, current_scen.shock_val);
    fprintf('====================================================\n');

    % 設定構造体の更新
    settings.shockssequence = [ current_scen.shock_val; zeros(settings.nperiods-1, 1) ];
    settings.irfshock = current_scen.irfshock;
    
    % --- 出力先フォルダの設定 (results > [shock] 直下に変更) ---
    output_root_folder = fullfile(mpath, 'results', current_scen.name); 
    if ~exist(output_root_folder, 'dir'), mkdir(output_root_folder); end

    % --- 【核心】ngdplt_solve_core.m の呼び出し ---
    [welfares, irf_results, Mbase_, oobase_] = ngdplt_solve_core(phi_gap_H_value, phi_level_H_value, settings);

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

    %% --- グラフの描画 ---
    fprintf('グラフを作成し、PDFに保存中 (%s)...\n', current_scen.name);

% 表示変数のリスト (chi_H および 統計ベース物価指数を補完)
    var_list = { ...
        'e_slash_star', ...
        'lambda_H', ...
        'lambda_H_slash_star', ...
        'lambda_F_slash_star', ...
        'lambda_F', ...
        'utility_H_with_delta', 'utility_F_with_delta', ...
        'utility_H_without_delta', 'utility_F_without_delta', ...
        'gamma_H', 'chi_H', ... % ★ chi_H (目標パス) を追加
        'p_H', 'p_F_star', ...
        'i_H_notional', ...
        'i_H', 'i_F', ...
        'c_H_W', 'c_F_W', ...
        'c_H_H', 'c_F_H', ...
        'c_H_F', 'c_F_F', ...
        'y_H', 'y_F', ...
        'l_H', 'l_F', ...
        'p_H_W', 'p_F_W_star', ...
        'p_H_W_bar', 'p_F_W_star_bar', ... % ★ 統計ベースの消費者物価指数を追加
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
        'eps_chi_H', ...
        'eps_i_F', ...
        'Delta_H', 'Delta_F', ...
    };

    % gnuplotを使用してPDFを出力
    graphics_toolkit("gnuplot");
    output_filename = fullfile(output_root_folder, 'graphs.pdf');
    if isfile(output_filename), delete(output_filename); end

    plots_per_page = 6;
    rows = 2;
    cols = 3;
    num_vars = length(var_list);
    nperiods = settings.nperiods;
    irfshock = settings.irfshock;
    shockssequence = settings.shockssequence;

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
                if strcmp(var_name, strtrim(irfshock))
                    plot(1:nperiods, shockssequence, 'b-', 'LineWidth', 1.5);
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
            
            % --- タイトル生成ロジック (省略なし) ---
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
            elseif strcmp(var_name, 'pi_H_W')
                latex_title = '\pi^{H \to W}';
            elseif strcmp(var_name, 'pi_F_W_star')
                latex_title = '\pi^{F \to W*}';
            elseif strcmp(var_name, 'c_H_W')
                latex_title = 'c^{H \to W}';
            elseif strcmp(var_name, 'c_F_W')
                latex_title = 'c^{F \to W*}';
            elseif strcmp(var_name, 'p_H_W')
                latex_title = 'p^{H \to W}';
            elseif strcmp(var_name, 'p_F_W_star')
                latex_title = 'p^{F \to W*}';
            elseif strcmp(var_name, 'c_H_H')
                latex_title = 'c^{H \to H}';
            elseif strcmp(var_name, 'c_F_H')
                latex_title = 'c^{F \to H}';
            elseif strcmp(var_name, 'c_H_F')
                latex_title = 'c^{H \to F*}';
            elseif strcmp(var_name, 'c_F_F')
                latex_title = 'c^{F \to F*}';
            elseif strcmp(var_name, 'p_H_tilde')
                latex_title = '\tilde{p}^H';
            elseif strcmp(var_name, 'p_F_star_tilde')
                latex_title = '\tilde{p}^{F*}';
            elseif strcmp(var_name, 'i_H_notional')
                latex_title = 'i^{H,notional}';
            elseif endsWith(var_name, '_bar')
                temp_name = regexprep(var_name, '_bar$', '');
                parts = strsplit(temp_name, '_');
                base = parts{1}; 
                if length(parts) > 1
                    superscript = strjoin(parts(2:end), ''); 
                    superscript = strrep(superscript, 'star', '*'); 
                    latex_title = [base, '\_bar^{', superscript, '}'];
                else
                    latex_title = [base, '\_bar'];
                end
            elseif strcmp(var_name, 'lambda_H_slash_star')
                latex_title = '\lambda^{H/*}'; 
            elseif strcmp(var_name, 'lambda_H')
                latex_title = '\lambda^H';
            elseif strcmp(var_name, 'lambda_F')
                latex_title = '\lambda^F';          
            elseif strcmp(var_name, 'lambda_F_slash_star')
                latex_title = '\lambda^{F/*}';
            elseif strcmp(var_name, 'pi_H')
                latex_title = '\pi^H';
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
                    % [重要] star を * に置換
                    superscript = strrep(superscript, 'star', '*'); 
                    latex_title = [base '^{' superscript '}'];
                else
                    latex_title = base;
                end
            end
            
            set(gca, 'FontSize', 8);
            title(latex_title, 'Interpreter', 'tex');

            if plot_index == 1
                legend('ZLB', 'Non-ZLB', 'Location', 'northeast');
            end
            
            % 厚生の表示
            if strcmp(var_name, 'utility_H_with_delta')
                xlim_vals = get(gca, 'XLim');
                ylim_vals = get(gca, 'YLim');
                x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
                y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
                y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
                text_H_pie = sprintf('Welfare (with Delta, ZLB): %.4f', welfare_H_with_delta_pie);
                text_H_lin = sprintf('Welfare (with Delta, Non-ZLB): %.4f', welfare_H_with_delta_lin);
                text(x_pos, y_pos_bottom, text_H_lin, 'Color', 'red', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
                text(x_pos, y_pos_top, text_H_pie, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            elseif strcmp(var_name, 'utility_H_without_delta')
                xlim_vals = get(gca, 'XLim');
                ylim_vals = get(gca, 'YLim');
                x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
                y_pos_top = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
                y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
                text_H_pie = sprintf('Welfare (without Delta, ZLB): %.4f', welfare_H_without_delta_pie);
                text_H_lin = sprintf('Welfare (without Delta, Non-ZLB): %.4f', welfare_H_without_delta_lin);
                text(x_pos, y_pos_bottom, text_H_lin, 'Color', 'red', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
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

    %% --- データの保存 ---
    fprintf('全てのシミュレーション結果を.matファイルに保存中 (%s)...\n', current_scen.name);

    % 保存用リスト
    vars_to_save = {};
    for i = 1:length(var_list)
        curr_v = var_list{i};
        if ~strncmp(curr_v, 'eps_', 4) 
            vars_to_save{end+1} = [curr_v, '_pie']; 
            vars_to_save{end+1} = [curr_v, '_lin']; 
        end
    end

    vars_to_save{end+1} = 'welfare_H_with_delta_pie';
    vars_to_save{end+1} = 'welfare_H_with_delta_lin';
    vars_to_save{end+1} = 'welfare_F_with_delta_pie';
    vars_to_save{end+1} = 'welfare_F_with_delta_lin';
    vars_to_save{end+1} = 'welfare_H_without_delta_pie';
    vars_to_save{end+1} = 'welfare_H_without_delta_lin';
    vars_to_save{end+1} = 'welfare_F_without_delta_pie';
    vars_to_save{end+1} = 'welfare_F_without_delta_lin';

    % Output フォルダを介さず直下の data_for_plotting.mat に保存
    full_save_path = fullfile(output_root_folder, 'data_for_plotting.mat');

    save(full_save_path, vars_to_save{:});
    fprintf('データ保存完了: %s\n', full_save_path);

    %% --- CSVへの書き出し ---
    fprintf('CSVファイルへの書き出しを開始します (%s)...\n', current_scen.name);

    try
        S = load(full_save_path);
        header = {}; 
        data_matrix = []; 
        
        for i = 1:length(vars_to_save)
            var_name = vars_to_save{i};
            if ~strncmp(var_name, 'welfare_', 8)
                header{end+1} = var_name;
                data_matrix = [data_matrix, S.(var_name)(:)]; 
            end
        end
        
        % 保存パスの構成 (results/[shock]/data_for_plotting.csv)
        csv_save_path = fullfile(output_root_folder, 'data_for_plotting.csv');

        fid = fopen(csv_save_path, 'w');
        if fid == -1, error('CSVファイルを開けませんでした'); end
        
        fprintf(fid, '%s', header{1}); 
        for i = 2:length(header), fprintf(fid, ',%s', header{i}); end
        fprintf(fid, '\n'); 
        fclose(fid);
        
        dlmwrite(csv_save_path, data_matrix, '-append', 'delimiter', ',', 'precision', '%.8g');
        
        % --- 厚生データのテキスト保存 ---
        scalar_save_path = fullfile(output_root_folder, 'data_for_plotting_welfare.txt');
        fid_scalar = fopen(scalar_save_path, 'w');
        if fid_scalar ~= -1
            fprintf(fid_scalar, 'Welfare Calculations (%s):\n\n', current_scen.name);
            
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
            fprintf('CSV/Text保存完了: %s\n', current_scen.name);
        end
        
    catch ME
        fprintf('CSVへの書き出し中にエラーが発生しました:\n');
        disp(ME.message);
    end
    
    close all; 

end

fprintf('\n全てのシナリオ（beta_shock, a_shock）の実行が完了しました。\n');