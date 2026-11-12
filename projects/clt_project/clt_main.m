% clt_main.m
% このスクリプトは results/optimal_params/clt_opt_all_data.mat から
% 最適パラメータを自動で読み込み、clt_solve_core.m を呼び出して
% beta_shock と a_shock の結果を詳細に出力します。

clear;
close all;
clc;

% --- パス設定 ---
mpath = fileparts(mfilename('fullpath'));
cd(mpath);
addpath(fullfile(mpath, 'src'));
% addpath('C:\dynare\Occbin_update-master\toolkit_files');

% 現在のプロジェクト名からポリシーキーを自動抽出 (例: clt_project -> clt)
[~, proj_folder_name] = fileparts(mpath);
policy_key = strrep(proj_folder_name, '_project', '');

global M_ oo_;

% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% --- パラメータ読み込み設定エリア ---
% ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
% true  : find_optimal_params.m で求めた最適パラメータ (.mat) を自動読み込み
% false : 下の manual_phi_* に書いた値を手動で使用
use_auto_opt_params = true;

% 手動実行時（use_auto_opt_params = false のとき）に使う値
manual_phi_gap_H   = 0.02;
manual_phi_level_H = 1.20;
% ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

opt_mat_path = fullfile(mpath, 'results', 'optimal_params', sprintf('%s_opt_all_data.mat', policy_key));

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
    param_source_str   = 'Manually specified in clt_main.m';
end

fprintf('==============================================================\n');
fprintf('  %s MAIN SIMULATION START\n', upper(policy_key));
fprintf('  Parameter Source : %s\n', param_source_str);
fprintf('  Applied Params   : phi_gap_H = %.6g, phi_level_H = %.8g\n', phi_gap_H_common, phi_level_H_common);
fprintf('==============================================================\n');

%% --- 共通シミュレーション設定 ---
settings.modnam = 'clt_model';
settings.modnamstar = 'clt_model_zlb';
settings.constraint = 'i_H_notional < -i_H_ss';
settings.constraint_relax = 'i_H_notional > -i_H_ss';

settings.nperiods = 400;
settings.maxiter = 100;

%% --- シナリオ定義 (両ショックに共通の最適パラメータを適用) ---
scenarios = struct();

% シナリオ1: Beta Shock (需要ショック)
scenarios(1).name        = 'beta_shock';
scenarios(1).irfshock    = 'eps_beta_H';
scenarios(1).shock_val   = 0.02;
scenarios(1).phi_gap_H   = phi_gap_H_common;
scenarios(1).phi_level_H = phi_level_H_common;

% シナリオ2: A Shock (供給ショック)
scenarios(2).name        = 'a_shock';
scenarios(2).irfshock    = 'eps_a_H';
scenarios(2).shock_val   = -0.1;
scenarios(2).phi_gap_H   = phi_gap_H_common;
scenarios(2).phi_level_H = phi_level_H_common;

% サマリー記録用構造体
summary_welfares = struct();

%% --- ループ実行開始 ---
for s_idx = 1:length(scenarios)
    
    current_scen = scenarios(s_idx);
    phi_gap_H_value   = current_scen.phi_gap_H;
    phi_level_H_value = current_scen.phi_level_H;

    fprintf('\n--------------------------------------------------------------\n');
    fprintf('CLT シナリオ実行中 (%d/%d): %s (%s = %.2f)\n', ...
        s_idx, length(scenarios), current_scen.name, current_scen.irfshock, current_scen.shock_val);
    fprintf('使用パラメータ: phi_gap_H = %.6g, phi_level_H = %.8g\n', phi_gap_H_value, phi_level_H_value);
    fprintf('--------------------------------------------------------------\n');

    % 設定構造体の更新
    settings.shockssequence = current_scen.shock_val;
    settings.irfshock = current_scen.irfshock;
    
    % --- 出力先フォルダの設定 (results > [shock]) ---
    output_root_folder = fullfile(mpath, 'results', current_scen.name); 
    if ~exist(output_root_folder, 'dir'), mkdir(output_root_folder); end

    % --- 【核心】clt_solve_core.m の呼び出し ---
    [welfares, irf_results, Mbase_, oobase_] = clt_solve_core(phi_gap_H_value, phi_level_H_value, settings);

    % --- 固有値・収束・振動の診断チェック ---
    stable_eigs = oobase_.dr.eigval(abs(oobase_.dr.eigval) < 1 - 1e-5);
    if any(real(stable_eigs) < -1e-5)
        warning('[%s] Negative stable eigenvalue (lambda < 0) detected!', current_scen.name);
    end
    i_H_ss_val = Mbase_.params(strcmp('i_H_ss', Mbase_.param_names));
    zlb_actual   = (irf_results.i_H_pie <= -i_H_ss_val + 1e-6);
    zlb_notional = (irf_results.i_H_notional_pie < -i_H_ss_val - 1e-6);
    if any(zlb_actual ~= zlb_notional)
        warning('[%s] OccBin did not converge within maxiter (regime cycling)!', current_scen.name);
    end
    di_lin = diff(irf_results.i_H_lin);
    zigzag_count = sum((di_lin(1:end-1) .* di_lin(2:end) < 0) & (abs(di_lin(1:end-1)) > 1e-6));
    if zigzag_count > 10
        warning('[%s] Oscillating solution (zigzag chattering) detected!', current_scen.name);
    end

    % 取得した厚生変数を個別に展開
    welfare_H_with_delta_pie    = welfares.welfare_H_with_delta_pie;
    welfare_H_with_delta_lin    = welfares.welfare_H_with_delta_lin;
    welfare_H_without_delta_pie = welfares.welfare_H_without_delta_pie;
    welfare_H_without_delta_lin = welfares.welfare_H_without_delta_lin;
    welfare_F_with_delta_pie    = welfares.welfare_F_with_delta_pie;
    welfare_F_with_delta_lin    = welfares.welfare_F_with_delta_lin;
    welfare_F_without_delta_pie = welfares.welfare_F_without_delta_pie;
    welfare_F_without_delta_lin = welfares.welfare_F_without_delta_lin;

    % サマリー用に記録
    summary_welfares.(current_scen.name) = welfares;

    % irf_results構造体の中身を現在のワークスペースに展開
    fnames = fieldnames(irf_results);
    for k = 1:length(fnames)
        eval([fnames{k}, ' = irf_results.', fnames{k}, ';']);
    end

    %% --- グラフの描画 ---
    fprintf('グラフを作成し、PDFに保存中 (%s)...\n', current_scen.name);

    var_list = { ...
        'e_slash_star', ...
        'lambda_H', ...
        'lambda_H_slash_star', ...
        'lambda_F_slash_star', ...
        'lambda_F', ...
        'utility_H_with_delta', 'utility_F_with_delta', ...
        'utility_H_without_delta', 'utility_F_without_delta', ...
        'gamma_H', 'chi_H', ...
        'p_H', 'p_F_star', ...
        'i_H_notional', ...
        'i_H', 'i_F', ...
        'c_H_W', 'c_F_W', ...
        'c_H_H', 'c_F_H', ...
        'c_H_F', 'c_F_F', ...
        'y_H', 'y_F', ...
        'l_H', 'l_F', ...
        'p_H_W', 'p_F_W_star', ...
        'p_H_W_bar', 'p_F_W_star_bar', ...
        'p_H_bar', 'p_F_star_bar', ...
        'p_H_bar_y_H', 'p_H_W_c_H_W', 'p_H_W_slash_p_H', ...
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

    graphics_toolkit("gnuplot");
    pdf_filename = sprintf('%s_graphs_%s.pdf', policy_key, current_scen.name);
    output_filename = fullfile(output_root_folder, pdf_filename);
    if isfile(output_filename), delete(output_filename); end

    plots_per_page = 6;
    rows = 2;
    cols = 3;
    num_vars = length(var_list);
    nperiods = settings.nperiods;
    irfshock = settings.irfshock;
    shock_plot_vec = [settings.shockssequence(:); zeros(nperiods - length(settings.shockssequence(:)), 1)];

    % ページ上部に表示する共通ヘッダー文字列
    if strcmp(current_scen.name, 'beta_shock')
        scen_label_tex = '\beta shock';
    elseif strcmp(current_scen.name, 'a_shock')
        scen_label_tex = 'a shock';
    else
        scen_label_tex = strrep(current_scen.name, '_', ' ');
    end
    page_header_str = sprintf('%s (%s):   \\phi_{gap}^H = %.6g,   \\phi_{level}^H = %.8g', ...
        upper(policy_key), scen_label_tex, phi_gap_H_value, phi_level_H_value);

    for i = 1:plots_per_page:num_vars
        h = figure('Visible', 'off', 'PaperType', 'a4', 'PaperOrientation', 'landscape', 'PaperUnits', 'centimeters');
        
        graph_width = 25;
        graph_height = 17.5;
        left_margin = (29.7 - graph_width) / 2;
        bottom_margin = (21.0 - graph_height) / 2 - 0.4;
        set(h, 'PaperPosition', [left_margin, bottom_margin, graph_width, graph_height]);
        
        ax_list = [];
        j_list  = [];

        for j = 1:plots_per_page
            plot_index = i + j - 1;
            if plot_index > num_vars, break; end
            ax = subplot(rows, cols, j);
            ax_list(end+1) = ax;
            j_list(end+1)  = j;
            
            var_name = var_list{plot_index};
            
            if strncmp(var_name, 'eps_', 4)
                if strcmp(var_name, strtrim(irfshock))
                    plot(1:nperiods, shock_plot_vec, 'b-', 'LineWidth', 1.5);
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
            
            % --- タイトル生成ロジック (_{bar}, _{tilde} で確実に全文字を下付き化) ---
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
                latex_title = 'p_{bar}^{H} y^{H}';
            elseif strcmp(var_name, 'p_H_W_c_H_W') 
                latex_title = 'p^{H \rightarrow W} c^{H \rightarrow W}';
            elseif strcmp(var_name, 'p_H_W_slash_p_H') 
                latex_title = 'p^{H \rightarrow W} / p^{H}';
            elseif strcmp(var_name, 'pi_H_W')
                latex_title = '\pi^{H \rightarrow W}';
            elseif strcmp(var_name, 'pi_F_W_star')
                latex_title = '\pi^{F \rightarrow W*}';
            elseif strcmp(var_name, 'c_H_W')
                latex_title = 'c^{H \rightarrow W}';
            elseif strcmp(var_name, 'c_F_W')
                latex_title = 'c^{F \rightarrow W*}';
            elseif strcmp(var_name, 'p_H_W')
                latex_title = 'p^{H \rightarrow W}';
            elseif strcmp(var_name, 'p_F_W_star')
                latex_title = 'p^{F \rightarrow W*}';
            elseif strcmp(var_name, 'p_H_W_bar')
                latex_title = 'p_{bar}^{H \rightarrow W}';
            elseif strcmp(var_name, 'p_F_W_star_bar')
                latex_title = 'p_{bar}^{F \rightarrow W*}';
            elseif strcmp(var_name, 'c_H_H')
                latex_title = 'c^{H \rightarrow H}';
            elseif strcmp(var_name, 'c_F_H')
                latex_title = 'c^{F \rightarrow H}';
            elseif strcmp(var_name, 'c_H_F')
                latex_title = 'c^{H \rightarrow F*}';
            elseif strcmp(var_name, 'c_F_F')
                latex_title = 'c^{F \rightarrow F*}';
            elseif strcmp(var_name, 'p_H_tilde')
                latex_title = 'p_{tilde}^{H}';
            elseif strcmp(var_name, 'p_F_star_tilde')
                latex_title = 'p_{tilde}^{F*}';
            elseif strcmp(var_name, 'i_H_notional')
                latex_title = 'i^{H,notional}';
            elseif endsWith(var_name, '_bar')
                % p_H_bar, p_F_star_bar, p_F_bar, p_H_star_bar
                temp_name = regexprep(var_name, '_bar$', '');
                parts = strsplit(temp_name, '_');
                base = parts{1}; 
                if length(parts) > 1
                    superscript = strjoin(parts(2:end), '_');
                    superscript = strrep(superscript, '_', '\rightarrow');
                    superscript = strrep(superscript, '\rightarrowstar', '*');
                    superscript = strrep(superscript, 'star', '*'); 
                    latex_title = [base, '_{bar}^{', superscript, '}'];
                else
                    latex_title = [base, '_{bar}'];
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
                    superscript = strrep(superscript, 'star', '*'); 
                    latex_title = [base '^{' superscript '}'];
                else
                    latex_title = base;
                end
            end
            
            set(ax, 'FontSize', 8);
            title(ax, latex_title, 'Interpreter', 'tex', 'FontSize', 10);

            if plot_index == 1
                legend('ZLB', 'Non-ZLB', 'Location', 'northeast');
            end
            
            % 厚生の表示
            if strcmp(var_name, 'utility_H_with_delta')
                xlim_vals = get(ax, 'XLim');
                ylim_vals = get(ax, 'YLim');
                x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
                y_pos_top    = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
                y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
                text_H_pie = sprintf('Welfare (with Delta, ZLB): %.4f', welfare_H_with_delta_pie);
                text_H_lin = sprintf('Welfare (with Delta, Non-ZLB): %.4f', welfare_H_with_delta_lin);
                text(x_pos, y_pos_bottom, text_H_lin, 'Color', 'red',  'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
                text(x_pos, y_pos_top,    text_H_pie, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            elseif strcmp(var_name, 'utility_H_without_delta')
                xlim_vals = get(ax, 'XLim');
                ylim_vals = get(ax, 'YLim');
                x_pos = xlim_vals(1) + 0.05 * (xlim_vals(2) - xlim_vals(1));
                y_pos_top    = ylim_vals(1) + 0.15 * (ylim_vals(2) - ylim_vals(1));
                y_pos_bottom = ylim_vals(1) + 0.05 * (ylim_vals(2) - ylim_vals(1));
                text_H_pie = sprintf('Welfare (without Delta, ZLB): %.4f', welfare_H_without_delta_pie);
                text_H_lin = sprintf('Welfare (without Delta, Non-ZLB): %.4f', welfare_H_without_delta_lin);
                text(x_pos, y_pos_bottom, text_H_lin, 'Color', 'red',  'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
                text(x_pos, y_pos_top,    text_H_pie, 'Color', 'blue', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
            end
        end

        % 全サブプロットを作成し終えた後に、一括で高さを少し縮めてヘッダーとの間隔を広げる
        for k = 1:length(ax_list)
            pos = get(ax_list(k), 'Position');
            if j_list(k) <= 3
                pos(2) = pos(2) - 0.025;
                pos(4) = pos(4) * 0.90;
            else
                pos(2) = pos(2) - 0.015;
                pos(4) = pos(4) * 0.90;
            end
            set(ax_list(k), 'Position', pos);
        end

        % ページ上部ヘッダーの描画
        annotation(h, 'textbox', [0.13, 0.955, 0.80, 0.04], ...
            'String', page_header_str, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'tex');
        
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

    vars_to_save = {};
    for i = 1:length(var_list)
        var_name = var_list{i};
        if ~strncmp(var_name, 'eps_', 4) 
            vars_to_save{end+1} = [var_name, '_pie']; 
            vars_to_save{end+1} = [var_name, '_lin']; 
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
    vars_to_save{end+1} = 'phi_gap_H_value';
    vars_to_save{end+1} = 'phi_level_H_value';

    mat_filename = sprintf('%s_data_for_plotting_%s.mat', policy_key, current_scen.name);
    full_save_path = fullfile(output_root_folder, mat_filename);

    save(full_save_path, vars_to_save{:});
    fprintf('データ保存完了: %s\n', full_save_path);

    %% --- CSVおよびテキストへの書き出し ---
    fprintf('CSVファイルへの書き出しを開始します (%s)...\n', current_scen.name);

    try
        S = load(full_save_path);
        header = {}; 
        data_matrix = []; 
        
        for i = 1:length(vars_to_save)
            var_name = vars_to_save{i};
            if ~strncmp(var_name, 'welfare_', 8) && ~strncmp(var_name, 'phi_', 4)
                header{end+1} = var_name;
                data_matrix = [data_matrix, S.(var_name)(:)]; 
            end
        end
        
        csv_filename = sprintf('%s_data_for_plotting_%s.csv', policy_key, current_scen.name);
        csv_save_path = fullfile(output_root_folder, csv_filename);

        fid = fopen(csv_save_path, 'w');
        if fid == -1, error('CSVファイルを開けませんでした'); end
        
        fprintf(fid, '%s', header{1}); 
        for i = 2:length(header), fprintf(fid, ',%s', header{i}); end
        fprintf(fid, '\n'); 
        fclose(fid);
        
        dlmwrite(csv_save_path, data_matrix, '-append', 'delimiter', ',', 'precision', '%.8g');
        
        % --- 厚生・使用パラメータデータのテキスト保存 ---
        txt_filename = sprintf('%s_data_for_plotting_welfare_%s.txt', policy_key, current_scen.name);
        scalar_save_path = fullfile(output_root_folder, txt_filename);
        fid_scalar = fopen(scalar_save_path, 'w');
        if fid_scalar ~= -1
            fprintf(fid_scalar, 'Welfare Calculations (%s):\n', current_scen.name);
            fprintf(fid_scalar, 'Parameters: phi_gap_H = %.6g, phi_level_H = %.8g (%s)\n\n', ...
                phi_gap_H_value, phi_level_H_value, param_source_str);
            
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

            % projects フォルダ直下の data_for_plotting フォルダへコピー
            projects_root = fileparts(mpath);
            central_data_folder = fullfile(projects_root, 'data_for_plotting');
            if ~exist(central_data_folder, 'dir')
                mkdir(central_data_folder);
            end
            
            copyfile(csv_save_path, fullfile(central_data_folder, csv_filename));
            copyfile(scalar_save_path, fullfile(central_data_folder, txt_filename));
            
            fprintf('一覧用フォルダにコピーしました: %s, %s\n', csv_filename, txt_filename);
        end
        
    catch ME
        fprintf('CSVへの書き出し中にエラーが発生しました:\n');
        disp(ME.message);
    end
    
    close all; 
end

%% --- 両ショック統合の実行サマリー表示 ---
w_H_beta_pie = summary_welfares.beta_shock.welfare_H_with_delta_pie;
w_H_a_pie    = summary_welfares.a_shock.welfare_H_with_delta_pie;
w_H_exp_pie  = 0.5 * (w_H_beta_pie + w_H_a_pie);

fprintf('\n==============================================================\n');
fprintf('  全てのシナリオ（beta_shock, a_shock）の実行が完了しました。\n');
fprintf('  [適用パラメータ] phi_gap_H = %.6g, phi_level_H = %.8g\n', phi_gap_H_common, phi_level_H_common);
fprintf('  [Home厚生 (Delta有り)] Beta: %.6f, A: %.6f => 期待厚生: %.6f\n', ...
    w_H_beta_pie, w_H_a_pie, w_H_exp_pie);
fprintf('==============================================================\n');