function [welfares, irf_results, Mbase_, oobase_] = tr_ppi_stage2_solve_core(phi_pi_H_value, phi_y_H_value, settings)

% src/stage2/tr_ppi_stage2_solve_core.m
% =========================================================================
% Solver実行、変数展開、Welfare計算 (TR PPI用 - 潜在産出量対応版)
% 
% 目的:
%   1. OccBinを使用して、ZLBを考慮したテイラー・ルール政策をシミュレートする。
%   2. 潜在産出量（Stage 1結果）を用いた本来の出力ギャップに基づく政策を評価。
%   3. 価格分散（Delta）を含む厳密な厚生損失を算出する。
% =========================================================================

global M_ oo_;

% --- [1] カレントディレクトリの保存と移動 ---
% .modファイルやビルドファイルが生成される src/stage2 フォルダへ一時的に移動します。
root_dir = pwd;
mpath = fileparts(mfilename('fullpath'));
cd(mpath);

% --- [2] パラメータの橋渡し ---
% 最適化スクリプトやメインスクリプトから渡された係数をDynareワークスペースへ。
assignin('base', 'phi_pi_H_value', phi_pi_H_value);
assignin('base', 'phi_y_H_value', phi_y_H_value);

% --- [3] シミュレーション設定の展開 ---
nperiods = settings.nperiods;
maxiter = settings.maxiter;
modnam = settings.modnam;         % 'tr_ppi_stage2_model'
modnamstar = settings.modnamstar; % 'tr_ppi_stage2_model_zlb'
constraint = settings.constraint;
constraint_relax = settings.constraint_relax;
shockssequence = settings.shockssequence; % 1列目: 主ショック, 2列目: 潜在産出量
irfshock = settings.irfshock;

% --- [4] OccBinソルバーの実行 ---
fprintf('OccBinソルバー (TR PPI Stage 2) を実行中...\n');
[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_] = ...
    solve_one_constraint(modnam, modnamstar, ...
    constraint, constraint_relax, ...
    shockssequence, irfshock, nperiods, maxiter);
fprintf('ソルバーの実行が完了しました。\n');

% 元のフォルダ (src) に戻る
cd(root_dir);

% --- [5] 結果をワークスペースの変数に展開 ---
% OccBinの出力を変数名ごとに展開し、ZLB考慮(_pie)と考慮なし(_lin)の両方を保持。
for i=1:Mbase_.endo_nbr
    eval([Mbase_.endo_names{i,:},'_pie=zdatapiecewise(:,i);']);
    eval([Mbase_.endo_names{i,:},'_lin=zdatalinear(:,i);']);
end

%% --- [6] パラメータと定常状態の一括取得 ---
N_val         = Mbase_.params(strcmp('N', Mbase_.param_names));
M_val         = Mbase_.params(strcmp('M', Mbase_.param_names));
phi_H_val     = Mbase_.params(strcmp('phi_H', Mbase_.param_names));
phi_F_val     = Mbase_.params(strcmp('phi_F', Mbase_.param_names));
beta_H_ss_val = Mbase_.params(strcmp('beta_H_ss', Mbase_.param_names));
beta_F_ss_val = Mbase_.params(strcmp('beta_F_ss', Mbase_.param_names));

% 各変数の定常状態（Steady State）の取得
e_slash_star_ss_val = oobase_.dr.ys(strcmp('e_slash_star', Mbase_.endo_names)); 
c_H_W_ss_val         = oobase_.dr.ys(strcmp('c_H_W', Mbase_.endo_names));
l_H_ss_val           = oobase_.dr.ys(strcmp('l_H', Mbase_.endo_names));
c_F_W_ss_val         = oobase_.dr.ys(strcmp('c_F_W', Mbase_.endo_names));
l_F_ss_val           = oobase_.dr.ys(strcmp('l_F', Mbase_.endo_names));

%% --- [7] 追加変数の計算 ---

% ▼▼▼ b_F_star の生成 (e_slash_star を使用して水準値を復元) ▼▼▼
e_slash_star_level_pie = e_slash_star_ss_val + e_slash_star_pie;
e_slash_star_level_lin = e_slash_star_ss_val + e_slash_star_lin;

b_H_level_pie = b_H_pie;
b_H_level_lin = b_H_lin;

b_F_star_level_pie = - (N_val/M_val) * (b_H_level_pie ./ e_slash_star_level_pie);
b_F_star_level_lin = - (N_val/M_val) * (b_H_level_lin ./ e_slash_star_level_lin);

b_F_star_pie = b_F_star_level_pie;
b_F_star_lin = b_F_star_level_lin;

% ▼▼▼ 効用 (Utility) と厚生 (Welfare) の計算 ▼▼▼

% 1. 各国の粘着性パラメータ等の準備
theta_H_val = Mbase_.params(strcmp('theta_H', Mbase_.param_names));
xi_H_val    = Mbase_.params(strcmp('xi_H', Mbase_.param_names));
theta_F_val = Mbase_.params(strcmp('theta_F', Mbase_.param_names));
xi_F_val    = Mbase_.params(strcmp('xi_F', Mbase_.param_names));

% 2. 8パターンの計算用配列の初期化 (ZLBケース: _pie)
Delta_H_pie = zeros(nperiods, 1);
Delta_F_pie = zeros(nperiods, 1);
utility_H_with_delta_pie = zeros(nperiods, 1);
utility_F_with_delta_pie = zeros(nperiods, 1);
utility_H_without_delta_pie = zeros(nperiods, 1);
utility_F_without_delta_pie = zeros(nperiods, 1);

for t = 1:nperiods
    if t == 1
        Delta_H_lag = 0; Delta_F_lag = 0;
    else
        Delta_H_lag = Delta_H_pie(t-1); Delta_F_lag = Delta_F_pie(t-1);
    end
    
    pi_H_val_t = pi_H_pie(t); 
    pi_F_star_val_t = pi_F_star_pie(t);

    % 価格分散 Delta の動学的更新 (二次の損失項)
    Delta_H_pie(t) = xi_H_val * Delta_H_lag + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_val_t^2);
    Delta_F_pie(t) = xi_F_val * Delta_F_lag + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_val_t^2);

    % 各国の不効用（労働）項の準備
    term_labor_H = (y_H_pie(t) - a_H_pie(t)); 
    term_labor_F = (y_F_pie(t) - a_F_pie(t));
    
    % パターンA: Deltaあり (価格分散による厚生損失を考慮)
    utility_H_with_delta_pie(t) = c_H_W_pie(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_pie(t);

    utility_F_with_delta_pie(t) = c_F_W_pie(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_pie(t);

    % パターンB: Deltaなし (純粋な消費と労働の効用)
    utility_H_without_delta_pie(t) = c_H_W_pie(t) ...
                                    - phi_H_val * (l_H_ss_val^2) * term_labor_H;

    utility_F_without_delta_pie(t) = c_F_W_pie(t) ...
                                    - phi_F_val * (l_F_ss_val^2) * term_labor_F;
end

% 3. 同様の計算を Non-ZLBケース (_lin) に対しても完全に実行
Delta_H_lin = zeros(nperiods, 1);
Delta_F_lin = zeros(nperiods, 1);
utility_H_with_delta_lin = zeros(nperiods, 1);
utility_F_with_delta_lin = zeros(nperiods, 1);
utility_H_without_delta_lin = zeros(nperiods, 1);
utility_F_without_delta_lin = zeros(nperiods, 1);

for t = 1:nperiods
    if t == 1
        Delta_H_lag_l = 0; Delta_F_lag_l = 0;
    else
        Delta_H_lag_l = Delta_H_lin(t-1); Delta_F_lag_l = Delta_F_lin(t-1);
    end
    
    pi_H_val_t_l = pi_H_lin(t);
    pi_F_star_val_t_l = pi_F_star_lin(t);

    Delta_H_lin(t) = xi_H_val * Delta_H_lag_l + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_val_t_l^2);
    Delta_F_lin(t) = xi_F_val * Delta_F_lag_l + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_val_t_l^2);

    term_labor_H_l = (y_H_lin(t) - a_H_lin(t));
    term_labor_F_l = (y_F_lin(t) - a_F_lin(t));

    utility_H_with_delta_lin(t) = c_H_W_lin(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H_l ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_lin(t);

    utility_F_with_delta_lin(t) = c_F_W_lin(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F_l ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_lin(t);

    utility_H_without_delta_lin(t) = c_H_W_lin(t) ...
                                    - phi_H_val * (l_H_ss_val^2) * term_labor_H_l;

    utility_F_without_delta_lin(t) = c_F_W_lin(t) ...
                                    - phi_F_val * (l_F_ss_val^2) * term_labor_F_l;
end

% --- [9] 生涯厚生（累積割引効用）の計算 ---
discounts_H = beta_H_ss_val.^(0:nperiods-1)';
discounts_F = beta_F_ss_val.^(0:nperiods-1)';

% ZLB考慮
welfares.welfare_H_with_delta_pie = sum(discounts_H .* utility_H_with_delta_pie);
welfares.welfare_H_without_delta_pie = sum(discounts_H .* utility_H_without_delta_pie);
welfares.welfare_F_with_delta_pie = sum(discounts_F .* utility_F_with_delta_pie);
welfares.welfare_F_without_delta_pie = sum(discounts_F .* utility_F_without_delta_pie);

% ZLB無視
welfares.welfare_H_with_delta_lin = sum(discounts_H .* utility_H_with_delta_lin);
welfares.welfare_H_without_delta_lin = sum(discounts_H .* utility_H_without_delta_lin);
welfares.welfare_F_with_delta_lin = sum(discounts_F .* utility_F_with_delta_lin);
welfares.welfare_F_without_delta_lin = sum(discounts_F .* utility_F_without_delta_lin);

% 合計
welfares.welfare_SUM_with_delta_pie = welfares.welfare_H_with_delta_pie + welfares.welfare_F_with_delta_pie;

%% --- [10] そのめる連結変数などの生成 ---
lambda_H_slash_star_pie = e_slash_star_pie + lambda_H_pie;
lambda_H_slash_star_lin = e_slash_star_lin + lambda_H_lin;
lambda_F_pie = lambda_H_slash_star_pie - e_slash_star_pie;
lambda_F_lin = lambda_F_slash_star_lin - e_slash_star_lin;

% 名目GDPおよび名目総消費の積（近似的な和）
p_H_bar_y_H_pie = p_H_bar_pie + y_H_pie;
p_H_bar_y_H_lin = p_H_bar_lin + y_H_lin;
p_H_W_c_H_W_pie = p_H_W_pie + c_H_W_pie;
p_H_W_c_H_W_lin = p_H_W_lin + c_H_W_lin; 

% --- [11] 戻り値の整理 (Packaging) ---
% who の引数を個別に取得して連結
vars = [who('*_pie'); who('*_lin'); who('Delta_H_*'); who('Delta_F_*'); who('utility_*')];
for k = 1:length(vars)
    irf_results.(vars{k}) = eval(vars{k});
end

end