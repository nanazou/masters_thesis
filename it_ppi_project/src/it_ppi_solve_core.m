function [welfares, irf_results, Mbase_, oobase_] = it_ppi_solve_core(phi_pi_H_value, settings)

% src/it_ppi_solve_core.m
% =========================================================================
% OccBinソルバーの実行、内生変数の展開、および厚生（Welfare）の再計算
%
% 注意：このスクリプトは src フォルダ内に配置され、
% 同フォルダ内の it_ppi_model.mod および it_ppi_model_zlb.mod を読み込みます。
% =========================================================================

global M_ oo_;

% --- [1] カレントディレクトリの保存と移動 ---
% シミュレーション実行前にモデルファイルが存在するディレクトリへ移動します
root_dir = pwd;
mpath = fileparts(mfilename('fullpath'));
cd(mpath);

% --- [2] パラメータの橋渡し ---
% Dynareの定常状態計算で使用するため、ベース・ワークスペースへ値を渡します
assignin('base', 'phi_pi_H_value', phi_pi_H_value);

% --- [3] シミュレーション設定の展開 ---
nperiods = settings.nperiods;
maxiter = settings.maxiter;
modnam = settings.modnam;         % 'it_ppi_model' (命名規則統一済み)
modnamstar = settings.modnamstar; % 'it_ppi_model_zlb' (命名規則統一済み)
constraint = settings.constraint;
constraint_relax = settings.constraint_relax;
shockssequence = settings.shockssequence;
irfshock = settings.irfshock;

% --- [4] OccBinソルバーの実行 ---
fprintf('OccBinソルバーを実行中...\n');
[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_] = ...
    solve_one_constraint(modnam, modnamstar, ...
    constraint, constraint_relax, ...
    shockssequence, irfshock, nperiods, maxiter);
fprintf('ソルバーの実行が完了しました。\n');

% 呼び出し元のフォルダ（プロジェクトルート等）に戻る
cd(root_dir);

% --- [5] 結果をワークスペースの変数に展開 ---
% zdatalinear (線形) と zdatapiecewise (OccBin) を変数名ごとに分解します
% ここで pi_F_W_star 等も自動的に展開されます
for i=1:Mbase_.endo_nbr
    eval([Mbase_.endo_names{i,:},'_pie=zdatapiecewise(:,i);']);
    eval([Mbase_.endo_names{i,:},'_lin=zdatalinear(:,i);']);
end

%% --- [6] パラメータと定常状態の一括取得 ---
% 厚生計算に使用するため、定常状態の値とパラメータを抽出します
N_val         = Mbase_.params(strcmp('N', Mbase_.param_names));
M_val         = Mbase_.params(strcmp('M', Mbase_.param_names));
phi_H_val     = Mbase_.params(strcmp('phi_H', Mbase_.param_names));
phi_F_val     = Mbase_.params(strcmp('phi_F', Mbase_.param_names));
beta_H_ss_val = Mbase_.params(strcmp('beta_H_ss', Mbase_.param_names));
beta_F_ss_val = Mbase_.params(strcmp('beta_F_ss', Mbase_.param_names));

% 各変数の定常状態水準値
e_slash_star_ss_val = oobase_.dr.ys(strcmp('e_slash_star', Mbase_.endo_names));
c_H_W_ss_val         = oobase_.dr.ys(strcmp('c_H_W', Mbase_.endo_names));
l_H_ss_val           = oobase_.dr.ys(strcmp('l_H', Mbase_.endo_names));
c_F_W_ss_val         = oobase_.dr.ys(strcmp('c_F_W', Mbase_.endo_names));
l_F_ss_val           = oobase_.dr.ys(strcmp('l_F', Mbase_.endo_names));

%% --- [7] 追加変数の計算 ---

% ▼▼▼ b_F_star の生成 (為替レートを考慮した水準ベースの計算) ▼▼▼
% 資産残高の評価には実効為替レートの水準が必要なため復元します
e_slash_star_level_pie = e_slash_star_ss_val + e_slash_star_pie;
e_slash_star_level_lin = e_slash_star_ss_val + e_slash_star_lin;

b_H_level_pie = b_H_pie;
b_H_level_lin = b_H_lin;

% 国際予算制約に基づき自国債券から外国債券残高を算出
b_F_star_level_pie = - (N_val/M_val) * (b_H_level_pie ./ e_slash_star_level_pie);
b_F_star_level_lin = - (N_val/M_val) * (b_H_level_lin ./ e_slash_star_level_lin);

b_F_star_pie = b_F_star_level_pie;
b_F_star_lin = b_F_star_level_lin;

%% --- [8] 効用 (Utility) と厚生 (Welfare) の計算 ---

% 1. パラメータの準備
theta_H_val = Mbase_.params(strcmp('theta_H', Mbase_.param_names));
xi_H_val    = Mbase_.params(strcmp('xi_H', Mbase_.param_names));
theta_F_val = Mbase_.params(strcmp('theta_F', Mbase_.param_names));
xi_F_val    = Mbase_.params(strcmp('xi_F', Mbase_.param_names));

% 2. 8パターン（Home/Foreign × Pie/Lin × With/Without Delta）の初期化
Delta_H_pie = zeros(nperiods, 1); Delta_F_pie = zeros(nperiods, 1);
utility_H_with_delta_pie = zeros(nperiods, 1); utility_F_with_delta_pie = zeros(nperiods, 1);
utility_H_without_delta_pie = zeros(nperiods, 1); utility_F_without_delta_pie = zeros(nperiods, 1);

% --- ZLB考慮ケース (_pie) のループ計算 ---
for t = 1:nperiods
    if t == 1
        Delta_H_lag = 0; Delta_F_lag = 0;
    else
        Delta_H_lag = Delta_H_pie(t-1); Delta_F_lag = Delta_F_pie(t-1);
    end
    
    % インフレ率の取得
    pi_H_val_t = pi_H_pie(t); 
    pi_F_star_val_t = pi_F_star_pie(t);

    % 価格分散 Delta の2次近似計算: Δ_t = ξΔ_{t-1} + 0.5 * (θξ/(1-ξ)) * π_t^2
    Delta_H_pie(t) = xi_H_val * Delta_H_lag + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_val_t^2);
    Delta_F_pie(t) = xi_F_val * Delta_F_lag + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_val_t^2);

    % 労働供給の歪み項 (y - a)
    term_labor_H = (y_H_pie(t) - a_H_pie(t)); 
    term_labor_F = (y_F_pie(t) - a_F_pie(t));
    
    % パターンA: Deltaあり (with_delta) - 価格粘着性による資源配分の歪みを考慮
    utility_H_with_delta_pie(t) = c_H_W_pie(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_pie(t);

    utility_F_with_delta_pie(t) = c_F_W_pie(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_pie(t);

    % パターンB: Deltaなし (without_delta) - 単純な実質変数のみの効用
    utility_H_without_delta_pie(t) = c_H_W_pie(t) ...
                                   - phi_H_val * (l_H_ss_val^2) * term_labor_H;

    utility_F_without_delta_pie(t) = c_F_W_pie(t) ...
                                   - phi_F_val * (l_F_ss_val^2) * term_labor_F;
end

% 3. 非ZLBケース (_lin) での計算
Delta_H_lin = zeros(nperiods, 1); Delta_F_lin = zeros(nperiods, 1);
utility_H_with_delta_lin = zeros(nperiods, 1); utility_F_with_delta_lin = zeros(nperiods, 1);
utility_H_without_delta_lin = zeros(nperiods, 1); utility_F_without_delta_lin = zeros(nperiods, 1);

for t = 1:nperiods
    if t == 1
        Delta_H_lag = 0; Delta_F_lag = 0;
    else
        Delta_H_lag = Delta_H_lin(t-1); Delta_F_lag = Delta_F_lin(t-1);
    end
    
    pi_H_val_t = pi_H_lin(t);
    pi_F_star_val_t = pi_F_star_lin(t);

    Delta_H_lin(t) = xi_H_val * Delta_H_lag + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_val_t^2);
    Delta_F_lin(t) = xi_F_val * Delta_F_lag + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_val_t^2);

    term_labor_H = (y_H_lin(t) - a_H_lin(t));
    term_labor_F = (y_F_lin(t) - a_F_lin(t));

    utility_H_with_delta_lin(t) = c_H_W_lin(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_lin(t);

    utility_F_with_delta_lin(t) = c_F_W_lin(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_lin(t);

    utility_H_without_delta_lin(t) = c_H_W_lin(t) ...
                                   - phi_H_val * (l_H_ss_val^2) * term_labor_H;

    utility_F_without_delta_lin(t) = c_F_W_lin(t) ...
                                   - phi_F_val * (l_F_ss_val^2) * term_labor_F;
end

% --- [9] 生涯効用（厚生）の合計計算 ---
% 定常状態の割引因子 beta で割り引いた累計を算出します
discounts_H = beta_H_ss_val.^(0:nperiods-1)';
discounts_F = beta_F_ss_val.^(0:nperiods-1)';

% 最終結果の格納
welfares.welfare_H_with_delta_pie = sum(discounts_H .* utility_H_with_delta_pie);
welfares.welfare_H_with_delta_lin = sum(discounts_H .* utility_H_with_delta_lin);
welfares.welfare_H_without_delta_pie = sum(discounts_H .* utility_H_without_delta_pie);
welfares.welfare_H_without_delta_lin = sum(discounts_H .* utility_H_without_delta_lin);

welfares.welfare_F_with_delta_pie = sum(discounts_F .* utility_F_with_delta_pie);
welfares.welfare_F_with_delta_lin = sum(discounts_F .* utility_F_with_delta_lin);
welfares.welfare_F_without_delta_pie = sum(discounts_F .* utility_F_without_delta_pie);
welfares.welfare_F_without_delta_lin = sum(discounts_F .* utility_F_without_delta_lin);

% 全世界合計（Global Welfare）
welfares.welfare_SUM_with_delta_pie = welfares.welfare_H_with_delta_pie + welfares.welfare_F_with_delta_pie;

%% --- [10] その他の連結変数の生成 ---
% 対数線形近似における積・商の関係（和・差）を計算します
lambda_H_slash_star_pie = e_slash_star_pie + lambda_H_pie;
lambda_H_slash_star_lin = e_slash_star_lin + lambda_H_lin;
lambda_F_pie = lambda_H_slash_star_pie - e_slash_star_pie;
lambda_F_lin = lambda_H_slash_star_lin - e_slash_star_lin;

% 名目GDPおよび名目総消費
p_H_bar_y_H_pie = p_H_bar_pie + y_H_pie;
p_H_bar_y_H_lin = p_H_bar_lin + y_H_lin;
p_H_W_c_H_W_pie = p_H_W_pie + c_H_W_pie;
p_H_W_c_H_W_lin = p_H_W_lin + c_H_W_lin; 

% --- [11] 戻り値の整理 ---
% whoコマンドを用いて計算結果をirf_results構造体にまとめます
vars = [who('*_pie'); who('*_lin'); who('Delta_H_*'); who('Delta_F_*'); who('utility_*')];
for k = 1:length(vars)
    irf_results.(vars{k}) = eval(vars{k});
end

end