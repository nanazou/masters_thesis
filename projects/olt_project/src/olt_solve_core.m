function [welfares, irf_results, Mbase_, oobase_] = olt_solve_core(phi_gap_H_value, phi_level_H_value, settings)

% src/olt_solve_core.m
% =========================================================================
% Solver実行、変数展開、Welfare計算 (OLT用)
% =========================================================================

global M_ oo_;

% --- カレントディレクトリの保存と移動 ---
% .modファイルがある src フォルダへ一時的に移動します
root_dir = pwd;
mpath = fileparts(mfilename('fullpath'));
cd(mpath);

% --- パラメータの橋渡し ---
assignin('base', 'phi_gap_H_value', phi_gap_H_value);
assignin('base', 'phi_level_H_value', phi_level_H_value);

% --- シミュレーション設定の展開 ---
nperiods = settings.nperiods;
maxiter = settings.maxiter;
modnam = settings.modnam;         % 'olt_model'
modnamstar = settings.modnamstar; % 'olt_model_zlb'
constraint = settings.constraint;
constraint_relax = settings.constraint_relax;
shockssequence = settings.shockssequence;
irfshock = settings.irfshock;

% --- OccBinソルバーの実行 ---
fprintf('OccBinソルバーを実行中...\n');
[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_] = ...
    solve_one_constraint(modnam, modnamstar, ...
    constraint, constraint_relax, ...
    shockssequence, irfshock, nperiods, maxiter);
fprintf('ソルバーの実行が完了しました。\n');

% 元のフォルダに戻る
cd(root_dir);

% 結果をワークスペースの変数に展開
% (前のループの変数が混ざらないよう上書き)
for i=1:Mbase_.endo_nbr
    eval([Mbase_.endo_names{i,:},'_pie=zdatapiecewise(:,i);']);
    eval([Mbase_.endo_names{i,:},'_lin=zdatalinear(:,i);']);
end

%% --- パラメータと定常状態の一括取得 ---
% 1. パラメータの取得
N_val         = Mbase_.params(strcmp('N', Mbase_.param_names));
M_val         = Mbase_.params(strcmp('M', Mbase_.param_names));
phi_H_val     = Mbase_.params(strcmp('phi_H', Mbase_.param_names));
phi_F_val     = Mbase_.params(strcmp('phi_F', Mbase_.param_names));
beta_H_ss_val = Mbase_.params(strcmp('beta_H_ss', Mbase_.param_names));
beta_F_ss_val = Mbase_.params(strcmp('beta_F_ss', Mbase_.param_names));

% 2. 定常状態（Steady State）の取得
e_slash_star_ss_val = oobase_.dr.ys(strcmp('e_slash_star', Mbase_.endo_names));
c_H_W_ss_val         = oobase_.dr.ys(strcmp('c_H_W', Mbase_.endo_names));
l_H_ss_val           = oobase_.dr.ys(strcmp('l_H', Mbase_.endo_names));
c_F_W_ss_val         = oobase_.dr.ys(strcmp('c_F_W', Mbase_.endo_names));
l_F_ss_val           = oobase_.dr.ys(strcmp('l_F', Mbase_.endo_names));

%% --- 追加変数の計算 ---

% ▼▼▼ b_F_star の生成 (水準経由で計算) ▼▼▼
e_slash_star_level_pie = e_slash_star_ss_val + e_slash_star_pie;
e_slash_star_level_lin = e_slash_star_ss_val + e_slash_star_lin;

b_H_level_pie = b_H_pie;
b_H_level_lin = b_H_lin;

b_F_star_level_pie = - (N_val/M_val) * (b_H_level_pie ./ e_slash_star_level_pie);
b_F_star_level_lin = - (N_val/M_val) * (b_H_level_lin ./ e_slash_star_level_lin);

b_F_star_pie = b_F_star_level_pie;
b_F_star_lin = b_F_star_level_lin;

% ▼▼▼ 効用 (Utility) と厚生 (Welfare) の計算 ▼▼▼

% 1. パラメータの準備
theta_H_val = Mbase_.params(strcmp('theta_H', Mbase_.param_names));
xi_H_val    = Mbase_.params(strcmp('xi_H', Mbase_.param_names));
theta_F_val = Mbase_.params(strcmp('theta_F', Mbase_.param_names));
xi_F_val    = Mbase_.params(strcmp('xi_F', Mbase_.param_names));

% 2. 価格分散 Delta_H, Delta_F の復元と効用の計算 (ZLBケース: _pie)
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

    % 価格分散 Delta の2次近似計算:
    % $$\Delta_t = \xi \Delta_{t-1} + 0.5 \frac{\theta \xi}{1 - \xi} \pi_t^2$$
    Delta_H_pie(t) = xi_H_val * Delta_H_lag + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_val_t^2);
    Delta_F_pie(t) = xi_F_val * Delta_F_lag + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_val_t^2);

    % 効用の計算 (y - a の項は betaショック/aショック両対応)
    term_labor_H = (y_H_pie(t) - a_H_pie(t)); 
    term_labor_F = (y_F_pie(t) - a_F_pie(t));
    
    % パターンA: Deltaあり (with_delta)
    utility_H_with_delta_pie(t) = c_H_W_pie(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_pie(t);

    utility_F_with_delta_pie(t) = c_F_W_pie(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_pie(t);

    % パターンB: Deltaなし (without_delta)
    utility_H_without_delta_pie(t) = c_H_W_pie(t) ...
                                   - phi_H_val * (l_H_ss_val^2) * term_labor_H;

    utility_F_without_delta_pie(t) = c_F_W_pie(t) ...
                                   - phi_F_val * (l_F_ss_val^2) * term_labor_F;
end

% 3. 価格分散 Delta_H, Delta_F の復元と効用の計算 (Non-ZLBケース: _lin)
Delta_H_lin = zeros(nperiods, 1);
Delta_F_lin = zeros(nperiods, 1);
utility_H_with_delta_lin = zeros(nperiods, 1);
utility_F_with_delta_lin = zeros(nperiods, 1);
utility_H_without_delta_lin = zeros(nperiods, 1);
utility_F_without_delta_lin = zeros(nperiods, 1);

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

    % パターンA: Deltaあり (with_delta)
    utility_H_with_delta_lin(t) = c_H_W_lin(t) ...
                         - phi_H_val * (l_H_ss_val^2) * term_labor_H ...
                         - phi_H_val * (l_H_ss_val^2) * Delta_H_lin(t);

    utility_F_with_delta_lin(t) = c_F_W_lin(t) ...
                         - phi_F_val * (l_F_ss_val^2) * term_labor_F ...
                         - phi_F_val * (l_F_ss_val^2) * Delta_F_lin(t);

    % パターンB: Deltaなし (without_delta)
    utility_H_without_delta_lin(t) = c_H_W_lin(t) ...
                                   - phi_H_val * (l_H_ss_val^2) * term_labor_H;

    utility_F_without_delta_lin(t) = c_F_W_lin(t) ...
                                   - phi_F_val * (l_F_ss_val^2) * term_labor_F;
end

% --- 生涯効用（厚生）の計算 ---
discounts_H = beta_H_ss_val.^(0:nperiods-1)';
welfare_H_with_delta_pie = sum(discounts_H .* utility_H_with_delta_pie);
welfare_H_with_delta_lin = sum(discounts_H .* utility_H_with_delta_lin);
welfare_H_without_delta_pie = sum(discounts_H .* utility_H_without_delta_pie);
welfare_H_without_delta_lin = sum(discounts_H .* utility_H_without_delta_lin);

discounts_F = beta_F_ss_val.^(0:nperiods-1)';
welfare_F_with_delta_pie = sum(discounts_F .* utility_F_with_delta_pie);
welfare_F_with_delta_lin = sum(discounts_F .* utility_F_with_delta_lin);
welfare_F_without_delta_pie = sum(discounts_F .* utility_F_without_delta_pie);
welfare_F_without_delta_lin = sum(discounts_F .* utility_F_without_delta_lin);

welfares.welfare_H_with_delta_pie = welfare_H_with_delta_pie;
welfares.welfare_H_with_delta_lin = welfare_H_with_delta_lin;
welfares.welfare_H_without_delta_pie = welfare_H_without_delta_pie;
welfares.welfare_H_without_delta_lin = welfare_H_without_delta_lin;
welfares.welfare_F_with_delta_pie = welfare_F_with_delta_pie;
welfares.welfare_F_with_delta_lin = welfare_F_with_delta_lin;
welfares.welfare_F_without_delta_pie = welfare_F_without_delta_pie;
welfares.welfare_F_without_delta_lin = welfare_F_without_delta_lin;
welfares.welfare_SUM_with_delta_pie = welfare_H_with_delta_pie + welfare_F_with_delta_pie;

%% --- その他の連結変数などの生成 ---
lambda_H_slash_star_pie = e_slash_star_pie + lambda_H_pie;
lambda_H_slash_star_lin = e_slash_star_lin + lambda_H_lin;
lambda_F_pie = lambda_H_slash_star_pie - e_slash_star_pie;
lambda_F_lin = lambda_H_slash_star_lin - e_slash_star_lin;

% --- 【新規追加: 変数名を積の形に変更】 ---
p_H_bar_y_H_pie = p_H_bar_pie + y_H_pie;
p_H_bar_y_H_lin = p_H_bar_lin + y_H_lin;
p_H_W_c_H_W_pie = p_H_W_pie + c_H_W_pie;
p_H_W_c_H_W_lin = p_H_W_lin + c_H_W_lin; 

% --- 戻り値の整理 ---
% who の引数を個別に取得して連結 (MATLAB/Octaveの構文エラー回避)
vars = [who('*_pie'); who('*_lin'); who('Delta_H_*'); who('Delta_F_*'); who('utility_*')];
for k = 1:length(vars)
    irf_results.(vars{k}) = eval(vars{k});
end

end