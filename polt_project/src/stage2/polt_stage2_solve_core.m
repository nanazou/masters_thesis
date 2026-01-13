function [welfares, irf_results, Mbase_, oobase_] = polt_stage2_solve_core(phi_gap_H_value, phi_level_H_value, settings)

% src/stage2/polt_stage2_solve_core.m
% =========================================================================
% Stage 2: Sticky Price Model Simulation with OccBin (ZLB consideration)
%
% 計算エンジンの全ロジック（シミュレーション・変数展開・厚生計算）を記述。
% =========================================================================

global M_ oo_ options_;

% --- [1] ディレクトリ管理 ---
% 呼び出し元のディレクトリ（src フォルダ）のパスを保存
root_dir = pwd;

% この関数ファイルが存在する場所（src/stage2/）のパスを取得
mpath = fileparts(mfilename('fullpath'));

% 自分自身の場所（stage2/）に移動
% これにより、Dynare のビルド補助ファイルが stage2 内に生成され、管理が容易になります
cd(mpath);

% --- [2] パラメータの橋渡し ---
% .mod ファイル内で参照される政策パラメータをベースワークスペースへ送り込む
assignin('base', 'phi_gap_H_value', phi_gap_H_value);
assignin('base', 'phi_level_H_value', phi_level_H_value);

% --- [3] OccBin ソルバーの実行 ---
% solve_one_constraint.m を呼び出してシミュレーションを実行
fprintf('OccBinソルバー (Stage 2) を実行中...\n');

[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_] = ...
    solve_one_constraint(settings.modnam, settings.modnamstar, ...
    settings.constraint, settings.constraint_relax, ...
    settings.shockssequence, settings.irfshock, settings.nperiods, settings.maxiter);

fprintf('ソルバーの実行が完了しました。\n');

% --- [4] 呼び出し元のディレクトリ（src）に戻る ---
cd(root_dir);

% --- [5] 結果の展開 (偏差として取得) ---
% OccBin の出力（zdatapiecewise, zdatalinear）を内生変数名に展開
for i = 1:Mbase_.endo_nbr
    % _pie: ZLB考慮ケース (Piecewise Linear)
    eval([Mbase_.endo_names{i,:}, '_pie = zdatapiecewise(:, i);']);
    % _lin: ZLB無視ケース (Pure Linear)
    eval([Mbase_.endo_names{i,:}, '_lin = zdatalinear(:, i);']);
end

% 描画・保存用マッピング
% settings.shockssequence の2列目が Stage 1 から引き継いだ y_potential パス
y_H_potential_pie = settings.shockssequence(:, 2);
y_H_potential_lin = settings.shockssequence(:, 2);

%% --- [6] パラメータと定常状態の一括取得 ---
% 計算に必要な定数情報をモデル構造体から抽出
N_val         = Mbase_.params(strcmp('N', Mbase_.param_names));
M_val         = Mbase_.params(strcmp('M', Mbase_.param_names));
phi_H_val     = Mbase_.params(strcmp('phi_H', Mbase_.param_names));
phi_F_val     = Mbase_.params(strcmp('phi_F', Mbase_.param_names));
beta_H_ss_val = Mbase_.params(strcmp('beta_H_ss', Mbase_.param_names));
beta_F_ss_val = Mbase_.params(strcmp('beta_F_ss', Mbase_.param_names));

% 各変数の定常状態（Steady State）の値
e_slash_star_ss_val = oobase_.dr.ys(strcmp('e_slash_star', Mbase_.endo_names));
c_H_W_ss_val         = oobase_.dr.ys(strcmp('c_H_W', Mbase_.endo_names));
l_H_ss_val           = oobase_.dr.ys(strcmp('l_H', Mbase_.endo_names));
c_F_W_ss_val         = oobase_.dr.ys(strcmp('c_F_W', Mbase_.endo_names));
l_F_ss_val           = oobase_.dr.ys(strcmp('l_F', Mbase_.endo_names));
b_H_ss_val           = oobase_.dr.ys(strcmp('b_H', Mbase_.endo_names)); 

%% --- [7] 追加変数の計算 (b_F_star の生成など) ---
% 水準（Level）データの復元（偏差＋定常状態）
e_slash_star_level_pie = e_slash_star_ss_val + e_slash_star_pie;
e_slash_star_level_lin = e_slash_star_ss_val + e_slash_star_lin;

b_H_level_pie = b_H_ss_val + b_H_pie;
b_H_level_lin = b_H_ss_val + b_H_lin;

% Home債券残高から Foreign債券残高を計算
b_F_star_level_pie = - (N_val/M_val) * (b_H_level_pie ./ e_slash_star_level_pie);
b_F_star_level_lin = - (N_val/M_val) * (b_H_level_lin ./ e_slash_star_level_lin);

b_F_star_ss = - (N_val/M_val) * (b_H_ss_val / e_slash_star_ss_val);

% 偏差データとしての b_F_star を再作成
b_F_star_pie = b_F_star_level_pie - b_F_star_ss;
b_F_star_lin = b_F_star_level_lin - b_F_star_ss;

%% --- [8] 効用 (Utility) と厚生 (Welfare) の計算 ---

% 1. 各国パラメータの準備
theta_H_val = Mbase_.params(strcmp('theta_H', Mbase_.param_names));
xi_H_val    = Mbase_.params(strcmp('xi_H', Mbase_.param_names));
theta_F_val = Mbase_.params(strcmp('theta_F', Mbase_.param_names));
xi_F_val    = Mbase_.params(strcmp('xi_F', Mbase_.param_names));

nperiods = settings.nperiods;

% 2. 8パターン（Home/Foreign × Pie/Lin × With/Without Delta）の初期化
Delta_H_pie = zeros(nperiods, 1); Delta_F_pie = zeros(nperiods, 1);
Delta_H_lin = zeros(nperiods, 1); Delta_F_lin = zeros(nperiods, 1);

utility_H_with_delta_pie = zeros(nperiods, 1); utility_F_with_delta_pie = zeros(nperiods, 1);
utility_H_without_delta_pie = zeros(nperiods, 1); utility_F_without_delta_pie = zeros(nperiods, 1);
utility_H_with_delta_lin = zeros(nperiods, 1); utility_F_with_delta_lin = zeros(nperiods, 1);
utility_H_without_delta_lin = zeros(nperiods, 1); utility_F_without_delta_lin = zeros(nperiods, 1);

% 3. 時系列ループ計算
for t = 1:nperiods
    if t == 1
        DH_lag_p = 0; DF_lag_p = 0; DH_lag_l = 0; DF_lag_l = 0;
    else
        DH_lag_p = Delta_H_pie(t-1); DF_lag_p = Delta_F_pie(t-1);
        DH_lag_l = Delta_H_lin(t-1); DF_lag_l = Delta_F_lin(t-1);
    end
    
    % 価格分散項 Delta の推移計算 (再帰式)
    Delta_H_pie(t) = xi_H_val * DH_lag_p + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_pie(t)^2);
    Delta_F_pie(t) = xi_F_val * DF_lag_p + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_pie(t)^2);
    Delta_H_lin(t) = xi_H_val * DH_lag_l + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_lin(t)^2);
    Delta_F_lin(t) = xi_F_val * DF_lag_l + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_lin(t)^2);

    % --- Pie (ZLB考慮) の効用計算 ---
    utility_H_with_delta_pie(t) = c_H_W_pie(t) - phi_H_val*(l_H_ss_val^2)*(y_H_pie(t)-a_H_pie(t)) - phi_H_val*(l_H_ss_val^2)*Delta_H_pie(t);
    utility_F_with_delta_pie(t) = c_F_W_pie(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t)) - phi_F_val*(l_F_ss_val^2)*Delta_F_pie(t);
    utility_H_without_delta_pie(t) = c_H_W_pie(t) - phi_H_val*(l_H_ss_val^2)*(y_H_pie(t)-a_H_pie(t));
    utility_F_without_delta_pie(t) = c_F_W_pie(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t));
    
    % --- Lin (線形近似) の効用計算 ---
    utility_H_with_delta_lin(t) = c_H_W_lin(t) - phi_H_val*(l_H_ss_val^2)*(y_H_lin(t)-a_H_lin(t)) - phi_H_val*(l_H_ss_val^2)*Delta_H_lin(t);
    utility_F_with_delta_lin(t) = c_F_W_lin(t) - phi_F_val*(l_F_ss_val^2)*(y_F_lin(t)-a_F_lin(t)) - phi_F_val*(l_F_ss_val^2)*Delta_F_lin(t);
    utility_H_without_delta_lin(t) = c_H_W_lin(t) - phi_H_val*(l_H_ss_val^2)*(y_H_lin(t)-a_H_lin(t));
    utility_F_without_delta_lin(t) = c_F_W_lin(t) - phi_F_val*(l_F_ss_val^2)*(y_F_lin(t)-a_F_lin(t));
end

% 4. 生涯厚生（生涯累積割引効用）の合算
disc_H = beta_H_ss_val.^(0:nperiods-1)'; 
disc_F = beta_F_ss_val.^(0:nperiods-1)';

welfares.welfare_H_with_delta_pie = sum(disc_H .* utility_H_with_delta_pie);
welfares.welfare_H_with_delta_lin = sum(disc_H .* utility_H_with_delta_lin);
welfares.welfare_H_without_delta_pie = sum(disc_H .* utility_H_without_delta_pie);
welfares.welfare_H_without_delta_lin = sum(disc_H .* utility_H_without_delta_lin);

welfares.welfare_F_with_delta_pie = sum(disc_F .* utility_F_with_delta_pie);
welfares.welfare_F_with_delta_lin = sum(disc_F .* utility_F_with_delta_lin);
welfares.welfare_F_without_delta_pie = sum(disc_F .* utility_F_without_delta_pie);
welfares.welfare_F_without_delta_lin = sum(disc_F .* utility_F_without_delta_lin);

% 全世界合計（Global Welfare）
welfares.welfare_SUM_with_delta_pie = welfares.welfare_H_with_delta_pie + welfares.welfare_F_with_delta_pie;

%% --- [9] その他の連結変数などの生成 ---
% 対数線形近似モデルにおける積・商の関係を和・差で表現
lambda_H_slash_star_pie = e_slash_star_pie + lambda_H_pie; 
lambda_H_slash_star_lin = e_slash_star_lin + lambda_H_lin;
lambda_F_pie = lambda_F_slash_star_pie - e_slash_star_pie; 
lambda_F_lin = lambda_F_slash_star_lin - e_slash_star_lin;

% 名目GDPおよび名目総消費の偏差
p_H_bar_y_H_pie = p_H_bar_pie + y_H_pie; 
p_H_bar_y_H_lin = p_H_bar_lin + y_H_lin;
p_H_W_c_H_W_pie = p_H_W_pie + c_H_W_pie;
p_H_W_c_H_W_lin = p_H_W_lin + c_H_W_lin; 

% --- [10] 戻り値の整理 (Packaging) ---
% who を使用して計算された全データを構造体にパッキングして返却
vars = [who('*_pie'); who('*_lin'); who('Delta_H_*'); who('Delta_F_*'); who('utility_*')];
for k = 1:length(vars)
    irf_results.(vars{k}) = eval(vars{k});
end

% 出力構造体の保持
Mbase_ = Mbase_;
oobase_ = oobase_;

end