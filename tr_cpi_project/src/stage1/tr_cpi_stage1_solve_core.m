function [welfares, irf_results, Mbase_, oobase_] = tr_cpi_stage1_solve_core(phi_pi_H_value, phi_pi_F_value, settings)

% src/stage1/tr_cpi_stage1_solve_core.m
% =========================================================================
% Stage 1: Flexible Price Model Simulation (Linear) - TR CPI Project
%
% この関数は、伸縮価格モデル（Stage 1）を実行し、
% 自然水準（Potential Output）のパスと厚生（Welfare）を計算します。
% 実行時に一時的に src/stage1/ フォルダに移動し、終了時に呼び出し元に戻ります。
% =========================================================================

global M_ oo_ options_;

% --- [1] ディレクトリ管理 ---
% 呼び出し元のディレクトリ（src フォルダ）のパスを保存
root_dir = pwd;

% この関数ファイルが存在する場所（src/stage1/）のパスを取得
mpath = fileparts(mfilename('fullpath'));

% 自分自身の場所（stage1/）に移動する
% これにより、Dynare のビルドファイルが stage1 内に生成され、ルートを汚しません
cd(mpath);

% --- [2] パラメータの橋渡し ---
% .mod ファイル内の Steady State 計算や式で使用するパラメータをベースワークスペースへ
% Stage 1 では標準的なインフレ反応係数（通常1.5など）を渡します
assignin('base', 'phi_pi_H_value', phi_pi_H_value);
assignin('base', 'phi_pi_F_value', phi_pi_F_value);

% --- [3] モデルの初期化と実行 ---
% settings.modnam = 'tr_cpi_stage1_model' を参照
% noclearall: 既存の変数を消さずに実行 / nolog: コマンドウィンドウへのログを抑制
eval(['dynare ', settings.modnam, ' noclearall nolog']);

% シミュレーション設定の展開
nperiods = settings.nperiods;
irfshock = settings.irfshock;
shock_magnitude = settings.shockssequence;

% --- [4] ショック行列の作成 ---
% 全期間ゼロで初期化されたショック行列を作成
ex_flex = zeros(nperiods, M_.exo_nbr);

% 指定されたショック名のインデックスを特定
shock_idx = strmatch(irfshock, M_.exo_names, 'exact');
if isempty(shock_idx)
    error('指定されたショック変数 %s がモデル内に見つかりません。', irfshock);
end

% 第1期（t=1）にショック値（大きさ）を代入
% shockssequenceは列ベクトルのため、1行目のみ抽出
if length(shock_magnitude) > 1
    ex_flex(1:nperiods, shock_idx) = shock_magnitude;
else
    ex_flex(1, shock_idx) = shock_magnitude;
end

% --- [5] シミュレーション実行 ---
% Stage 1 は伸縮価格モデルのため、1次近似（Linear）の simult_ を使用します
fprintf('Stage 1 (Flexible Price) シミュレーションを実行中... シナリオ: %s\n', irfshock);

y0 = oo_.dr.ys; % 定常状態の値を取得
dr = oo_.dr;   % 決定ルール（Decision Rules）
iorder = 1;    % 1次近似（Linear approximation）

% シミュレーション実行（1期からnperiods期まで）
y_flex_sim = simult_(M_, options_, y0, dr, ex_flex, iorder);

% --- [6] 呼び出し元のディレクトリ（src）に戻る ---
cd(root_dir);

% --- [7] 結果の展開 (偏差として取得) ---
% Dynareの生出力 y_flex_sim から定常状態 y0 を引き、偏差データを生成
for i=1:M_.endo_nbr
    % _pie 命名規則 (偏差データ)
    eval([M_.endo_names{i,:},'_pie = y_flex_sim(i, 2:end)'' - y0(i);']);
    % _lin も同一値として作成（Stage 1 では比較対象がないため）
    eval([M_.endo_names{i,:},'_lin = y_flex_sim(i, 2:end)'' - y0(i);']);
end

% グラフ描画および Stage 2 への受け渡し用としてマッピング
% Stage 1 の y_H (伸縮価格下) こそが Stage 2 における Potential Output となる
y_H_potential_pie = y_H_pie;
y_H_potential_lin = y_H_lin;

%% --- [8] パラメータと定常状態の一括取得 ---
% 厚生計算や補助変数の生成に必要な定数値をモデルから抽出
N_val         = M_.params(strcmp('N', M_.param_names));
M_val         = M_.params(strcmp('M', M_.param_names));
phi_H_val     = M_.params(strcmp('phi_H', M_.param_names));
phi_F_val     = M_.params(strcmp('phi_F', M_.param_names));
beta_H_ss_val = M_.params(strcmp('beta_H_ss', M_.param_names));
beta_F_ss_val = M_.params(strcmp('beta_F_ss', M_.param_names));

% 各変数の定常状態水準値
e_slash_star_ss_val = y0(strcmp('e_slash_star', M_.endo_names));
c_H_W_ss_val         = y0(strcmp('c_H_W', M_.endo_names));
l_H_ss_val           = y0(strcmp('l_H', M_.endo_names));
c_F_W_ss_val         = y0(strcmp('c_F_W', M_.endo_names));
l_F_ss_val           = y0(strcmp('l_F', M_.endo_names));
b_H_ss_val           = y0(strcmp('b_H', M_.endo_names)); 

% 金融政策の参照用に、正規化されたCPIの定常状態値も取得 (命名規則の修正)
p_H_W_bar_ss_val      = y0(strcmp('p_H_W_bar', M_.endo_names));
p_F_W_star_bar_ss_val = y0(strcmp('p_F_W_star_bar', M_.endo_names)); % 【修正箇所】命名規則を一致

%% --- [9] 追加変数の計算 ---
% 水準（Level）データの復元（偏差＋定常状態）
e_slash_star_level_pie = e_slash_star_ss_val + e_slash_star_pie;
e_slash_star_level_lin = e_slash_star_ss_val + e_slash_star_lin;

b_H_level_pie = b_H_ss_val + b_H_pie;
b_H_level_lin = b_H_ss_val + b_H_lin;

% 為替レートと Home 債券残高から Foreign 債券残高を計算
b_F_star_level_pie = - (N_val/M_val) * (b_H_level_pie ./ e_slash_star_level_pie);
b_F_star_level_lin = - (N_val/M_val) * (b_H_level_lin ./ e_slash_star_level_lin);

b_F_star_ss = - (N_val/M_val) * (b_H_ss_val / e_slash_star_ss_val);

% 再度、偏差データとしての b_F_star を作成
b_F_star_pie = b_F_star_level_pie - b_F_star_ss;
b_F_star_lin = b_F_star_level_lin - b_F_star_ss;

%% --- [10] 効用 (Utility) と厚生 (Welfare) の計算 ---

% 1. 各国パラメータの準備
theta_H_val = M_.params(strcmp('theta_H', M_.param_names));
xi_H_val    = M_.params(strcmp('xi_H', M_.param_names));
theta_F_val = M_.params(strcmp('theta_F', M_.param_names));
xi_F_val    = M_.params(strcmp('xi_F', M_.param_names));

% 2. 8パターン（Home/Foreign × Pie/Lin × With/Without Delta）の初期化
Delta_H_pie = zeros(nperiods, 1); Delta_F_pie = zeros(nperiods, 1);
Delta_H_lin = zeros(nperiods, 1); Delta_F_lin = zeros(nperiods, 1);

utility_H_with_delta_pie = zeros(nperiods, 1); utility_F_with_delta_pie = zeros(nperiods, 1);
utility_H_without_delta_pie = zeros(nperiods, 1); utility_F_without_delta_pie = zeros(nperiods, 1);
utility_H_with_delta_lin = zeros(nperiods, 1); utility_F_with_delta_lin = zeros(nperiods, 1);
utility_H_without_delta_lin = zeros(nperiods, 1); utility_F_without_delta_lin = zeros(nperiods, 1);

% 伸縮価格モデル（Stage 1）ではインフレ変数がゼロになる場合があるため初期化
if ~exist('pi_H_pie','var'), pi_H_pie = zeros(nperiods,1); end
if ~exist('pi_F_star_pie','var'), pi_F_star_pie = zeros(nperiods,1); end
if ~exist('pi_H_lin','var'), pi_H_lin = zeros(nperiods,1); end
if ~exist('pi_F_star_lin','var'), pi_F_star_lin = zeros(nperiods,1); end

% 3. 時系列ループ計算
for t = 1:nperiods
    if t == 1
        DH_lag_p = 0; DF_lag_p = 0; DH_lag_l = 0; DF_lag_l = 0;
    else
        DH_lag_p = Delta_H_pie(t-1); DF_lag_p = Delta_F_pie(t-1);
        DH_lag_l = Delta_H_lin(t-1); DF_lag_l = Delta_F_lin(t-1);
    end
    
    % 価格分散項 Delta (Stage 1 では xi=0 のため常に 0 となるが、定義式を記述)
    Delta_H_pie(t) = xi_H_val * DH_lag_p + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_pie(t)^2);
    Delta_F_pie(t) = xi_F_val * DF_lag_p + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_pie(t)^2);
    Delta_H_lin(t) = xi_H_val * DH_lag_l + 0.5 * (theta_H_val * xi_H_val / (1 - xi_H_val)) * (pi_H_lin(t)^2);
    Delta_F_lin(t) = xi_F_val * DF_lag_l + 0.5 * (theta_F_val * xi_F_val / (1 - xi_F_val)) * (pi_F_star_lin(t)^2);

    % --- パターンA: Deltaあり (with_delta) ---
    utility_H_with_delta_pie(t) = c_H_W_pie(t) - phi_H_val*(l_H_ss_val^2)*(y_H_pie(t)-a_H_pie(t)) - phi_H_val*(l_H_ss_val^2)*Delta_H_pie(t);
    utility_F_with_delta_pie(t) = c_F_W_pie(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t)) - phi_F_val*(l_F_ss_val^2)*Delta_F_pie(t);
    utility_H_with_delta_lin(t) = c_H_W_lin(t) - phi_H_val*(l_H_ss_val^2)*(y_H_lin(t)-a_H_lin(t)) - phi_H_val*(l_H_ss_val^2)*Delta_H_lin(t);
    utility_F_with_delta_lin(t) = c_F_W_lin(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t)) - phi_F_val*(l_F_ss_val^2)*Delta_F_lin(t);

    % --- パターンB: Deltaなし (without_delta) ---
    utility_H_without_delta_pie(t) = c_H_W_pie(t) - phi_H_val*(l_H_ss_val^2)*(y_H_pie(t)-a_H_pie(t));
    utility_F_without_delta_pie(t) = c_F_W_pie(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t));
    utility_H_without_delta_lin(t) = c_H_W_lin(t) - phi_H_val*(l_H_ss_val^2)*(y_H_lin(t)-a_H_lin(t));
    utility_F_without_delta_lin(t) = c_F_W_lin(t) - phi_F_val*(l_F_ss_val^2)*(y_F_pie(t)-a_F_pie(t));
end

% 4. 割引現在価値の合計（生涯厚生 Welfare）
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

%% --- [11] その他の連結変数などの生成 ---
% 為替レートと他変数との積・商の関係にある変数を対数近似（和・差）で計算
lambda_H_slash_star_pie = e_slash_star_pie + lambda_H_pie; 
lambda_H_slash_star_lin = e_slash_star_lin + lambda_H_lin;
lambda_F_pie = lambda_H_slash_star_pie - e_slash_star_pie; 
lambda_F_lin = lambda_F_slash_star_lin - e_slash_star_lin;

% 名目GDPおよび名目総消費
p_H_bar_y_H_pie = p_H_bar_pie + y_H_pie; 
p_H_bar_y_H_lin = p_H_bar_lin + y_H_lin;
% 名目総消費の計算（理論定義 pHW * cHW を維持）
p_H_W_c_H_W_pie = p_H_W_pie + c_H_W_pie;
p_H_W_c_H_W_lin = p_H_W_lin + c_H_W_lin; 

% --- [12] 戻り値の整理 (Packaging) ---
% who を使用してワークスペース上の関連変数を一括で構造体にパッキング
vars = [who('*_pie'); who('*_lin'); who('Delta_H_*'); who('Delta_F_*'); who('utility_*')];
for k = 1:length(vars)
    irf_results.(vars{k}) = eval(vars{k});
end

% Dynare内部構造体の保持
Mbase_ = M_;
oobase_ = oo_;

end