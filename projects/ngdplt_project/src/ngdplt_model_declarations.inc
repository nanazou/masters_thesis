// src/ngdplt_model_declarations.inc

// --- 1. 内生変数の宣言 ---
var y_H l_H c_H_W c_H_H c_H_F t_H p_H_bar p_F_bar p_H_W p_H_W_bar pi_H pi_H_W lambda_H b_H v_H w_H p_H_tilde beta_H a_H tau_H p_H gamma_H chi_H Delta_H
    y_F l_F c_F_W c_F_F c_F_H t_F p_F_star_bar p_H_star_bar p_F_W_star p_F_W_star_bar pi_F_star pi_F_W_star lambda_F_slash_star v_F w_F p_F_star_tilde beta_F a_F tau_F p_F_star Delta_F
    e_slash_star
    i_H i_F i_H_notional;

// --- 2. 外生ショック変数の宣言 ---
varexo eps_beta_H eps_a_H eps_tau_H eps_chi_H
       eps_beta_F eps_a_F eps_tau_F eps_i_F;

// --- 3. パラメータの宣言 ---
parameters theta_H theta_F phi_H phi_F alpha_H alpha_F N M
           xi_H xi_F
           beta_H_ss rho_beta_H
           beta_F_ss rho_beta_F
           a_H_ss rho_a_H
           a_F_ss rho_a_F
           tau_H_ss rho_tau_H
           tau_F_ss rho_tau_F
           phi_gap_H phi_level_H rho_chi_H chi_H_ss // NGDPLT用のパラメータ
           phi_pi_F                                 // 外国(テイラー・ルール)用のパラメータ
           y_H_ss pi_H_ss i_H_ss pi_H_W_ss          // CPI定常状態を追加
           y_F_ss pi_F_star_ss i_F_ss pi_F_W_star_ss;

// --- 4. パラメータのキャリブレーション ---
beta_H_ss  = 0.995;
rho_beta_H = 0.977;
beta_F_ss  = 0.995;
rho_beta_F = 0.9;
theta_H    = 11.0;
theta_F    = 11.0;
tau_H_ss    = 0.2;
rho_tau_H  = 0.9;
tau_F_ss    = 0.2;
rho_tau_F  = 0.9;
phi_H      = 8/11;
phi_F      = 8/11;
a_H_ss     = 1.0;
rho_a_H    = 0.95;
a_F_ss     = 1.0;
rho_a_F    = 0.95;
alpha_H    = 0.85;
alpha_F    = 0.85;
xi_H       = 0.99;
xi_F       = 0.75;
N          = 1.0;
M          = 1.0;

// 外国の政策パラメータ
phi_pi_F   = 1.5;

// 自国(H)の名目GDP水準目標用パラメータ
// MATLABのワークスペースから値を読み込む
phi_gap_H    = phi_gap_H_value;
phi_level_H  = phi_level_H_value;
rho_chi_H  = 0.9;
chi_H_ss   = 1.0;

// 定常状態の値
pi_H_ss = 1.0;
pi_F_star_ss = 1.0;
pi_H_W_ss = 1.0;           // 自国CPIインフレ定常状態
pi_F_W_star_ss = 1.0;      // 外国CPIインフレ定常状態
i_H_ss = 1/beta_H_ss * pi_H_ss - 1;
i_F_ss = 1/beta_F_ss * pi_F_star_ss - 1;
y_H_ss = 1.0;
y_F_ss = 1.0;