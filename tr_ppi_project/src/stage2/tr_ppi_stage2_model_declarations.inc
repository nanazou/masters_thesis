// src/stage2/tr_ppi_stage2_model_declarations.inc

// --- 1. 内生変数の宣言 ---
var y_H l_H c_H_W c_H_H c_H_F t_H p_H_bar p_F_bar p_H_W p_H_W_bar pi_H pi_H_W lambda_H b_H v_H w_H p_H_tilde beta_H a_H tau_H p_H Delta_H
    y_F l_F c_F_W c_F_F c_F_H t_F p_F_star_bar p_H_star_bar p_F_W_star p_F_W_star_bar pi_F_star pi_F_W_star lambda_F_slash_star v_F w_F p_F_star_tilde beta_F a_F tau_F p_F_star Delta_F
    e_slash_star
    i_H i_F i_H_notional
    y_H_potential; // Stage 2用: Stage 1（伸縮価格モデル）から引き継ぐ潜在産出量のパス

// --- 2. 外生ショック変数の宣言 ---
varexo eps_beta_H eps_a_H eps_tau_H eps_i_H
       eps_beta_F eps_a_F eps_tau_F eps_i_F
       eps_y_H_potential; // Stage 2用: Stage 1で計算された「伸縮価格下のy_H」の変動分を外部投入する

// --- 3. パラメータの宣言 ---
parameters theta_H theta_F phi_H phi_F alpha_H alpha_F N M
           xi_H xi_F
           beta_H_ss rho_beta_H
           beta_F_ss rho_beta_F
           a_H_ss rho_a_H
           a_F_ss rho_a_F
           tau_H_ss rho_tau_H
           tau_F_ss rho_tau_F
           phi_pi_H phi_y_H rho_i_H // 自国(H)のテイラー・ルール反応係数
           phi_pi_F                   // 外国(F)のインフレ目標反応係数
           y_H_ss pi_H_ss i_H_ss pi_H_W_ss
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

// 自国(H)のテイラー・ルール (phi_pi_H_value, phi_y_H_value を最適化ループ等から受け取る)
rho_i_H    = 0.8;
phi_pi_H   = phi_pi_H_value;
phi_y_H    = phi_y_H_value;

// 外国(F)のテイラー・ルール (インフレ目標)
phi_pi_F   = 1.5;

// 定常状態の値
pi_H_ss = 1.0;
pi_F_star_ss = 1.0;
pi_H_W_ss = 1.0;
pi_F_W_star_ss = 1.0;
i_H_ss = 1/beta_H_ss * pi_H_ss - 1;
i_F_ss = 1/beta_F_ss * pi_F_star_ss - 1;
y_H_ss = 1.0;
y_F_ss = 1.0;