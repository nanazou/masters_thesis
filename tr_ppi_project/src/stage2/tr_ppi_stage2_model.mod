// src/stage2/tr_ppi_stage2_model.mod

// Stage 2 専用の宣言ファイルを読み込む (修正を反映)
@#include "tr_ppi_stage2_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む (修正: パスを ../tr_ppi_model_equations_common.inc に)
    @#include "../tr_ppi_model_equations_common.inc"

    // --- IV. 潜在産出量の定義 ---
    // Stage 1 で算出した潜在産出量パスを、外部ショック eps_y_H_potential として反映
    y_H_potential = y_H_ss + eps_y_H_potential;

    // --- V. 金融政策ルール ---
    // 外国は生産者物価物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国はテイラー・ルール (TR)
    i_H_notional = rho_i_H*i_H_notional(-1) + (1-rho_i_H)*(i_H_ss + phi_pi_H*log(pi_H) + phi_y_H*(log(y_H)-log(y_H_potential))) + eps_i_H;
    
    i_H = i_H_notional;
end;

steady_state_model;
    beta_H = beta_H_ss;
    beta_F = beta_F_ss;
    a_H = a_H_ss;
    a_F = a_F_ss;
    tau_H = tau_H_ss;
    tau_F = tau_F_ss;
    
    // Stage 2 でも定常状態では履歴効果やギャップはゼロ
    gamma_H = 0;
    
    y_H = y_H_ss;
    y_F = y_F_ss;
    y_H_potential = y_H_ss; // 定常状態では実力通り
    pi_H = pi_H_ss;
    pi_F_star = pi_F_star_ss;

    // 消費者物価インフレ率の定常状態 (pi_F_W_star に統一)
    pi_H_W = pi_H_W_ss;
    pi_F_W_star = pi_F_W_star_ss;

    i_H = i_H_ss;
    i_F = i_F_ss;
    i_H_notional = i_H_ss;
    b_H = 0;
    
    e_slash_star = 1;

    l_H = y_H / a_H;
    l_F = y_F / a_F;
    
    p_H_bar = 1;
    p_F_star_bar = 1;
    
    p_F_bar = e_slash_star * p_F_star_bar;
    p_H_star_bar = p_H_bar / e_slash_star;
    
    // ★ 追加した変数の定常状態の計算 (p_H_W_bar, p_F_W_star_bar)
    p_H_W_bar = p_H_bar^alpha_H * (e_slash_star * p_F_star_bar)^(1-alpha_H);
    p_F_W_star_bar = p_F_star_bar^alpha_F * (p_H_bar / e_slash_star)^(1-alpha_F);

    p_H = (N^(1/(1-theta_H)))*p_H_bar;
    p_F_star = (M^(1/(1-theta_F)))*p_F_star_bar;
    
    p_H_W = p_H^alpha_H * (e_slash_star*p_F_star)^(1-alpha_H);
    p_F_W_star = p_F_star^alpha_F * (p_H/e_slash_star)^(1-alpha_F);
    
    c_H_W = (p_H/p_H_W)*y_H;
    c_F_W = (p_F_star/p_F_W_star)*y_F;

    c_H_H = alpha_H * (p_H_W / p_H) * c_H_W;
    c_H_F = (1-alpha_H) * (p_H_W / (e_slash_star*p_F_star)) * c_H_W;
    c_F_F = alpha_F * (p_F_W_star / p_F_star) * c_F_W;
    c_F_H = (1-alpha_F) * (p_F_W_star / (p_H/e_slash_star)) * c_F_W;
    
    lambda_H = 1 / (p_H_W * c_H_W);
    lambda_F_slash_star = 1 / (p_F_W_star * c_F_W);
    t_H = tau_H * p_H * y_H;
    t_F = tau_F * p_F_star * y_F;
    
    p_H_tilde = p_H_bar;
    p_F_star_tilde = p_F_star_bar;
    
    v_H = (phi_H * N^2 * (y_H/a_H)^2 * p_H^(2*theta_H)) / (1 - beta_H*xi_H);
    w_H = (lambda_H * N * (1-tau_H) * y_H * p_H^(theta_H)) / (1 - beta_H*xi_H);
    
    v_F = (phi_F * M^2 * (y_F/a_F)^2 * p_F_star^(2*theta_F)) / (1 - beta_F*xi_F);
    w_F = (lambda_F_slash_star * M * (1-tau_F) * y_F * p_F_star^(theta_F)) / (1 - beta_F*xi_F);

    Delta_H = 1;
    Delta_F = 1;
end;

// 粘着価格設定 (POLT Stage 2 と同様に設定)
xi_H = 0.99;
xi_F = 0.75;

check;
steady;
// デフォルトはorder=2だが、occbinはorder=1しか使用しない
stoch_simul(order=1, irf=0);