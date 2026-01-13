// src/stage2/polt_stage2_model.mod

// Stage 2 専用の宣言ファイルを読み込む
@#include "polt_stage2_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む
    @#include "../polt_model_equations_common.inc"

    // --- IV. 潜在生産とギャップの累積 ---
    // Stage 1 で計算した生産パスを外生ショック eps_y_H_potential として受け取る
    y_H_potential = y_H_ss + eps_y_H_potential;

    // 過去の生産ギャップの累積 (履歴効果 / コミットメント項)
    gamma_H = gamma_H(-1) + (log(y_H(-1)) - log(y_H_potential(-1)));

    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;

    // 自国は潜在生産水準目標 (POLT)
    i_H_notional = i_H_ss + phi_gap_H*(log(y_H) - log(y_H_potential)) + phi_level_H*gamma_H;
    i_H = i_H_notional;
end;

steady_state_model;
    beta_H = beta_H_ss; beta_F = beta_F_ss;
    a_H = a_H_ss; a_F = a_F_ss;
    tau_H = tau_H_ss; tau_F = tau_F_ss;
    gamma_H = 0;
    y_H = y_H_ss; y_F = y_F_ss;
    y_H_potential = y_H_ss;
    pi_H = pi_H_ss; pi_F_star = pi_F_star_ss;
    
    // 消費者物価インフレ率の定常状態
    pi_H_W = pi_H_W_ss;
    pi_F_W_star = pi_F_W_star_ss;

    i_H = i_H_ss; i_F = i_F_ss; i_H_notional = i_H_ss;
    b_H = 0;
    
    e_slash_star = 1;
    
    l_H = y_H / a_H; l_F = y_F / a_F;

    // 命名規則に従い修正: star_bar の順序に統一
    p_H_bar = 1; p_F_star_bar = 1; 
    
    p_F_bar = e_slash_star * p_F_star_bar; 
    p_H_star_bar = p_H_bar / e_slash_star;
    
    p_H = (N^(1/(1-theta_H)))*p_H_bar; 
    p_F_star = (M^(1/(1-theta_F)))*p_F_star_bar;
    
    // 理論的な消費者物価指数（生計費指数）
    p_H_W = p_H^alpha_H * (e_slash_star*p_F_star)^(1-alpha_H);
    p_F_W_star = p_F_star^alpha_F * (p_H/e_slash_star)^(1-alpha_F);

    // 【重要】統計ベースの消費者物価指数 (値札ベース) の定常状態を追加
    p_H_W_bar = p_H_bar^alpha_H * (e_slash_star*p_F_star_bar)^(1-alpha_H);
    p_F_W_star_bar = p_F_star_bar^alpha_F * (p_H_bar/e_slash_star)^(1-alpha_F);
    
    c_H_W = (p_H/p_H_W)*y_H; c_F_W = (p_F_star/p_F_W_star)*y_F;
    c_H_H = alpha_H * (p_H_W / p_H) * c_H_W; 
    
    c_H_F = (1-alpha_H) * (p_H_W / (e_slash_star*p_F_star)) * c_H_W;
    c_F_F = alpha_F * (p_F_W_star / p_F_star) * c_F_W; 
    c_F_H = (1-alpha_F) * (p_F_W_star / (p_H/e_slash_star)) * c_F_W;
    
    lambda_H = 1 / (p_H_W * c_H_W); lambda_F_slash_star = 1 / (p_F_W_star * c_F_W);
    t_H = tau_H * p_H * y_H; t_F = tau_F * p_F_star * y_F;

    // 命名規則に従い修正: star_tilde の順序に統一
    p_H_tilde = p_H_bar; 
    p_F_star_tilde = p_F_star_bar;
    
    v_H = (phi_H * N^2 * (y_H/a_H)^2 * p_H^(2*theta_H)) / (1 - beta_H*xi_H);
    w_H = (lambda_H * N * (1-tau_H) * y_H * p_H^(theta_H)) / (1 - beta_H*xi_H);
    v_F = (phi_F * M^2 * (y_F/a_F)^2 * p_F_star^(2*theta_F)) / (1 - beta_F*xi_F);
    w_F = (lambda_F_slash_star * M * (1-tau_F) * y_F * p_F_star^(theta_F)) / (1 - beta_F*xi_F);
    
    Delta_H = 1; Delta_F = 1;
end;

// ★ 粘着価格パラメータの設定 (Stage 1 と異なり価格硬直性を有効にする)
xi_H = 0.99;
xi_F = 0.75;

check;
steady;
// デフォルトはorder=2だが、occbinはorder=1しか使用しない
stoch_simul(order=1, irf=0);