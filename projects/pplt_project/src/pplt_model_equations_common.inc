// src/pplt_model_equations_common.inc

// --- I. 家計の最適化条件 (5本) ---
lambda_H = 1 / (p_H_W * c_H_W);
lambda_F_slash_star = 1 / (p_F_W_star * c_F_W);
lambda_H = beta_H * (1+i_H) * lambda_H(+1);
lambda_F_slash_star = beta_F * (1+i_F) * lambda_F_slash_star(+1);
lambda_H * e_slash_star = beta_H * (1+i_F) * lambda_H(+1) * e_slash_star(+1);

// --- II. 価格設定（非線形モデル） (8本) ---
v_H = phi_H * N^2 * (y_H/a_H)^2 * p_H^(2*theta_H) + beta_H * xi_H * v_H(+1);
w_H = lambda_H * N * (1-tau_H) * y_H * p_H^(theta_H) + beta_H * xi_H * w_H(+1);
p_H_tilde^(1+theta_H) = (theta_H/(theta_H-1)) * (v_H / w_H);
p_H_bar^(1-theta_H) = (1-xi_H)*p_H_tilde^(1-theta_H) + xi_H*p_H_bar(-1)^(1-theta_H);

v_F = phi_F * M^2 * (y_F/a_F)^2 * p_F_star^(2*theta_F) + beta_F * xi_F * v_F(+1);
w_F = lambda_F_slash_star * M * (1-tau_F) * y_F * p_F_star^(theta_F) + beta_F * xi_F * w_F(+1);
p_F_star_tilde^(1+theta_F) = (theta_F/(theta_F-1)) * (v_F / w_F);
p_F_star_bar^(1-theta_F) = (1-xi_F)*p_F_star_tilde^(1-theta_F) + xi_F*p_F_star_bar(-1)^(1-theta_F);

// --- III. 生産、需要、各種定義式 (25本) ---
p_H = (N^(1/(1-theta_H)))*p_H_bar;
p_F_star = (M^(1/(1-theta_F)))*p_F_star_bar;

// 生産関数に Delta を追加
y_H = a_H*l_H / Delta_H;
y_F = a_F*l_F / Delta_F;

// 価格分散の動学方程式
Delta_H = (1-xi_H) * (p_H_tilde / p_H_bar)^(-theta_H) + xi_H * (pi_H)^theta_H * Delta_H(-1);
Delta_F = (1-xi_F) * (p_F_star_tilde / p_F_star_bar)^(-theta_F) + xi_F * (pi_F_star)^theta_F * Delta_F(-1);

// 理論上の消費者物価指数（生計費指数：予算制約や効用計算に使用）
p_H_W = p_H^alpha_H * (e_slash_star*p_F_star)^(1-alpha_H);
p_F_W_star = p_F_star^alpha_F * (p_H/e_slash_star)^(1-alpha_F);

// 正規化された消費者物価指数（統計局ベース：金融政策のターゲットに使用）
p_H_W_bar = p_H_bar^alpha_H * (e_slash_star*p_F_star_bar)^(1-alpha_H);
p_F_W_star_bar = p_F_star_bar^alpha_F * (p_H_bar/e_slash_star)^(1-alpha_F);

c_H_H = alpha_H * (p_H_W / p_H) * c_H_W;
c_H_F = (1-alpha_H) * (p_H_W / (e_slash_star*p_F_star)) * c_H_W;
c_F_F = alpha_F * (p_F_W_star / p_F_star) * c_F_W;
c_F_H = (1-alpha_F) * (p_F_W_star / (p_H/e_slash_star)) * c_F_W;

t_H = tau_H*p_H*y_H;
t_F = tau_F*p_F_star*y_F;
p_F_bar = e_slash_star*p_F_star_bar;
p_H_star_bar = p_H_bar/e_slash_star;

// インフレ率の定義（正規化された指数に基づく）
pi_H = p_H_bar/p_H_bar(-1);
pi_F_star = p_F_star_bar/p_F_star_bar(-1);
pi_H_W = p_H_W_bar / p_H_W_bar(-1);
pi_F_W_star = p_F_W_star_bar / p_F_W_star_bar(-1);

y_H = c_H_H + (M/N)*c_F_H;
y_F = c_F_F + (N/M)*c_H_F;

p_H_W*c_H_W + b_H = p_H*y_H + (1+i_F(-1))*(e_slash_star/e_slash_star(-1))*b_H(-1);

// --- IV. ショックのプロセス (6本) ---
log(beta_H) = (1-rho_beta_H)*log(beta_H_ss) + rho_beta_H*log(beta_H(-1)) + eps_beta_H;
log(a_H)    = (1-rho_a_H)*log(a_H_ss)       + rho_a_H*log(a_H(-1))       + eps_a_H;
log(tau_H)  = (1-rho_tau_H)*log(tau_H_ss)   + rho_tau_H*log(tau_H(-1))   + eps_tau_H;
log(beta_F) = (1-rho_beta_F)*log(beta_F_ss) + rho_beta_F*log(beta_F(-1)) + eps_beta_F;
log(a_F)    = (1-rho_a_F)*log(a_F_ss)       + rho_a_F*log(a_F(-1))       + eps_a_F;
log(tau_F)  = (1-rho_tau_F)*log(tau_F_ss)   + rho_tau_F*log(tau_F(-1))   + eps_tau_F;