// src/stage2/tr_cpi_stage2_model_zlb.mod

@#include "tr_cpi_stage2_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む (修正: パスを ../tr_cpi_model_equations_common.inc に)
    @#include "../tr_cpi_model_equations_common.inc"

    // --- IV. 潜在産出量の定義 ---
    // Stage 1（伸縮価格モデル）で算出した産出量パスを外生ショックとして受け取る
    y_H_potential = y_H_ss + eps_y_H_potential;

    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は消費者物価テイラールール（TR CPI）
    i_H_notional = rho_i_H*i_H_notional(-1) + (1-rho_i_H)*(i_H_ss + phi_pi_H*log(pi_H_W) + phi_y_H*(log(y_H)-log(y_H_potential))) + eps_i_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;