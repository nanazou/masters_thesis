// src/stage2/tr_ppi_stage2_model_zlb.mod

// Stage 2 専用の宣言ファイルを読み込む (修正: tr_ppi_stage2_model_declarations.inc)
@#include "tr_ppi_stage2_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む (修正: パスを ../tr_ppi_model_equations_common.inc に)
    @#include "../tr_ppi_model_equations_common.inc"

    // --- IV. 潜在産出量の定義 ---
    // Stage 1（伸縮価格モデル）で算出した産出量パスを外生ショックとして受け取る
    y_H_potential = y_H_ss + eps_y_H_potential;

    // --- V. 金融政策ルール ---
    // 外国は生産者物価物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は生産者物価テイラー・ルール（TR PPI）
    i_H_notional = rho_i_H*i_H_notional(-1) + (1-rho_i_H)*(i_H_ss + phi_pi_H*log(pi_H) + phi_y_H*(log(y_H)-log(y_H_potential))) + eps_i_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;