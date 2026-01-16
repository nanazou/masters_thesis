// src/it_cpi_model_zlb.mod

@#include "it_cpi_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む
    @#include "it_cpi_model_equations_common.inc"
    
    // --- V. 金融政策ルール ---
    // 外国は生産者物価物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;

    // 自国は消費者物価インフレ目標 (IT CPI) 
    i_H_notional = i_H_ss + phi_pi_H*(pi_H_W - pi_H_W_ss) + eps_i_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;