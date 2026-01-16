// src/it_ppi_model_zlb.mod

@#include "it_ppi_model_declarations.inc"

model;
    @#include "it_ppi_model_equations_common.inc"
    
    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;

    // 自国は生産者物価インフレ目標 (IT PPI) で計算
    i_H_notional = i_H_ss + phi_pi_H*(pi_H - pi_H_ss) + eps_i_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;