// src/olt_model_zlb.mod

@#include "olt_model_declarations.inc"

model;
    @#include "olt_model_equations_common.inc"
    
    // --- IV. ショックと追加のプロセス ---
    // 目標パス (AR(1)プロセス)
    log(chi_H) = (1-rho_chi_H)*log(chi_H_ss) + rho_chi_H*log(chi_H(-1)) + eps_chi_H;

    // 過去の産出量ギャップの累積 (履歴効果/コミットメント)
    // 産出量 = y_H
    gamma_H = gamma_H(-1) + (log(y_H(-1)) - log(chi_H(-1)));

    // --- V. 金融政策ルール ---
    // 外国は生産者物価物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は生産水準目標 (OLT)
    i_H_notional = i_H_ss + phi_gap_H*(log(y_H)-log(chi_H)) + phi_level_H*gamma_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;