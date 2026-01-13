// src/clt_model_zlb.mod

@#include "clt_model_declarations.inc"

model;
    @#include "clt_model_equations_common.inc"
    
    // --- IV. ショックと追加のプロセス ---
    // 目標パス
    log(chi_H) = (1-rho_chi_H)*log(chi_H_ss) + rho_chi_H*log(chi_H(-1)) + eps_chi_H;

    // 過去の実質消費ギャップの累積 (コミットメント)
    gamma_H = gamma_H(-1) + (log(c_H_W(-1)) - log(chi_H(-1)));

    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は総消費水準目標 (CLT)
    i_H_notional = i_H_ss + phi_gap_H*(log(c_H_W)-log(chi_H)) + phi_level_H*gamma_H; 
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;