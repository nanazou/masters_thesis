// src/cplt_model_zlb.mod

@#include "cplt_model_declarations.inc"

model;
    @#include "cplt_model_equations_common.inc"
    
    // --- IV. ショックと追加의 プロセス ---
    // 目標パス
    log(chi_H) = (1-rho_chi_H)*log(chi_H_ss) + rho_chi_H*log(chi_H(-1)) + eps_chi_H;

    // 過去のCPIギャップの累積 (コミットメント)
    gamma_H = gamma_H(-1) + (log(p_H_W(-1)) - log(chi_H(-1)));

    // --- V. 金融政策ルール ---
    // 外国は生産者物価物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は消費者物価水準目標 (CPLT)で計算
    i_H_notional = i_H_ss + phi_gap_H*(log(p_H_W)-log(chi_H)) + phi_level_H*gamma_H; 
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;