// src/nclt_model_zlb.mod

@#include "nclt_model_declarations.inc"

model;
    @#include "nclt_model_equations_common.inc"

    // --- IV. ショックと追加のプロセス ---
    // 目標パス (AR(1)プロセス)
    log(chi_H)  = (1-rho_chi_H)*log(chi_H_ss)   + rho_chi_H*log(chi_H(-1))   + eps_chi_H; 
    
    // 過去の「名目消費」ギャップの累積
    // 名目消費 = p_H_W * c_H_W
    gamma_H = gamma_H(-1) + (log(p_H_W(-1)*c_H_W(-1)) - log(chi_H(-1)));
    
    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;     
    
    // 自国は名目総消費水準目標（NCLT）
    i_H_notional = i_H_ss + phi_gap_H * (log(p_H_W*c_H_W) - log(chi_H)) + phi_level_H * gamma_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;