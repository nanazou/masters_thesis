// src/stage2/polt_stage2_model_zlb.mod

// Stage 2 専用の宣言ファイルを読み込む
@#include "polt_stage2_model_declarations.inc"

model;
    // 共通数式ファイルを読み込む
    @#include "../polt_model_equations_common.inc"
    
    // --- IV. ショックと追加のプロセス ---
    // 潜在産出量 (Stage 1 で計算された自然産出量のデータをショックとして受け取る)
    y_H_potential = y_H_ss + eps_y_H_potential;

    // 過去の産出量ギャップの累積 (履歴効果 / Make-up コミットメント)
    // 実際の産出量 (y_H) と 目標となる自然産出量 (y_H_potential) の乖離を蓄積
    gamma_H = gamma_H(-1) + (log(y_H(-1)) - log(y_H_potential(-1)));

    // --- V. 金融政策ルール ---
    // 外国は生産者物価インフレ目標 (IT PPI)
    i_F = i_F_ss + phi_pi_F*(pi_F_star - pi_F_star_ss) + eps_i_F;
    
    // 自国は潜在生産水準目標 (POLT) 
    i_H_notional = i_H_ss + phi_gap_H*(log(y_H)-log(y_H_potential)) + phi_level_H*gamma_H;
    
    // ZLB時は、実際の金利(i_H)をゼロに固定
    i_H = 0;
end;