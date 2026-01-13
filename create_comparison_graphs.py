import scipy.io
import matplotlib.pyplot as plt
import japanize_matplotlib # 日本語表示のためにインポート
import numpy as np
import os
import sys

# =========================================================================
# --- 1. 設定項目 ---
# =========================================================================

# 比較する政策のリスト定義
policies = [
    {'key': 'clt',      'label': '総消費水準目標（CLT）',    'color': 'orange',      'ls': ':'},
    {'key': 'cplt',     'label': '消費者物価水準目標（CPLT）', 'color': 'magenta',     'ls': '-'},
    {'key': 'it_cpi',   'label': '消費者物価インフレ目標（IT-CPI）', 'color': 'darkred',     'ls': '-.'},
    {'key': 'it_ppi',   'label': '生産者物価インフレ目標（IT-PPI）', 'color': 'brown',       'ls': ':'},
    {'key': 'nclt',     'label': '名目総消費水準目標（NCLT）',   'color': 'blue',          'ls': '-'},
    {'key': 'ngdplt',   'label': '名目GDP水準目標（NGDPLT）',  'color': 'green',       'ls': '-.'},
    {'key': 'olt',      'label': '生産水準目標（OLT）',      'color': '#FFD700',     'ls': '-.'},
    {'key': 'polt',     'label': '潜在生産水準目標（POLT）',   'color': 'black',        'ls': '--'},
    {'key': 'pplt',      'label': '生産者物価水準目標（PPLT）',    'color': 'cyan',         'ls': ':'},
    {'key': 'tr_cpi',   'label': '消費者物価テイラールール（TR-CPI）',      'color': '#DC143C',      'ls': '--'},
    {'key': 'tr_ppi',   'label': '生産者物価テイラールール（TR-PPI）',      'color': 'red',          'ls': '--'},
    {'key': 'flex_ppi', 'label': '完全伸縮価格（IT-PPI）',          'color': 'gray',         'ls': '-.'}
]

# ショックの種類定義（ログ用とグラフ用を分離）
shocks = [
    {
        'id': 'beta_shock', 
        'log_name': 'beta ショック（需要ショック）', 
        'plot_name': r' $\beta$ ショック（需要ショック）', 
        'scenario': 'beta'
    },
    {
        'id': 'a_shock',    
        'log_name': 'a ショック（供給ショック）', 
        'plot_name': r' $a$ ショック（供給ショック）', 
        'scenario': 'a'
    }
]

# 現在のスクリプトのある場所をベースディレクトリとする
base_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = base_dir

# グラフを保存するルートディレクトリ
output_root_dir = os.path.join(base_dir, 'comparison_graphs')

# プロットする期間（100期まで）
plot_periods = 100
x_axis = np.arange(plot_periods)

# グラフ化する変数のリスト (MATLAB側の var_list と完全に同期)
var_list = [ 
    'e_slash_star',
    'lambda_H', 'lambda_H_slash_star',
    'lambda_F_slash_star', 'lambda_F',             
    'utility_H_with_delta', 'utility_F_with_delta',
    'utility_H_without_delta', 'utility_F_without_delta', 
    'gamma_H',
    'p_H', 'p_F_star',
    'i_H_notional', 'i_H', 'i_F', 
    'c_H_W', 'c_F_W', 'c_H_H', 'c_F_H', 'c_H_F', 'c_F_F',
    'y_H', 'y_H_potential', 'y_F',
    'l_H', 'l_F', 
    'p_H_W', 'p_F_W_star',
    'p_H_W_bar', 'p_F_W_star_bar', 
    'p_H_bar', 'p_F_star_bar', 
    'p_H_bar_y_H', 'p_H_W_c_H_W', 
    'pi_H', 'pi_F_star', 
    'pi_H_W', 'pi_F_W_star', 
    'p_H_tilde', 'p_F_star_tilde',
    'v_H', 'v_F',
    'w_H', 'w_F', 
    'p_F_bar', 'p_H_star_bar',
    't_H', 't_F',
    'b_H', 'b_F_star', 
    'beta_H', 'beta_F',
    'a_H', 'a_F',
    'tau_H', 'tau_F',
    'chi_H',
    'eps_beta_H', 'eps_beta_F',
    'eps_a_H', 'eps_a_F',
    'eps_tau_H', 'eps_tau_F',
    'eps_i_H', 'eps_i_F',
    'eps_chi_H', 'eps_y_H_potential',
    'Delta_H', 'Delta_F'
]

# LaTeXのタイトルをマッピング
title_map = {
    'e_slash_star': r'名目為替レート（ $e^{/*}$ ）',
    'lambda_H': r'自国の所得の限界効用（ $\lambda^H$ ）',
    'lambda_H_slash_star': r'自国の外国通貨建て所得の限界効用（ $\lambda^{H/*}$ ）',
    'lambda_F_slash_star': r'外国の所得の限界効用（ $\lambda^{F/*}$ ）',
    'lambda_F': r'外国の自国通貨建て所得の限界効用（ $\lambda^F$ ）',
    'utility_H_with_delta': r'自国の効用（ $\Delta$ 有り）',
    'utility_F_with_delta': r'外国の効用（ $\Delta$ 有り）', 
    'utility_H_without_delta': r'自国の効用（ $\Delta$ 無し）', 
    'utility_F_without_delta': r'外国の効用（ $\Delta$ 無し）', 
    'gamma_H': r'自国金融政策の目標未達分の累積（ $\gamma^H$ ）',
    'i_H_notional': r'自国の潜在利子率（ $i^{H,notional}$ ）',
    'i_H': r'自国の利子率（ $i^H$ ）',
    'i_F': r'外国の利子率（ $i^F$ ）',
    'c_H_W': r'自国の総消費指数（ $c^{H \to W}$ ）',
    'c_F_W': r'外国の総消費指数（ $c^{F \to W}$ ）',
    'c_H_H': r'自国の自国財消費指数（ $c^{H \to H}$ ）',
    'c_F_H': r'外国の自国財消費指数（ $c^{F \to H}$ ）',
    'c_H_F': r'自国の外国財消費指数（ $c^{H \to F}$ ）',
    'c_F_F': r'外国の外国財消費指数（ $c^{F \to F}$ ）',
    'y_H': r'自国の生産（ $y^H$ ）',
    'y_H_potential': r'自国の潜在産出（ $y^{H,pot}$ ）',
    'y_F': r'外国の生産（ $y^F$ ）',
    'l_H': r'自国の労働（ $l^H$ ）',
    'l_F': r'外国の労働（ $l^F$ ）',
    'p_H': r'自国の生産者物価指数（PPI）（ $p^H$ ）',
    'p_H_star_bar': r'自国の外国通貨建て生産者物価指数（PPI）（ $\bar{p}^{H*}$ ）',
    'p_F_star': r'外国の生産者物価指数（PPI）（ $p^{F*}$ ）',
    'p_F_bar': r'外国の自国通貨建て生産者物価指数（PPI）（ $\bar{p}^F$ ）',
    'p_H_bar': r'自国の正規化された生産者物価指数（PPI）（ $\bar{p}^H$ ）',
    'p_F_star_bar': r'外国の正規化された生産者物価指数（PPI）（ $\bar{p}^{F*}$ ）',
    'p_H_W': r'自国の消費者物価指数（CPI）（ $p^{H \to W}$ ）',
    'p_H_W_bar': r'自国の正規化された消費者物価指数（CPI）（ $\bar{p}^{H \to W}$ ）',
    'p_F_W_star': r'外国の消費者物価指数（CPI）（ $p^{F \to W*}$ ）',
    'p_F_W_star_bar': r'外国の正規化された消費者物価指数（CPI）（ $\bar{p}^{F \to W*}$ ）',
    'p_H_bar_y_H': r'自国の名目GDP（ $\bar{p}^H y^H$ ）',
    'p_H_W_c_H_W': r'自国の名目総消費（ $p^{H \to W} c^{H \to W}$ ）',
    'p_H_tilde': r'自国の最適生産者物価（ $\tilde{p}^H$ ）',
    'p_F_star_tilde': r'外国の最適生産者物価（ $\tilde{p}^{F*}$ ）',
    'v_H': r'自国の最適生産者物価の補助変数1（ $v^H$ ）',
    'v_F': r'外国の最適生産者物価の補助変数1（ $v^F$ ）',
    'w_H': r'自国の最適生産者物価の補助変数2（ $w^H$ ）',
    'w_F': r'外国の最適生産者物価の補助変数2（ $w^F$ ）',
    'pi_H': r'自国の生産者物価指数（PPI）を用いたグロス・インフレ率（ $\pi^H$ ）',
    'pi_H_W': r'自国の消費者物価指数（CPI）を用いたグロス・インフレ率（ $\pi^{H \to W}$ ）',
    'pi_F_star': r'外国の生産者物価指数（PPI）を用いたグロス・インフレ率（ $\pi^{F*}$ ）',
    'pi_F_W_star': r'外国の消費者物価（CPI）を用いたグロス・インフレ率（ $\pi^{F \to W*}$ ）',    
    't_H': r'自国の一括移転（ $t^H$ ）',
    't_F': r'外国の一括移転（ $t^F$ ）',
    'b_H': r'自国の対外純資産（ $b^H$ ）',
    'b_F_star': r'外国の対外純資産（ $b^{F*}$ ）',
    'beta_H': r'自国の主観的割引因子（ $\beta^H$ ）',
    'beta_F': r'外国の主観的割引因子（ $\beta^F$ ）',
    'a_H': r'自国の生産性（ $a^H$ ）',
    'a_F': r'外国の生産性（ $a^F$ ）',
    'tau_H': r'自国の税率（ $\tau^H$ ）',
    'tau_F': r'外国の税率（ $\tau^F$ ）',
    'chi_H': r'自国金融政策の目標パス（ $\chi^H$ ）',
    'eps_beta_H': r'需要ショック（ $\epsilon^{\beta^H}$ ）',
    'eps_beta_F': r'外国環境需要ショック（ $\epsilon^{\beta^F}$ ）',
    'eps_a_H': r'生産性ショック（ $\epsilon^{a^H}$ ）',
    'eps_a_F': r'外国生産性ショック（ $\epsilon^{a^F}$ ）',
    'eps_tau_H': r'税率ショック（ $\epsilon^{\tau^H}$ ）',
    'eps_tau_F': r'外国税率ショック（ $\epsilon^{tau^F}$ ）',
    'eps_i_H': r'利子率ショック（ $\epsilon^{i^H}$ ）',
    'eps_i_F': r'外国利子率ショック（ $\epsilon^{i^F}$ ）',
    'eps_chi_H': r'目標パスショック（ $\epsilon^{\chi^H}$ ）',
    'eps_y_H_potential': r'潜在産出量ショック（ $\epsilon^{y^{H,pot}}$ ）',
    'Delta_H': r'自国の価格分散（ $\Delta^H$ ）',
    'Delta_F': r'外国の価格分散（ $\Delta^F$ ）',
}


# --- メイン処理ループ (ショックごと) ---

for shock in shocks:
    shock_id = shock['id']
    shock_log_name = shock['log_name']
    shock_plot_name = shock['plot_name']
    shock_scenario = shock['scenario']
    
    print(f"\n========== {shock_log_name} の比較処理を開始します ==========")
    
    current_output_dir = os.path.join(output_root_dir, shock_id)
    if not os.path.exists(current_output_dir):
        os.makedirs(current_output_dir)

    # --- 2. データの読み込み ---
    results = {}
    welfare_values = {}
    welfare_no_delta_values = {} 

    for policy in policies:
        key = policy['key']
        
        # パス設定の分岐
        if key == 'flex_ppi':
            mat_rel_path = os.path.join('polt_project', 'results', 'stage1', shock_id, 'data_for_plotting.mat')
        elif key in ['polt', 'tr_ppi', 'tr_cpi']:
            mat_rel_path = os.path.join(f'{key}_project', 'results', 'stage2', shock_id, 'data_for_plotting.mat')
        else:
            mat_rel_path = os.path.join(f'{key}_project', 'results', shock_id, 'data_for_plotting.mat')
            
        mat_path = os.path.normpath(os.path.join(root_dir, mat_rel_path))
        results[key] = {}
        
        if not os.path.exists(mat_path):
            print(f"  [警告] ファイルが見つかりません： {mat_path}")
            results[key] = None
            continue

        try:
            data = scipy.io.loadmat(mat_path)
            
            for var in var_list:
                # 命名規則に基づき読み込み
                if var == 'utility_H_with_delta': mat_v = 'utility_H_with_delta_pie'
                elif var == 'utility_F_with_delta': mat_v = 'utility_F_with_delta_pie'
                elif var == 'utility_H_without_delta': mat_v = 'utility_H_without_delta_pie'
                elif var == 'utility_F_without_delta': mat_v = 'utility_F_without_delta_pie'
                else: mat_v = f"{var}_pie"
                
                if mat_v in data:
                    results[key][var] = data[mat_v].flatten()

            # 厚生データの読み込み (スカラー値)
            if 'welfare_H_with_delta_pie' in data:
                 welfare_values[f'{key}_H'] = data['welfare_H_with_delta_pie'].flatten()[0]
            if 'welfare_F_with_delta_pie' in data:
                 welfare_values[f'{key}_F'] = data['welfare_F_with_delta_pie'].flatten()[0]
            if 'welfare_H_without_delta_pie' in data:
                 welfare_no_delta_values[f'{key}_H'] = data['welfare_H_without_delta_pie'].flatten()[0]
            if 'welfare_F_without_delta_pie' in data:
                 welfare_no_delta_values[f'{key}_F'] = data['welfare_F_without_delta_pie'].flatten()[0]

        except Exception as e:
            print(f"  [エラー] 読み込み中にエラーが発生 ({key})： {e}")
            results[key] = None

    # --- 3. グラフの描画ループ ---
    print(f"  グラフを作成中...")
    
    count_created = 0
    for var_name in var_list:
        
        is_h_utility = var_name == 'utility_H_with_delta'
        is_f_utility = var_name == 'utility_F_with_delta'
        is_h_without_delta_utility = var_name == 'utility_H_without_delta'
        is_f_without_delta_utility = var_name == 'utility_F_without_delta'

        has_data = False
        for policy in policies:
            key = policy['key']
            if results.get(key) is not None and var_name in results[key]:
                has_data = True
                break
        
        if not has_data:
            continue

        fig, ax = plt.subplots(figsize=(12, 8)) 
        
        title_text = title_map.get(var_name, var_name)
        ax.set_title(f'{title_text}の比較（{shock_plot_name}）', fontsize=16)

        legend_handles = []
        
        for policy in policies:
            key = policy['key']
            if results.get(key) is not None and var_name in results[key]:
                data_series = results[key][var_name]
                plot_len = min(plot_periods, len(data_series))
                
                line, = ax.plot(x_axis[:plot_len], data_series[:plot_len], 
                                color=policy['color'], linestyle=policy['ls'], 
                                label=policy['label'], linewidth=2, alpha=0.8)
                legend_handles.append(line)

        ax.set_xlabel('期（四半期）', fontsize=12)
        ax.set_ylabel('定常状態からの乖離', fontsize=12)
        ax.grid(True)
        ax.set_xlim(0, plot_periods)
        ax.legend(handles=legend_handles, loc='best', fontsize=9)


        # --- 厚生値のテキスト表示設定 ---
        is_welfare_plot = is_h_utility or is_f_utility or is_h_without_delta_utility or is_f_without_delta_utility
        
        if is_welfare_plot:
            is_home = is_h_utility or is_h_without_delta_utility
            
            f_size = 20
            pos_x = 0.98
            pos_y = 0.02
            v_align = 'bottom'
            h_align = 'right'
            # テキスト内部の揃え：デフォルトは右揃え
            m_align = 'right'

            if shock_id == 'beta_shock':
                if is_home:
                    # ① Betaショック × 自国(H)
                    f_size = 22
                    pos_x = 0.98
                    pos_y = 0.02
                    v_align = 'bottom'
                    h_align = 'right'
                else:
                    # ② Betaショック × 外国(F)
                    # 【修正点】サイズを 13 から 14 へ大きくし、内部を右揃えに設定
                    f_size = 14  
                    pos_x = 0.02
                    pos_y = 0.02
                    v_align = 'bottom'
                    h_align = 'left'
                    m_align = 'right' 
            
            elif shock_id == 'a_shock':
                if is_home:
                    # ③ Aショック × 自国(H)
                    f_size = 22
                    pos_x = 0.98
                    pos_y = 0.02
                    v_align = 'bottom'
                    h_align = 'right'
                else:
                    # ④ Aショック × 外国(F)
                    # 【修正点】サイズを 22 にし、内部を右揃えに設定
                    f_size = 22  
                    pos_x = 0.98
                    pos_y = 0.02
                    v_align = 'bottom'
                    h_align = 'right'
                    m_align = 'right'

            # 表示データの準備
            if is_home:
                w_dict = welfare_values if is_h_utility else welfare_no_delta_values
                delta_txt = r'$\Delta$ 有り' if is_h_utility else r'$\Delta$ 無し'
                suffix = '_H'
            else:
                w_dict = welfare_values if is_f_utility else welfare_no_delta_values
                delta_txt = r'$\Delta$ 有り' if is_f_utility else r'$\Delta$ 無し'
                suffix = '_F'
            
            shock_name_jp = r'$\beta$ ショック（需要ショック）' if shock_id == 'beta_shock' else r'$a$ ショック（供給ショック）'
            title_line = f"厚生{suffix[1]}（{delta_txt}, {shock_name_jp}）："
            
            # Y軸の拡張
            ymin, ymax = ax.get_ylim()
            if pos_y < 0.5:
                y_range = ymax - ymin
                ax.set_ylim(ymin - y_range * 0.6, ymax) 
            
            # 厚生値テキストの作成
            text_lines = [title_line]
            for policy in policies:
                key = policy['key']
                w_key = f"{key}{suffix}"
                if w_key in w_dict:
                    val = w_dict[w_key]
                    text_lines.append(f"{policy['label']}： {val:.4f}")
            
            text_str = "\n".join(text_lines)
            props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
            
            # multialignment=m_align を指定して内部テキストを揃える
            ax.text(pos_x, pos_y, text_str, transform=ax.transAxes, fontsize=f_size, 
                    verticalalignment=v_align, horizontalalignment=h_align, 
                    multialignment=m_align, bbox=props)

        else:
            # その他の変数のパディング調整
            ymin, ymax = ax.get_ylim()
            padding = (ymax - ymin) * 0.05 if ymax != ymin else 0.1
            ax.set_ylim(ymin - padding, ymax + padding)

        # 保存処理
        plt.tight_layout() 
        file_name = f'{shock_id}_{var_name}.png'
        save_path = os.path.join(current_output_dir, file_name)
        plt.savefig(save_path) 
        plt.close(fig)
        count_created += 1

    # =========================================================================
    # --- POLT限定： 生産追従確認グラフ (y_H vs y_H_potential) ---
    # =========================================================================
    if results.get('polt') is not None and 'y_H' in results['polt'] and 'y_H_potential' in results['polt']:
        print(f"  POLT専用の生産追従グラフを作成中...")
        fig_polt, ax_polt = plt.subplots(figsize=(12, 8))
        y_h_data = results['polt']['y_H'][:plot_periods]
        y_pot_data = results['polt']['y_H_potential'][:plot_periods]
        
        ax_polt.plot(x_axis[:len(y_h_data)], y_h_data, color='black', label=r'実際の生産（ $y^H$ ）', linewidth=2.5)
        ax_polt.plot(x_axis[:len(y_pot_data)], y_pot_data, color='red', linestyle='--', label=r'潜在生産（ $y^{H,pot}$ ）', linewidth=2.5)
        
        ax_polt.set_title(f'潜在生産水準目標（POLT）における目標追従の確認\n（{shock_plot_name}）', fontsize=16)
        ax_polt.set_xlabel('期（四半期）', fontsize=12)
        ax_polt.set_ylabel('定常状態からの乖離', fontsize=12)
        ax_polt.grid(True)
        ax_polt.set_xlim(0, plot_periods)
        ax_polt.legend(loc='best', fontsize=12)
        
        polt_specific_path = os.path.join(current_output_dir, f'{shock_id}_polt_y_H_and_y_H_potential.png')
        plt.tight_layout()
        plt.savefig(polt_specific_path)
        plt.close(fig_polt)

    print(f"  -> 合計 {count_created} 個の比較グラフを保存しました： {current_output_dir}")

print("\n全ての処理が完了しました。")