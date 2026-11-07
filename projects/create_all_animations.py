"""
【使用方法】
1. このスクリプト（create_all_animations.py）を projects ディレクトリ直下に配置してください。
2. ターミナルを開き任意のカレントディレクトリで実行してください。
   実行コマンド: python path/to/create_all_animations.py
"""

import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Windows標準の日本語フォントを指定（マイナス記号の文字化けも防止）
plt.rcParams["font.family"] = "MS Gothic"
plt.rcParams["axes.unicode_minus"] = False

# 比較する政策のリスト定義
policies = [
    {'key': 'clt',      'label': '総消費水準目標（CLT）'},
    {'key': 'cplt',     'label': '消費者物価水準目標（CPLT）'},
    {'key': 'it_cpi',   'label': '消費者物価インフレ目標（IT-CPI）'},
    {'key': 'it_ppi',   'label': '生産者物価インフレ目標（IT-PPI）'},
    {'key': 'nclt',     'label': '名目総消費水準目標（NCLT）'},
    {'key': 'ngdplt',   'label': '名目GDP水準目標（NGDPLT）'},
    {'key': 'olt',      'label': '生産水準目標（OLT）'},
    {'key': 'polt',     'label': '潜在生産水準目標（POLT）'},
    {'key': 'pplt',     'label': '生産者物価水準目標（PPLT）'},
    {'key': 'tr_cpi',   'label': '消費者物価テイラールール（TR-CPI）'},
    {'key': 'tr_ppi',   'label': '生産者物価テイラールール（TR-PPI）'},
    {'key': 'flex_ppi', 'label': '完全伸縮価格（IT-PPI）'}
]

shocks = [
    {'id': 'beta_shock', 'log_name': 'beta ショック（需要ショック）', 'plot_name': 'β ショック（需要ショック）'},
    {'id': 'a_shock',    'log_name': 'a ショック（供給ショック）',    'plot_name': 'a ショック（供給ショック）'}
]

base_dir = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------
# まとめて保存するための animations フォルダを作成
# ---------------------------------------------------
animations_dir = os.path.join(base_dir, 'animations')
os.makedirs(animations_dir, exist_ok=True)

# 経済モデルのパラメータ設定
kappa_H = 0.05
alpha_H = 0.85
MAX_FRAMES = 100

# ★描画対象データの選択: '_pie' (ZLB考慮・OccBin解) または '_lin' (ZLB無視・線形解)
DATA_SUFFIX = '_pie'

# 曲線の傾き計算
SLOPE_YS = 2 * kappa_H
SLOPE_YD = -1.0
SLOPE_CS = (2 * kappa_H) / alpha_H + (1 - alpha_H) / alpha_H
SLOPE_CD = -1.0


def get_auto_limits(data, num_frames, padding=0.2):
    """データ自身の最小値と最大値を取得し、マージンを加える。
       極端なズームイン（1e-9等）を防ぐため、最低でも 0.005 の余白を確保する。"""
    d_min, d_max = np.min(data[:num_frames]), np.max(data[:num_frames])
    margin = max((d_max - d_min) * padding, 0.005)
    return (d_min - margin, d_max + margin)


def resolve_csv_path(base_dir, key, shock_id):
    """新しいファイル命名規則に基づいてCSVファイルのパスを解決する（フォールバック付き）"""
    if key == 'flex_ppi':
        csv_filename = f'polt_stage1_data_for_plotting_{shock_id}.csv'
        primary_path = os.path.join(base_dir, 'polt_project', 'results', 'stage1', shock_id, csv_filename)
        return primary_path, os.path.dirname(primary_path)

    csv_filename = f'{key}_data_for_plotting_{shock_id}.csv'
    if key in ['polt', 'tr_ppi', 'tr_cpi']:
        primary_dir = os.path.join(base_dir, f'{key}_project', 'results', 'stage2', shock_id)
    else:
        primary_dir = os.path.join(base_dir, f'{key}_project', 'results', shock_id)

    primary_path = os.path.join(primary_dir, csv_filename)
    if os.path.exists(primary_path):
        return primary_path, primary_dir

    # フォールバック: projects/data_for_plotting フォルダ内の集約CSVも確認
    fallback_path = os.path.join(base_dir, 'data_for_plotting', csv_filename)
    if os.path.exists(fallback_path):
        os.makedirs(primary_dir, exist_ok=True)
        return fallback_path, primary_dir

    return primary_path, primary_dir


# =========================================================================
# --- メイン処理ループ ---
# =========================================================================

for shock in shocks:
    shock_id = shock['id']
    print(f"\n========== {shock['log_name']} のアニメーション作成を開始します ==========")

    for policy in policies:
        key = policy['key']

        csv_path, output_dir = resolve_csv_path(base_dir, key, shock_id)
        csv_path = os.path.normpath(csv_path)

        if not os.path.exists(csv_path):
            print(f"  [スキップ] ファイルが見つかりません： {csv_path}")
            continue

        try:
            df = pd.read_csv(csv_path, skipinitialspace=True)
            df.columns = [col.strip() for col in df.columns]

            # Dynareの微小な計算誤差（1e-15等）を切り捨てて「0」と認識させる
            df = df.round(10)

            y = df[f'y_H{DATA_SUFFIX}'].values
            p = df[f'p_H{DATA_SUFFIX}'].values
            c = df[f'c_H_W{DATA_SUFFIX}'].values
            pw = df[f'p_H_W{DATA_SUFFIX}'].values

            num_frames = min(MAX_FRAMES, len(df))

            # 丸め処理後の実際のデータの変動幅をもとに各軸の描画範囲を自動決定する
            xlim_y = get_auto_limits(y, num_frames)
            ylim_p = get_auto_limits(p, num_frames)
            xlim_c = get_auto_limits(c, num_frames)
            ylim_pw = get_auto_limits(pw, num_frames)

            # 図の初期化
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

            line_ys, = ax1.plot([], [], color='red', label='YS Curve')
            line_yd, = ax1.plot([], [], color='blue', label='YD Curve')
            point_y_p, = ax1.plot([], [], 'ko', markersize=8, label='Equilibrium')
            trace_y_p, = ax1.plot([], [], color='gray', linestyle='--', alpha=0.5)

            line_cs, = ax2.plot([], [], color='red', label='CS Curve')
            line_cd, = ax2.plot([], [], color='blue', label='CD Curve')
            point_c_pw, = ax2.plot([], [], 'ko', markersize=8, label='Equilibrium')
            trace_c_pw, = ax2.plot([], [], color='gray', linestyle='--', alpha=0.5)

            # 軸の設定
            ax1.set_xlim(xlim_y)
            ax1.set_ylim(ylim_p)
            ax1.set_title('生産・生産者物価平面 (YS-YD Curve)', fontsize=14)
            ax1.set_xlabel(r'生産の乖離 ($y_t^H$)', fontsize=12)
            ax1.set_ylabel(r'PPIの乖離 ($p_t^H$)', fontsize=12)
            ax1.axhline(0, color='black', linewidth=0.5)
            ax1.axvline(0, color='black', linewidth=0.5)
            ax1.grid(True, linestyle=':', alpha=0.6)
            ax1.legend()

            ax2.set_xlim(xlim_c)
            ax2.set_ylim(ylim_pw)
            ax2.set_title('消費・消費者物価平面 (CS-CD Curve)', fontsize=14)
            ax2.set_xlabel(r'総消費指数の乖離 ($c_t^{H \to W}$)', fontsize=12)
            ax2.set_ylabel(r'CPIの乖離 ($p_t^{H \to W}$)', fontsize=12)
            ax2.axhline(0, color='black', linewidth=0.5)
            ax2.axvline(0, color='black', linewidth=0.5)
            ax2.grid(True, linestyle=':', alpha=0.6)
            ax2.legend()

            fig.tight_layout(rect=[0, 0, 1, 0.95])

            def update(frame):
                curr_y, curr_p = y[frame], p[frame]

                # 直線を描画するためのX座標の配列
                x_vals_ax1 = np.linspace(xlim_y[0], xlim_y[1], 100)
                ys_vals = SLOPE_YS * (x_vals_ax1 - curr_y) + curr_p
                yd_vals = SLOPE_YD * (x_vals_ax1 - curr_y) + curr_p

                line_ys.set_data(x_vals_ax1, ys_vals)
                line_yd.set_data(x_vals_ax1, yd_vals)
                point_y_p.set_data([curr_y], [curr_p])
                trace_y_p.set_data(y[:frame+1], p[:frame+1])

                curr_c, curr_pw = c[frame], pw[frame]

                x_vals_ax2 = np.linspace(xlim_c[0], xlim_c[1], 100)
                cs_vals = SLOPE_CS * (x_vals_ax2 - curr_c) + curr_pw
                cd_vals = SLOPE_CD * (x_vals_ax2 - curr_c) + curr_pw

                line_cs.set_data(x_vals_ax2, cs_vals)
                line_cd.set_data(x_vals_ax2, cd_vals)
                point_c_pw.set_data([curr_c], [curr_pw])
                trace_c_pw.set_data(c[:frame+1], pw[:frame+1])

                fig.suptitle(f"{policy['label']} ({shock['plot_name']}) - Time: $t={frame}$", fontsize=16)

                return line_ys, line_yd, point_y_p, trace_y_p, line_cs, line_cd, point_c_pw, trace_c_pw

            # 各プロジェクトのresultsフォルダ内にも統一命名規則で保存
            indiv_gif_name = f"{key}_animation_{shock_id}.gif"
            output_filename = os.path.join(output_dir, indiv_gif_name)

            print(f"  [{key}] アニメーションを作成中... -> {output_filename}")

            ani = FuncAnimation(fig, update, frames=num_frames, interval=100, blit=False)
            ani.save(output_filename, writer='pillow', fps=10)

            plt.close(fig)

            # --- animationsフォルダへのコピー処理 ---
            copy_dest = os.path.join(animations_dir, indiv_gif_name)
            shutil.copy2(output_filename, copy_dest)
            print(f"  [{key}] 完了: 一覧用フォルダにコピーしました -> {copy_dest}")

        except Exception as e:
            print(f"  [エラー] {key} の処理中にエラーが発生しました： {e}")

print("\n全ての処理が完了しました。")