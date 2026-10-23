import os
from datetime import datetime

# --- 設定 ---
TARGET_EXTENSIONS = {'.m', '.mod', '.inc'}
IGNORE_DIRS = {'results', 'memos', 'project_all_code'} # 出力先フォルダや除外するフォルダ

def create_all_consolidated_files():
    # スクリプトがある場所（projectsフォルダ）を基準とする
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 1. 出力用フォルダを作成 (projects/project_all_code/)
    output_dir = os.path.join(script_dir, 'project_all_code')
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. projectsフォルダ直下にある「_project」で終わるフォルダをすべて取得
    project_folders = [d for d in os.listdir(script_dir) 
                       if os.path.isdir(os.path.join(script_dir, d)) and d.endswith('_project')]
    
    if not project_folders:
        print("対象となるプロジェクトフォルダ（'_project'で終わるフォルダ）が見つかりません。")
        return

    print(f"発見されたプロジェクト数: {len(project_folders)}")

    # 各プロジェクトに対して処理をループ
    for project_name in project_folders:
        target_root = os.path.join(script_dir, project_name)
        
        # 3. フォルダ名からポリシー名（例: clt_project -> clt）を抽出
        policy_key = project_name.replace('_project', '')
        
        # 新しい出力ファイル名を生成してパスを結合
        output_filename = f"project_all_code_{policy_key}.txt"
        output_path = os.path.join(output_dir, output_filename)
        
        print(f"\nScanning target: {target_root}")
        print(f"Output destination: {output_path}")

        with open(output_path, 'w', encoding='utf-8') as outfile:
            outfile.write(f"--- Consolidated Project Files: {project_name} ---\n")
            outfile.write(f"Target Directory: {target_root}\n")
            outfile.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            for root, dirs, files in os.walk(target_root):
                # --- フォルダの除外ロジック (ネスト対応) ---
                path_parts = os.path.normpath(root).split(os.sep)
                
                # 1. 'src' がパスに含まれる場合
                if 'src' in path_parts:
                    src_idx = path_parts.index('src')
                    # 'src' の直下より深い階層をチェック
                    if len(path_parts) > src_idx + 1:
                        current_folder = path_parts[src_idx + 1]
                        # 'src/stage...' までは許可するが、その下の階層(src/stage1/foo)は除外
                        if not current_folder.startswith('stage') or len(path_parts) > src_idx + 2:
                            dirs[:] = []
                            continue

                # 2. ルート直下に 'stage' がある場合 (src を介さないケース)
                elif any(p.startswith('stage') for p in path_parts):
                    stage_idx = next(i for i, p in enumerate(path_parts) if p.startswith('stage'))
                    if len(path_parts) > stage_idx + 1:
                        dirs[:] = []
                        continue

                # 3. 汎用的な除外設定
                dirs[:] = [d for d in dirs if d not in IGNORE_DIRS and not d.startswith('+')]

                # --- ファイルの処理 ---
                for file in sorted(files):
                    file_path = os.path.join(root, file)
                    if file_path == os.path.abspath(__file__) or file_path == os.path.abspath(output_path) or file.startswith('.'):
                        continue

                    if os.path.splitext(file)[1].lower() in TARGET_EXTENSIONS:
                        relative_path = os.path.relpath(file_path, target_root)
                        print(f"  > Including: {relative_path}")
                        outfile.write("\n" + "="*80 + "\n" + f"FILE: {relative_path}\n" + "="*80 + "\n\n")
                        try:
                            with open(file_path, 'r', encoding='utf-8') as infile:
                                outfile.write(infile.read())
                        except Exception as e:
                            outfile.write(f"<< Error reading file: {str(e)} >>\n")
                        outfile.write("\n")

    print(f"\nDone! 全てのプロジェクトのコードが {output_dir} に出力されました。")

if __name__ == "__main__":
    create_all_consolidated_files()