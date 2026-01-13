% =========================================================================
% env_setup.m
% 目的: Dynare および OccBin ツールキットのパス設定
%
% --- ソフトウェア互換性ガイド (Software Compatibility Guide) ---
% Dynare と Octave の組み合わせには厳密な互換性がある。
% 自身のインストールしているバージョンが以下に適合しているか確認すること。
%
% [互換性対応表]
% Dynare Version | Supported Octave (Windows) | Required Matlab
% -------------------------------------------------------------
% 6.4            | 10.2                       | 9.5 (R2018b)
% 6.3            | 9.4                        | 9.5 (R2018b)
% 6.2            | 9.2                        | 9.5 (R2018b)
% 6.1            | 9.1                        | 9.5 (R2018b)
% 6.0            | 8.4.0                      | 9.5 (R2018b)
% 5.5            | 8.3.0                      | 8.3 (R2014a)
% 5.4            | 8.1                        | 8.3 (R2014a)
% 5.3            | 7.3.0                      | 8.3 (R2014a)
% 5.0 - 5.2      | 6.4.0                      | 8.3 (R2014a)
% 4.6.4          | 6.2.0                      | 7.9 (R2009b)
% 4.6.0 - 4.6.3  | 5.2                        | 7.9 (R2009b)
% 4.5.x          | 4.4.1                      | 7.5 (R2007b)
%
% 詳細は Dynare リリースノート、および以下の公式フォーラムを参照：
% https://forum.dynare.org/t/compatibility-of-softwares/25986/2
% =========================================================================

% --- 1. Dynare へのパス設定 ---
% 自身の PC にインストールされている Dynare の場所に合わせて書き換える
addpath('C:\dynare\6.3\matlab');

% --- 2. OccBin ツールキットへのパス設定 ---
% 展開した Occbin_update-master 内の toolkit_files を指定する
addpath('C:\dynare\Occbin_update-master\toolkit_files');

% --- 3. 設定確認メッセージ ---
fprintf('\nSetup complete. You can now run your project scripts.\n');