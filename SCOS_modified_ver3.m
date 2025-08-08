%% SCOS/DCS 統合評価シミュレーションコード (MCXベース) - v10.0
%
% 概要:
% 複数の血流レベル（基準の1～5倍）を設定し、各血流指標（BFI）の
% 「相対的な変化」の追従精度を評価する。
% 光子追跡(MCX)は最初に1回のみ実行し、計算を効率化。
% SDSごとに個別のグラフを出力し、各手法の線形性を評価する。
% 最終更新日: 2025/06/27
close all; clear; clc;

%% =========================================================================
% ステージ1: MCXによる光伝播シミュレーション (血流に非依存)
% =========================================================================
fprintf('--- ステージ1: MCXによる光伝播シミュレーション開始 (1回のみ実行) ---\n');
addpath(genpath('E:\モンテカルロシミュレーション/MCXStudio-linux64-v2025/MCXStudio/MATLAB/iso2mesh'));
addpath(genpath('E:\モンテカルロシミュレーション/MCXStudio-linux64-v2025'));
addpath('E:\モンテカルロシミュレーション/FW__code_request');
cfg = struct(); cfg.nphoton = 5e6; cfg.unitinmm = 0.5; cfg.seed = randi([1 2^31-1]);
vol_dim = [120, 120, 60]; cfg.vol = uint8(zeros(vol_dim));
scalp_th = 3; skull_th = 5;
idx_z_scalp = 1:ceil(scalp_th / cfg.unitinmm);
idx_z_skull = (1+ceil(scalp_th/cfg.unitinmm)):ceil((scalp_th+skull_th)/cfg.unitinmm);
idx_z_brain = (1+ceil((scalp_th+skull_th)/cfg.unitinmm)):vol_dim(3);
cfg.vol(:, :, idx_z_scalp) = 1; cfg.vol(:, :, idx_z_skull) = 2; cfg.vol(:, :, idx_z_brain) = 3;
cfg.prop = [0.0,0.0,1.0,1.0; 0.016,8.0,0.9,1.37; 0.01,7.0,0.9,1.37; 0.02,10.0,0.9,1.37];
cfg.srcpos = [vol_dim(1)/2, vol_dim(2)/2, 1]; cfg.srcdir = [0, 0, 1];
sds_mm = [10, 20, 30]; sds_grid = sds_mm / cfg.unitinmm;
detector_radius_mm = 2; detector_radius_grid = detector_radius_mm / cfg.unitinmm;
cfg.detpos = [];
for i = 1:length(sds_grid)
    cfg.detpos(i, :) = [cfg.srcpos(1) + sds_grid(i), cfg.srcpos(2), 1, detector_radius_grid];
end
num_detectors = size(cfg.detpos, 1);
cfg.tstart = 0; cfg.tend = 5e-9; cfg.tstep = 5e-9;
cfg.savedetflag = 'dsp';
cfg.isreflect = 1; cfg.autopilot = 1; cfg.gpuid = 1;
[flux, detpos] = mcxlab(cfg);
fprintf('MCXシミュレーション完了。検出された総光子数: %d\n', length(detpos.detid));
fprintf('-------------------------------------------------------\n\n');

%% =========================================================================
% ステージ2: 複数のBFIレベルにおける解析
% =========================================================================
fprintf('--- ステージ2: 複数のBFIレベルにおける解析開始 ---\n');
% --- 解析用パラメータ ---
lambda = 785e-6; n_brain = cfg.prop(4,4); k0 = 2*pi*n_brain/lambda;
mus_prime_brain = cfg.prop(4, 2) * (1 - cfg.prop(4, 3));
T = 10; beta = 0.5;

% --- テストするBFIのリストを定義 ---
alpha_brain = 0.5;
baseline_Db = 1e-6; % 基準となる拡散係数 [mm^2/s]
relative_levels = [1, 2, 3, 4, 5]; % 基準に対する倍率
db_list = baseline_Db * relative_levels;
bfi_list = alpha_brain * db_list; % テストするBFIのリスト

% --- 結果保存用の構造体を初期化 ---
all_results = struct();

% --- BFIレベルごとにループ ---
for bfi_idx = 1:length(bfi_list)
    current_true_BFI = bfi_list(bfi_idx);
    fprintf('\n===== BFIレベル %d/%d (True BFI = %.2e) の解析中 =====\n', bfi_idx, length(bfi_list), current_true_BFI);
    
    % --- SDSごとにループ ---
    for sds_idx = 1:num_detectors
        sds_label = sprintf('sds_%dmm', sds_mm(sds_idx));
        fprintf('>>> 検出器 (SDS = %d mm) の処理中...\n', sds_mm(sds_idx));
        
        idx_all_photons = find(detpos.detid == sds_idx);
        if isempty(idx_all_photons); continue; end
        
        % g1(τ)の計算
        path_brain_all = detpos.ppath(idx_all_photons, 3);
        idx_deep_photons = find(path_brain_all > 0);
        idx_shallow_photons = find(path_brain_all == 0);
        N_deep = length(idx_deep_photons); N_shallow = length(idx_shallow_photons); N_total = N_deep + N_shallow;
        path_brain_deep = path_brain_all(idx_deep_photons);
        
        tau_vector_ms = linspace(0, T, 200); tau_vector_s = tau_vector_ms / 1000;
        g1_vector = zeros(size(tau_vector_s));
        if N_deep > 0
            g1_deep_decay = zeros(size(tau_vector_s));
            for j = 1:length(tau_vector_s)
                decay_term = exp(-2 * k0^2 * current_true_BFI * (mus_prime_brain * path_brain_deep) * tau_vector_s(j));
                g1_deep_decay(j) = mean(decay_term);
            end
            g1_vector = (N_shallow * 1 + N_deep * g1_deep_decay) / N_total;
        else
            g1_vector(:) = 1.0;
        end
        
        % --- 3手法によるBFI推定 ---
        integrand = (1 - tau_vector_ms(:)/T) .* (g1_vector(:).^2);
        K_simulated = sqrt(beta * (2/T) * trapz(tau_vector_ms, integrand));
        
        % ① SCOS (近似)
        est_BFI_SCOS_Approx = 1 / K_simulated^2;
        
        % ② SCOS (厳密)
        try
            eqn_scos = @(tau_c) calculate_K_from_tau_c(tau_c, T, beta) - K_simulated;
            options = optimset('Display','off');
            estimated_tau_c_ms = fzero(eqn_scos, 0.1, options);
            estimated_gamma_scos = 1 / (estimated_tau_c_ms / 1000);
            mean_transport_path_length_brain = mean(mus_prime_brain * path_brain_deep);
            est_BFI_SCOS_Rigorous = estimated_gamma_scos / (2 * k0^2 * mean_transport_path_length_brain);
        catch
            est_BFI_SCOS_Rigorous = NaN;
        end
        
        % ③ DCS (フィット)
        try
            fit_model_dcs = fittype('a * exp(-b * x) + c');
            offset_guess = N_shallow / N_total;
            gamma_guess = 2 * k0^2 * current_true_BFI * mean(mus_prime_brain * path_brain_deep);
            b_guess = gamma_guess / (1-offset_guess) / 1000;
            fit_options = fitoptions(fit_model_dcs); fit_options.StartPoint = [1-offset_guess, b_guess, offset_guess];
            [fit_obj, ~] = fit(tau_vector_ms(:), g1_vector(:), fit_model_dcs, fit_options);
            fit_gamma = fit_obj.b * (1-fit_obj.c) * 1000;
            mean_transport_path_length_brain = mean(mus_prime_brain * path_brain_deep);
            est_BFI_DCS = fit_gamma / (2 * k0^2 * mean_transport_path_length_brain);
        catch
             est_BFI_DCS = NaN;
        end
        
        % --- 結果の保存 ---
        all_results.(sds_label).True_BFI(bfi_idx) = current_true_BFI;
        all_results.(sds_label).BFI_SCOS_Approx(bfi_idx) = est_BFI_SCOS_Approx;
        all_results.(sds_label).BFI_SCOS_Rigorous(bfi_idx) = est_BFI_SCOS_Rigorous;
        all_results.(sds_label).BFI_DCS(bfi_idx) = est_BFI_DCS;
    end
end

%% =========================================================================
% ステージ3: 相対BFIの精度評価とプロット
% =========================================================================
fprintf('\n--- ステージ3: 相対BFIの精度評価とプロット ---\n');

true_relative_bfi = bfi_list / bfi_list(1); % [1, 2, 3, 4, 5]

% --- SDSごとにFigureを生成 ---
for sds_idx = 1:num_detectors
    sds_label = sprintf('sds_%dmm', sds_mm(sds_idx));
    
    % 対応するSDSの結果テーブルを取得
    sds_results_table = all_results.(sds_label);
    
    % --- 各手法の相対BFIを計算 ---
    % 基準となるBFI (リストの最初の値) で各推定値を割る
    baseline_scos_approx = sds_results_table.BFI_SCOS_Approx(1);
    relative_scos_approx = sds_results_table.BFI_SCOS_Approx / baseline_scos_approx;
    
    baseline_scos_rigorous = sds_results_table.BFI_SCOS_Rigorous(1);
    relative_scos_rigorous = sds_results_table.BFI_SCOS_Rigorous / baseline_scos_rigorous;
    
    baseline_dcs = sds_results_table.BFI_DCS(1);
    relative_dcs = sds_results_table.BFI_DCS / baseline_dcs;
    
    % --- プロット ---
    figure('Name', sprintf('Relative BFI Accuracy at SDS = %d mm', sds_mm(sds_idx)));
    hold on;
    
    % 理想的な線 (y=x)
    plot(true_relative_bfi, true_relative_bfi, 'k--', 'LineWidth', 2, 'DisplayName', '理想値 (y=x)');
    
    % 各手法の結果
    plot(true_relative_bfi, relative_scos_approx, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'SCOS (近似: 1/K^2)');
    plot(true_relative_bfi, relative_scos_rigorous, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'SCOS (厳密: 1/\tau_c)');
    plot(true_relative_bfi, relative_dcs, '^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DCS (フィット)');
    
    hold off;
    
    % グラフの装飾
    title(sprintf('相対BFIの精度評価 (SDS = %d mm)', sds_mm(sds_idx)));
    xlabel('真の相対BFI (基準値からの倍率)');
    ylabel('推定された相対BFI (基準値からの倍率)');
    legend('show', 'Location', 'northwest');
    grid on;
    axis equal; % x軸とy軸のスケールを合わせる
    xlim([0, max(relative_levels) * 1.1]);
    ylim([0, max(relative_levels) * 1.1]);
    xticks(relative_levels);
    yticks(relative_levels);
end

%% =========================================================================
% ヘルパー関数
% =========================================================================
function K = calculate_K_from_tau_c(tau_c, T, beta)
    if tau_c == 0; K = sqrt(beta); return; end
    x = T / tau_c;
    if x == 0; K_sq = beta; else; K_sq = beta * ( (2/x) * (1 - (1/x)*(1-exp(-x))) ); end
    K = sqrt(K_sq);
end