% Corrected group-comparison statistics for wDICE, SNR, and PVE
%
% Applies the same corrected procedure used in ODF_based_repro_comparisons_v2.m:
%   - Dataset comparison (7TS vs 7TL vs 3T, at each resolution): Kruskal-Wallis
%     omnibus + Dunn's test post-hoc (Holm-Bonferroni corrected). Replaces the
%     default Tukey-Kramer post-hoc used in the original wDICE/SNR scripts,
%     which is the wrong pairing for a rank-based omnibus test.
%   - Resolution comparison (1.5 vs 2.0 vs 2.5mm, within each dataset):
%     Friedman omnibus (repeated-measures, since every subject was scanned at
%     all 3 resolutions) + Wilcoxon signed-rank post-hoc (Holm-Bonferroni).
%
% NOTE on PVE: no dataset- or resolution-comparison code for PVE (non-WM
% volume fraction) was found anywhere in the deposited analysis repository
% (github.com/ChenNanHao/DSI-reproducibility-analysis) -- only the
% correlation/regression code (PVE_vs_repro.m, SNR_PVE_regression.m) touches
% PVE. This script's PVE sections are therefore NEW code, not a fix to a
% pre-existing analysis, using the same cached PVE_single/PVE_crossing data
% (G:\dsi_data_7T_20200901_try\PVE\PVE results.mat) and the same statistical
% procedure applied to the other metrics for consistency.

clearvars; clc; close all
warning('off','all');

idx_7TS = [1:3 5:6];
idx_7TL = 7:11;
idx_3T  = [12:13 15:18];
dsnames = {'7TS','7TL','3T'};
resnames = {'1.5mm','2.0mm','2.5mm'};
ds_pairs = [1 2; 1 3; 2 3];
res_pairs = [1 2; 1 3; 2 3];

function out = run_metric(mname, G7TS, G7TL, G3T, dsnames, resnames, ds_pairs, res_pairs)
    % G7TS/G7TL/G3T are [n_subjects x 3 resolutions]
    groups = {G7TS, G7TL, G3T};
    fprintf('\n============================================================\n');
    fprintf('%s\n', mname);
    fprintf('============================================================\n');

    %% Dataset comparison: Kruskal-Wallis + Dunn/Holm
    p_dataset = zeros(1,3);
    posthoc_dataset = cell(1,3);
    for rr = 1:3
        x = cat(1, G7TS(:,rr), G7TL(:,rr), G3T(:,rr));
        group = [ones(size(G7TS,1),1); 2*ones(size(G7TL,1),1); 3*ones(size(G3T,1),1)];
        [p, ~, stats] = kruskalwallis(x', group, 'off');
        p_dataset(rr) = p;

        N = sum(stats.n); R = stats.meanranks; n = stats.n; sumt = stats.sumt;
        p_raw = zeros(3,1);
        for k = 1:3
            i = ds_pairs(k,1); j = ds_pairs(k,2);
            se = sqrt( ((N*(N+1)/12) - (sumt/(12*(N-1)))) * (1/n(i) + 1/n(j)) );
            z = (R(i) - R(j)) / se;
            p_raw(k) = 2*(1 - normcdf(abs(z)));
        end
        [p_sorted, order] = sort(p_raw);
        kk = numel(p_sorted);
        p_adj_sorted = zeros(kk,1);
        p_adj_sorted(1) = min(1, kk * p_sorted(1));
        for jj = 2:kk
            p_adj_sorted(jj) = max(p_adj_sorted(jj-1), min(1, (kk-jj+1) * p_sorted(jj)));
        end
        p_holm = zeros(kk,1);
        p_holm(order) = p_adj_sorted;

        posthoc_dataset{rr} = table(dsnames(ds_pairs(:,1))', dsnames(ds_pairs(:,2))', p_raw, p_holm, ...
            'VariableNames', {'Group_A','Group_B','p_raw_dunn','p_holm'});
    end
    fprintf('\n-- Dataset comparison (7TS vs 7TL vs 3T), Kruskal-Wallis + Dunn/Holm --\n');
    for rr = 1:3
        fprintf('  %s: omnibus p = %.4f\n', resnames{rr}, p_dataset(rr));
        disp(posthoc_dataset{rr})
    end

    %% Resolution comparison: Friedman + Wilcoxon/Holm
    p_resolution = zeros(1,3);
    posthoc_resolution = cell(1,3);
    for dd = 1:3
        CC = groups{dd};
        p_fried = friedman(CC, 1, 'off');
        p_resolution(dd) = p_fried;

        p_raw = zeros(3,1);
        for pp = 1:3
            p_raw(pp) = signrank(CC(:,res_pairs(pp,1)), CC(:,res_pairs(pp,2)));
        end
        [p_sorted, order] = sort(p_raw);
        kk = numel(p_sorted);
        p_adj_sorted = zeros(kk,1);
        p_adj_sorted(1) = min(1, kk * p_sorted(1));
        for jj = 2:kk
            p_adj_sorted(jj) = max(p_adj_sorted(jj-1), min(1, (kk-jj+1) * p_sorted(jj)));
        end
        p_holm = zeros(kk,1);
        p_holm(order) = p_adj_sorted;

        posthoc_resolution{dd} = table(resnames(res_pairs(:,1))', resnames(res_pairs(:,2))', p_raw, p_holm, ...
            'VariableNames', {'Res_A','Res_B','p_raw_signrank','p_holm'});
    end
    fprintf('\n-- Resolution comparison (1.5 vs 2.0 vs 2.5mm), Friedman + Wilcoxon/Holm --\n');
    for dd = 1:3
        fprintf('  %s: Friedman omnibus p = %.4f\n', dsnames{dd}, p_resolution(dd));
        disp(posthoc_resolution{dd})
    end

    out.p_dataset = p_dataset;
    out.posthoc_dataset = posthoc_dataset;
    out.p_resolution = p_resolution;
    out.posthoc_resolution = posthoc_resolution;
end

Results = struct();

%% --- wDICE (whole, homogeneous, heterogeneous) ---
cd('G:\dsi_data_7T_20200901_try')
w = load('wDICE CST results.mat');  % wDICE, wDICE_homo, wDICE_hetero: [3 res x 18 subj]

wnames = {'wDICE_whole','wDICE_homo','wDICE_hetero'};
wdata = {w.wDICE, w.wDICE_homo, w.wDICE_hetero};
for i = 1:3
    M = wdata{i}';  % transpose to [18 subj x 3 res]
    G7TS = M(idx_7TS,:); G7TL = M(idx_7TL,:); G3T = M(idx_3T,:);
    Results.(wnames{i}) = run_metric(wnames{i}, G7TS, G7TL, G3T, dsnames, resnames, ds_pairs, res_pairs);
end

%% --- SNR ---
s = load('SNR_results_20220527.mat');  % SNR_7TS/7TL/3T already [n_subj x 3]
Results.SNR = run_metric('SNR', s.SNR_7TS, s.SNR_7TL, s.SNR_3T, dsnames, resnames, ds_pairs, res_pairs);

%% --- PVE (single-fiber, crossing-fiber) ---
cd('G:\dsi_data_7T_20200901_try\PVE')
pv = load('PVE results.mat');  % PVE_single, PVE_crossing: [18 subj x 3 res]

pvnames = {'PVE_single','PVE_crossing'};
pvdata = {pv.PVE_single, pv.PVE_crossing};
for i = 1:2
    M = pvdata{i};
    G7TS = M(idx_7TS,:); G7TL = M(idx_7TL,:); G3T = M(idx_3T,:);
    Results.(pvnames{i}) = run_metric(pvnames{i}, G7TS, G7TL, G3T, dsnames, resnames, ds_pairs, res_pairs);
end

%% --- PVE by CST segment (whole tract, homogeneous, heterogeneous) ---
% Tract-masked PVE (see Compute_PVE_CST_segments.m), same segmentation as
% wDICE_whole/homo/hetero above.
cd('G:\dsi_data_7T_20200901_try')
pvc = load('PVE_CST_segment_results.mat');  % PVE_whole, PVE_homo, PVE_hetero: [18 subj x 3 res]

pvcnames = {'PVE_whole','PVE_homo','PVE_hetero'};
pvcdata = {pvc.PVE_whole, pvc.PVE_homo, pvc.PVE_hetero};
for i = 1:3
    M = pvcdata{i};
    G7TS = M(idx_7TS,:); G7TL = M(idx_7TL,:); G3T = M(idx_3T,:);
    Results.(pvcnames{i}) = run_metric(pvcnames{i}, G7TS, G7TL, G3T, dsnames, resnames, ds_pairs, res_pairs);
end

cd('G:\dsi_data_7T_20200901_try')
save('SNR_PVE_wDICE_stats_corrected.mat', 'Results');
fprintf('\n\n=== ALL METRICS COMPLETE. Saved to SNR_PVE_wDICE_stats_corrected.mat ===\n');
