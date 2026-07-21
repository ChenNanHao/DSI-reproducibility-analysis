% ODF-based reproducibility statistics — corrected version (v2)
%
% Fixes applied relative to the original ODF_based_repro_comparisons.m:
%   1. Dataset-comparison post-hoc: Tukey-Kramer -> Dunn's test (manual),
%      Holm-Bonferroni corrected.
%      (Tukey-Kramer is a means-based ANOVA procedure, the wrong pairing
%      for a rank-based Kruskal-Wallis omnibus test. MATLAB's multcompare
%      does provide a 'dunn-sidak' CriticalValueType option, but it
%      returns NaN p-values / infinite CIs specifically when the input
%      stats struct has source='kruskalwallis' -- confirmed by testing
%      the identical call against anova1 stats, where it works correctly.
%      This appears to be a MATLAB Statistics Toolbox limitation, not a
%      usage error. Dunn's test (Dunn, 1964) is implemented directly here
%      using the mean-ranks/n/sumt fields kruskalwallis already returns,
%      which sidesteps multcompare entirely for this comparison.)
%   2. Resolution comparison: previously hardcoded CC(1,:),CC(2,:),CC(3,:),
%      silently dropping subjects beyond the third regardless of true
%      group size. Replaced with Friedman (proper repeated-measures
%      omnibus test, since every subject was scanned at all 3
%      resolutions) + Wilcoxon signed-rank pairwise post-hoc with
%      Holm-Bonferroni correction, using ALL subjects.
%   3. Every metric/fiber-group combination now gets uniquely named
%      output variables (no more p_SIMODF / c_SIMODF being silently
%      overwritten when the script moves from single-fiber to
%      crossing-fiber), and everything is saved to disk.

clearvars; clc; close all
warning('off','all');

cd('G:\dsi_data_7T_20200901_try\2p5mm_results_20220411')
load('ODF results.mat', 'Single_SIMODF', 'Crossing_SIMODF', 'Single_Angledev', 'Crossing_Angledev');

idx_7TS = [1:3 5:6];
idx_7TL = 7:11;
idx_3T  = [12:13 15:18];
dsnames = {'7TS','7TL','3T'};
resnames = {'1.5mm','2.0mm','2.5mm'};

metrics = struct( ...
    'name', {'SIMODF_single','SIMODF_crossing','AngleDev_single','AngleDev_crossing'}, ...
    'data', {Single_SIMODF, Crossing_SIMODF, Single_Angledev, Crossing_Angledev});

Results = struct();

for m = 1:numel(metrics)
    mname = metrics(m).name;
    M = metrics(m).data;   % n_subjects_total x 3 resolutions

    G7TS = M(idx_7TS,:);
    G7TL = M(idx_7TL,:);
    G3T  = M(idx_3T,:);
    groups = {G7TS, G7TL, G3T};

    fprintf('\n============================================================\n');
    fprintf('%s\n', mname);
    fprintf('============================================================\n');

    %% --- (A) Dataset comparison: across 7TS/7TL/3T, at each resolution ---
    % Correct as originally coded (subject-level means, all subjects) --
    % only the post-hoc method changes: Tukey-Kramer -> Dunn's test (manual),
    % Holm-Bonferroni corrected.
    p_dataset = zeros(1,3);
    posthoc_dataset = cell(1,3);
    ds_pairs = [1 2; 1 3; 2 3];
    for rr = 1:3
        x = cat(1, G7TS(:,rr), G7TL(:,rr), G3T(:,rr));
        group = [ones(size(G7TS,1),1); 2*ones(size(G7TL,1),1); 3*ones(size(G3T,1),1)];
        [p, ~, stats] = kruskalwallis(x', group, 'off');
        p_dataset(rr) = p;

        N = sum(stats.n);
        R = stats.meanranks;
        n = stats.n;
        sumt = stats.sumt;
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

    %% --- (B) Resolution comparison: across 1.5/2.0/2.5mm, within each dataset ---
    % Corrected: Friedman omnibus (repeated measures) + Wilcoxon signed-rank
    % post-hoc with Holm-Bonferroni correction, using ALL subjects.
    p_resolution = zeros(1,3);
    posthoc_resolution = cell(1,3);
    pairs = [1 2; 1 3; 2 3];
    for dd = 1:3
        CC = groups{dd};   % n_subjects x 3 resolutions
        p_fried = friedman(CC, 1, 'off');
        p_resolution(dd) = p_fried;

        p_raw = zeros(3,1);
        for pp = 1:3
            p_raw(pp) = signrank(CC(:,pairs(pp,1)), CC(:,pairs(pp,2)));
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

        posthoc_resolution{dd} = table(resnames(pairs(:,1))', resnames(pairs(:,2))', p_raw, p_holm, ...
            'VariableNames', {'Res_A','Res_B','p_raw_signrank','p_holm'});
    end
    fprintf('\n-- Resolution comparison (1.5 vs 2.0 vs 2.5mm), Friedman + Wilcoxon/Holm --\n');
    for dd = 1:3
        fprintf('  %s: Friedman omnibus p = %.4f\n', dsnames{dd}, p_resolution(dd));
        disp(posthoc_resolution{dd})
    end

    % store uniquely per metric — no overwriting between single/crossing
    Results.(mname).p_dataset = p_dataset;
    Results.(mname).posthoc_dataset = posthoc_dataset;
    Results.(mname).p_resolution = p_resolution;
    Results.(mname).posthoc_resolution = posthoc_resolution;
end

cd('G:\dsi_data_7T_20200901_try\2p5mm_results_20220411')
save('ODF_stats_corrected_v2.mat', 'Results');
fprintf('\n\n=== ALL METRICS COMPLETE. Saved to ODF_stats_corrected_v2.mat ===\n');
