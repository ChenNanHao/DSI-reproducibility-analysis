% SNR-only version of the Objective 2 analysis, for comparison against the
% SNR+PVE model. Motivation: PVE's unique contribution was consistently
% near zero once SNR was accounted for (Aim2_R2_decomposition.m,
% Resolution_stratified_SNR_PVE_analysis.m), and PVE shows a counter-
% intuitive pattern (elevated in 7TSTE vs 3T despite matched nominal
% resolution - see Figure 5b) that is hard to interpret cleanly. This script
% checks whether a simpler SNR-only model (a) explains nearly as much
% variance as the SNR+PVE model, and (b) gives a cleaner/more consistent
% story, both pooled and resolution-stratified.
%
% Part A: pooled model, Y ~ SNR (+Dataset) + (1|Subject)
% Part B: resolution-stratified, Y ~ SNR (+Dataset), n=16 per resolution

clearvars; clc; close all
warning('off','all');

idx_7TS = [1:3 5:6];
idx_7TL = 7:11;
idx_3T  = [12:13 15:18];
resnames = {'1.5mm','2.0mm','2.5mm'};
res_mm = [1.5 2.0 2.5];

cd('G:\dsi_data_7T_20200901_try')
odf = load('2p5mm_results_20220411\ODF results.mat');
snr = load('SNR_results_20220527.mat');
cd('G:\dsi_data_7T_20200901_try')

metrics = struct( ...
    'name', {'SIMODF_single','SIMODF_crossing','AngleDev_single','AngleDev_crossing'}, ...
    'data', {odf.Single_SIMODF, odf.Crossing_SIMODF, odf.Single_Angledev, odf.Crossing_Angledev});

SNR_all = nan(18,3);
SNR_all(idx_7TS,:) = snr.SNR_7TS;
SNR_all(idx_7TL,:) = snr.SNR_7TL;
SNR_all(idx_3T,:)  = snr.SNR_3T;

keep_idx = [idx_7TS idx_7TL idx_3T];
DatasetLabel = cell(18,1);
for s = idx_7TS, DatasetLabel{s} = '7TS'; end
for s = idx_7TL, DatasetLabel{s} = '7TL'; end
for s = idx_3T,  DatasetLabel{s} = '3T';  end

function R2 = marginal_R2(lme, T)
    yhat_fixed = predict(lme, T, 'Conditional', false);
    var_fixed = var(yhat_fixed, 'omitnan');
    [psi, mse] = covarianceParameters(lme);
    var_random = psi{1}(1,1);
    var_resid = mse;
    R2 = var_fixed / (var_fixed + var_random + var_resid);
end

%% ================= Part A: pooled SNR-only model =================
fprintf('\n########## PART A: POOLED, SNR-ONLY (N=48) ##########\n');
snr_only_pooled = struct();
for m = 1:numel(metrics)
    mname = metrics(m).name;
    Mdata = metrics(m).data;

    Subject = []; Dataset = {}; SNRv = []; Y = [];
    for s = 1:18
        if ismember(s, idx_7TS), dsl = '7TS';
        elseif ismember(s, idx_7TL), dsl = '7TL';
        elseif ismember(s, idx_3T), dsl = '3T';
        else, continue;
        end
        for r = 1:3
            Subject(end+1,1) = s; %#ok<*AGROW>
            Dataset{end+1,1} = dsl;
            SNRv(end+1,1) = SNR_all(s,r);
            Y(end+1,1) = Mdata(s,r);
        end
    end
    T = table(categorical(Subject), categorical(Dataset), SNRv, Y, ...
        'VariableNames', {'Subject','Dataset','SNR','Y'});

    lme_snr  = fitlme(T, 'Y ~ SNR + (1|Subject)');
    lme_full = fitlme(T, 'Y ~ SNR + Dataset + (1|Subject)');

    R2_snr_marg = marginal_R2(lme_snr, T);
    R2_full_marg = marginal_R2(lme_full, T);
    yhat_fixed = predict(lme_full, T, 'Conditional', false);
    var_fixed = var(yhat_fixed, 'omitnan');
    [psi, mse] = covarianceParameters(lme_full);
    var_random = psi{1}(1,1); var_resid = mse;
    R2_full_cond = (var_fixed + var_random) / (var_fixed + var_random + var_resid);

    cmp = compare(lme_snr, lme_full);
    p_dataset = cmp.pValue(2);

    beta_snr = lme_snr.Coefficients.Estimate(strcmp(lme_snr.Coefficients.Name,'SNR'));
    p_snr = lme_snr.Coefficients.pValue(strcmp(lme_snr.Coefficients.Name,'SNR'));

    fprintf('\n-- %s --\n', mname);
    fprintf('  SNR-only model:  beta=%.4f, p=%.4g, R2_marginal=%.4f\n', beta_snr, p_snr, R2_snr_marg);
    fprintf('  +Dataset model:  R2_marginal=%.4f, R2_conditional=%.4f\n', R2_full_marg, R2_full_cond);
    fprintf('  Dataset LRT p (beyond SNR alone) = %.4f\n', p_dataset);

    snr_only_pooled.(mname) = struct('beta_snr',beta_snr,'p_snr',p_snr, ...
        'R2_marginal_snronly',R2_snr_marg,'R2_marginal_full',R2_full_marg, ...
        'R2_conditional_full',R2_full_cond,'p_dataset',p_dataset);
end

%% ================= Part B: resolution-stratified SNR-only =================
fprintf('\n\n########## PART B: RESOLUTION-STRATIFIED, SNR-ONLY (n=16 per resolution) ##########\n');
snr_only_strat = struct();
for m = 1:numel(metrics)
    mname = metrics(m).name;
    Mdata = metrics(m).data;
    fprintf('\n============================================================\n%s\n============================================================\n', mname);

    for r = 1:3
        Y = Mdata(keep_idx, r);
        SNRv = SNR_all(keep_idx, r);
        Dataset = categorical(DatasetLabel(keep_idx));
        T = table(Dataset, SNRv, Y, 'VariableNames', {'Dataset','SNR','Y'});

        mdl_snr  = fitlm(T, 'Y ~ SNR');
        mdl_full = fitlm(T, 'Y ~ SNR + Dataset');

        R2_snr = mdl_snr.Rsquared.Ordinary;
        beta_snr = mdl_snr.Coefficients.Estimate(strcmp(mdl_snr.Coefficients.Properties.RowNames,'SNR'));
        p_snr = mdl_snr.Coefficients.pValue(strcmp(mdl_snr.Coefficients.Properties.RowNames,'SNR'));

        anova_tbl = anova(mdl_full);
        p_dataset = anova_tbl.pValue(strcmp(anova_tbl.Properties.RowNames,'Dataset'));

        fprintf('  -- %s (n=16) --  beta=%.4f, p=%.4f, R2=%.4f | Dataset (beyond SNR) p=%.4f\n', ...
            resnames{r}, beta_snr, p_snr, R2_snr, p_dataset);

        snr_only_strat.(mname).(sprintf('res%d',r)) = struct('beta_snr',beta_snr,'p_snr',p_snr,'R2',R2_snr,'p_dataset',p_dataset);
    end
end

cd('G:\dsi_data_7T_20200901_try')
save('SNR_only_results.mat', 'snr_only_pooled', 'snr_only_strat');
fprintf('\n=== DONE. Saved SNR_only_results.mat ===\n');
