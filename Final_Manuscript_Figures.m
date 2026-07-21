%% Final_Manuscript_Figures.m
%
% Regenerates Figures 3-7 of the DSI reproducibility manuscript using the
% CORRECTED statistics (Dunn's test + Holm-Bonferroni for dataset
% comparisons, Friedman test for resolution comparisons, linear
% mixed-effects models for SNR/PVE associations), matching the text now
% in Results paragraphs 53-58 of the manuscript.
%
% NOT covered here:
%   - Figure 2 (voxel-wise SIMODF/AngleDev maps) - this is a qualitative
%     map/slice visualization, not a statistical summary figure, and
%     does not depend on any of the corrected group-comparison stats.
%
% Figure 5 is now SNR-only (single panel). The PVE panels formerly (b)/(c)
% were removed along with all PVE reporting from the manuscript; the
% tract-masked PVE pipeline (Compute_PVE_CST_segments.m,
% PVE_CST_segment_results.mat) is no longer used by this script.
%
% Style (bar color, marker shape, resolution color, fonts) is matched to
% the original plotting scripts in pipeline\ so the regenerated figures
% look consistent with Figures 1 and 2 (which are unaffected by this
% correction):
%   pipeline\ODF_based_repro_comparisons.m          (Figure 3 style)
%   pipeline\Tracts_shape_similarity_DICE_CST_20220412.m (Figure 4 style)
%   pipeline\SNR calculation 20220510.m, SNR_PVE_regression.m (Figure 5 style)
%   pipeline\PVE_vs_repro.m                         (Figure 6-7 style)
%
% Output: TIFF files at 600 dpi (exceeds NeuroImage/Elsevier's 500 dpi
% minimum for combination line+color artwork) saved to
% D:\DSI_manuscript\figures\Figure3.tif ... Figure7.tif
% PNG copies (also Elsevier-acceptable) are saved alongside for quick
% on-screen review.

clearvars; clc; close all
warning('off','all');
cd('G:\dsi_data_7T_20200901_try')

outdir = 'D:\DSI_manuscript\figures';
if ~exist(outdir,'dir'), mkdir(outdir); end

idx_7TS = [1:3 5:6]; idx_7TL = 7:11; idx_3T = [12:13 15:18];
dsnames  = {'7TS','7TL','3T'};
resnames = {'1.5 mm','2.0 mm','2.5 mm'};
bar_face  = [0 0 0; 0.7 0.7 0.7; 0.5 0.5 0.5];      % 7TS / 7TL / 3T
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; % 1.5/2.0/2.5 mm
dataset_marker = {'o','*','s'};                      % 7TS / 7TL / 3T

%% ===================== load all data =====================
odf  = load('2p5mm_results_20220411\ODF results.mat');
w    = load('wDICE CST results.mat');
snr  = load('SNR_results_20220527.mat');
cd('PVE'); pv = load('PVE results.mat'); cd('G:\dsi_data_7T_20200901_try')

SNR_all = nan(18,3);
SNR_all(idx_7TS,:) = snr.SNR_7TS; SNR_all(idx_7TL,:) = snr.SNR_7TL; SNR_all(idx_3T,:) = snr.SNR_3T;

metric_data = struct( ...
    'SIMODF_single',    odf.Single_SIMODF, ...
    'SIMODF_crossing',  odf.Crossing_SIMODF, ...
    'AngleDev_single',  odf.Single_Angledev, ...
    'AngleDev_crossing',odf.Crossing_Angledev, ...
    'wDICE_whole',      w.wDICE', ...
    'wDICE_homo',       w.wDICE_homo', ...
    'wDICE_hetero',     w.wDICE_hetero', ...
    'SNR',              SNR_all, ...
    'PVE_single',       pv.PVE_single, ...
    'PVE_crossing',     pv.PVE_crossing);

odf_stats   = load('2p5mm_results_20220411\ODF_stats_corrected_v2.mat').Results;
wps_stats   = load('SNR_PVE_wDICE_stats_corrected.mat').Results;
stats = struct(); fn1 = fieldnames(odf_stats); for k=1:numel(fn1), stats.(fn1{k}) = odf_stats.(fn1{k}); end
fn2 = fieldnames(wps_stats); for k=1:numel(fn2), stats.(fn2{k}) = wps_stats.(fn2{k}); end

a2 = load('Aim2_analysis_results.mat');

%% ===================== FIGURE 3: ODF-based group comparisons =====================
% Panels: SIMODF single-fiber, SIMODF crossing-fiber,
%         AngleDev single-fiber, AngleDev crossing-fiber
fig = figure('Position',[100 100 1450 1050],'Color','w');
% Panel height trimmed to reserve a header band above/between rows for the
% "(a)/(b)" row label + "Single-fiber"/"Crossing-fiber" column headers.
% Bands widened further to accommodate the larger header/label font sizes.
avg_pos3 = [ ...
    0.0552 0.53  0.4099 0.35; ...   % top-left    (SIMODF single)
    0.5350 0.53  0.4099 0.35; ...   % top-right   (SIMODF crossing)
    0.0552 0.06  0.4099 0.35; ...   % bottom-left (AngleDev single)
    0.5350 0.06  0.4099 0.35];      % bottom-right(AngleDev crossing)
panels3 = {'SIMODF_single','SIMODF_crossing','AngleDev_single','AngleDev_crossing'};
ylabels3 = {' SIM_O_D_F ',' SIM_O_D_F ',' Angle_D_e_v ',' Angle_D_e_v '};
ylims3  = {[0 3],[0 3],[0 20],[0 60]};
yticks3 = {0:1:3, 0:1:3, 0:5:20, 0:20:60};
for pIdx = 1:4
    ax = axes(fig,'Position',avg_pos3(pIdx,:));
    name = panels3{pIdx};
    M = metric_data.(name);
    means = [mean(M(idx_7TS,:))' mean(M(idx_7TL,:))' mean(M(idx_3T,:))'];
    stds  = [std(M(idx_7TS,:))'  std(M(idx_7TL,:))'  std(M(idx_3T,:))'];
    b = plot_group_bars(ax, means, stds, bar_face, resnames, ylabels3{pIdx}, ylims3{pIdx}, yticks3{pIdx}, 18, 18);
    add_dataset_sig_brackets(ax, means, stds, stats.(name).posthoc_dataset, ylims3{pIdx}, 18);
    if pIdx == 1
        legend(b, {' 7T_S_T_E',' 7T_L_T_E',' 3T'}, 'FontName','calibri','FontSize',18,'NumColumns',3,'Orientation','horizontal','Location','northwest','TextColor','k'); legend boxoff
    end
end
% Row labels "(a)"/"(b)" (bold italic, 28pt) and column headers "Single-fiber"/
% "Crossing-fiber" (italic, not bold, 24pt) - matches the style of the
% original manuscript Figure 3.
axoverlay = axes(fig,'Position',[0 0 1 1],'Visible','off','XLim',[0 1],'YLim',[0 1]); hold(axoverlay,'on');
colx = [0.0552+0.4099/2, 0.5350+0.4099/2];
rowy_header = [0.92, 0.45];
rowy_label  = rowy_header;
rowlabels3 = {'(a)','(b)'};
coltitles3 = {'Single-fiber','Crossing-fiber'};
for r = 1:2
    text(axoverlay, 0.015, rowy_label(r), rowlabels3{r}, 'FontName','calibri','FontSize',28,'FontWeight','bold','FontAngle','italic','Color','k','HorizontalAlignment','left','VerticalAlignment','middle');
    for c = 1:2
        text(axoverlay, colx(c), rowy_header(r), coltitles3{c}, 'FontName','calibri','FontSize',24,'FontWeight','normal','Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
save_figure(fig, fullfile(outdir,'Figure3'));

%% ===================== FIGURE 4: tract-based (wDICE) group comparisons =====================
fig = figure('Position',[100 100 1500 750],'Color','w');
% Layout/typography matched to Figure 3: no figure title, italic (non-bold,
% 24pt) subtitles above each panel, dataset-comparison asterisks only,
% shared legend moved inside the first panel.
avg_pos4 = [ ...
    0.0500 0.10 0.2733 0.77; ...   % Whole tract
    0.3733 0.10 0.2733 0.77; ...   % Homogeneous region
    0.6967 0.10 0.2733 0.77];      % Heterogeneous region
panels4 = {'wDICE_whole','wDICE_homo','wDICE_hetero'};
titles4 = {'Whole tract','Homogeneous region','Heterogeneous region'};
for pIdx = 1:3
    ax = axes(fig,'Position',avg_pos4(pIdx,:));
    name = panels4{pIdx};
    M = metric_data.(name);
    means = [mean(M(idx_7TS,:))' mean(M(idx_7TL,:))' mean(M(idx_3T,:))'];
    stds  = [std(M(idx_7TS,:))'  std(M(idx_7TL,:))'  std(M(idx_3T,:))'];
    b = plot_group_bars(ax, means, stds, bar_face, resnames, ' wDICE coefficient (%)', [0 120], 0:40:120, 18, 18);
    hold(ax,'on'); plot(ax, linspace(0.5,3.5,100), ones(1,100)*72, 'r-', 'linewidth',2);
    add_dataset_sig_brackets(ax, means, stds, stats.(name).posthoc_dataset, [0 120], 18);
    if pIdx == 1
        legend(b, {' 7T_S_T_E',' 7T_L_T_E',' 3T'}, 'FontName','calibri','FontSize',18,'NumColumns',3,'Orientation','horizontal','Location','northwest','TextColor','k'); legend boxoff
    end
end
axoverlay4 = axes(fig,'Position',[0 0 1 1],'Visible','off','XLim',[0 1],'YLim',[0 1]); hold(axoverlay4,'on');
colx4 = avg_pos4(:,1) + avg_pos4(:,3)/2;
for pIdx = 1:3
    text(axoverlay4, colx4(pIdx), 0.91, titles4{pIdx}, 'FontName','calibri','FontSize',24,'FontWeight','normal','Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
end
save_figure(fig, fullfile(outdir,'Figure4'));

%% ===================== FIGURE 5: SNR group comparisons =====================
% Single panel: SNR only. PVE panels (formerly (b)/(c)) have been removed
% from the manuscript entirely (SNR-only Objective 2, no PVE reporting) -
% see Compute_PVE_CST_segments.m / PVE_CST_segment_results.mat for the
% now-unused tract-masked PVE pipeline this panel used to draw on.
fig = figure('Position',[100 100 650 600],'Color','w');
ax = axes(fig,'Position',[0.16 0.14 0.80 0.72]);
M = metric_data.SNR;
means = [mean(M(idx_7TS,:))' mean(M(idx_7TL,:))' mean(M(idx_3T,:))'];
stds  = [std(M(idx_7TS,:))'  std(M(idx_7TL,:))'  std(M(idx_3T,:))'];
b = plot_group_bars(ax, means, stds, bar_face, resnames, ' SNR (a.u.) ', [0 45], 0:15:45, 18, 18);
add_dataset_sig_brackets(ax, means, stds, stats.SNR.posthoc_dataset, [0 45], 18);
legend(b, {' 7T_S_T_E',' 7T_L_T_E',' 3T'}, 'FontName','calibri','FontSize',18,'NumColumns',3,'Orientation','horizontal','Location','northwest','TextColor','k'); legend boxoff
save_figure(fig, fullfile(outdir,'Figure5'));

%% ===================== FIGURE 6: ODF-based reproducibility vs SNR/PVE (LMM) =====================
% Rows: SIMODF-vs-SNR, SIMODF-vs-PVE, AngleDev-vs-SNR, AngleDev-vs-PVE
% Cols: single-fiber (left), crossing-fiber (right) - matches original caption convention
fig = figure('Position',[50 50 1100 1800],'Color','w');
tlo6 = tiledlayout(fig,4,2,'TileSpacing','loose','Padding','compact');
plotspecs6 = { ...
    struct('metric','SIMODF_single',  'pred','SNR', 'ylab',' SIM_O_D_F ','ylims',[0 3],  'row',1,'col',1), ...
    struct('metric','SIMODF_crossing','pred','SNR', 'ylab',' SIM_O_D_F ','ylims',[0 3],  'row',1,'col',2), ...
    struct('metric','SIMODF_single',  'pred','PVE', 'ylab',' SIM_O_D_F ','ylims',[0 3],  'row',2,'col',1), ...
    struct('metric','SIMODF_crossing','pred','PVE', 'ylab',' SIM_O_D_F ','ylims',[0 3],  'row',2,'col',2), ...
    struct('metric','AngleDev_single',  'pred','SNR', 'ylab',' Angle_D_e_v ','ylims',[0 20],'row',3,'col',1), ...
    struct('metric','AngleDev_crossing','pred','SNR', 'ylab',' Angle_D_e_v ','ylims',[15 40],'row',3,'col',2), ...
    struct('metric','AngleDev_single',  'pred','PVE', 'ylab',' Angle_D_e_v ','ylims',[0 20],'row',4,'col',1), ...
    struct('metric','AngleDev_crossing','pred','PVE', 'ylab',' Angle_D_e_v ','ylims',[15 40],'row',4,'col',2) ...
};
for k = 1:numel(plotspecs6)
    s = plotspecs6{k};
    ax = nexttile(tlo6, (s.row-1)*2 + s.col);
    Y = metric_data.(s.metric);
    if strcmp(s.pred,'SNR')
        X = SNR_all; xlab = ' SNR (a.u.) '; xlims = [0 40];
    else
        % PVE uses fiber-type-matched predictor: PVE_single for *_single metrics, PVE_crossing for *_crossing metrics
        if endsWith(s.metric,'single'), X = metric_data.PVE_single; else, X = metric_data.PVE_crossing; end
        xlab = ' Non-WM volume fraction (%) '; xlims = [0 40];
    end
    mdl = a2.lme_results.(s.metric).model;
    [beta,pv] = get_lmm_coef(mdl, s.pred);
    r2m = a2.lme_results.(s.metric).R2_marginal;
    r2c = a2.lme_results.(s.metric).R2_conditional;
    plot_lmm_scatter(ax, X, Y, idx_7TS, idx_7TL, idx_3T, resol_color, dataset_marker, xlab, s.ylab, xlims, s.ylims, beta, pv, r2m, r2c);
    title(ax, sprintf('%s vs %s', strrep(s.metric,'_',' '), s.pred), 'FontName','calibri','FontSize',11,'FontWeight','normal');
end
title(tlo6, 'Figure 6. Associations between ODF-based reproducibility and SNR/PVE (linear mixed-effects model)', 'FontName','calibri','FontSize',13);
add_scatter_legend(fig, resol_color, dataset_marker, dsnames, resnames, 0.965);
save_figure(fig, fullfile(outdir,'Figure6'));

%% ===================== FIGURE 7: tract-based reproducibility vs SNR/PVE (LMM) =====================
% Rows: vs SNR, vs PVE (using PVE_single as the tract-relevant PVE proxy)
% Cols: whole tract, homogeneous, heterogeneous
fig = figure('Position',[50 50 1500 1000],'Color','w');
tlo7 = tiledlayout(fig,2,3,'TileSpacing','loose','Padding','compact');
tract_metrics = {'wDICE_whole','wDICE_homo','wDICE_hetero'};
tract_titles  = {'Whole tract','Homogeneous region','Heterogeneous region'};
for col = 1:3
    name = tract_metrics{col};
    Y = metric_data.(name);
    mdl = a2.lme_results.(name).model;
    r2m = a2.lme_results.(name).R2_marginal;
    r2c = a2.lme_results.(name).R2_conditional;

    ax = nexttile(tlo7, col);
    [beta,pv] = get_lmm_coef(mdl,'SNR');
    plot_lmm_scatter(ax, SNR_all, Y, idx_7TS, idx_7TL, idx_3T, resol_color, dataset_marker, ' SNR (a.u.) ', ' wDICE (%) ', [0 40], [40 100], beta, pv, r2m, r2c);
    title(ax, [tract_titles{col} ' vs SNR'], 'FontName','calibri','FontSize',12,'FontWeight','normal');

    ax = nexttile(tlo7, col+3);
    [beta,pv] = get_lmm_coef(mdl,'PVE');
    plot_lmm_scatter(ax, metric_data.PVE_single, Y, idx_7TS, idx_7TL, idx_3T, resol_color, dataset_marker, ' Non-WM volume fraction (%) ', ' wDICE (%) ', [0 30], [40 100], beta, pv, r2m, r2c);
    title(ax, [tract_titles{col} ' vs PVE'], 'FontName','calibri','FontSize',12,'FontWeight','normal');
end
title(tlo7, 'Figure 7. Associations between tract-based reproducibility and SNR/PVE (linear mixed-effects model)', 'FontName','calibri','FontSize',13);
add_scatter_legend(fig, resol_color, dataset_marker, dsnames, resnames, 0.955);
save_figure(fig, fullfile(outdir,'Figure7'));

fprintf('\nAll figures saved to %s\n', outdir);


%% ===================== local functions =====================
function b = plot_group_bars(ax, means, stds, bar_face, resnames, ylab, ylims_, yticks_, fsize, ysize)
    if nargin < 9, fsize = 13; end
    if nargin < 10, ysize = 14; end
    axes(ax); %#ok<LAXES>
    X = categorical(resnames); X = reordercats(X, resnames);
    ngroups = size(means,1); nbars = size(means,2);
    groupwidth = min(0.8, nbars/(nbars+1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
        errorbar(x, means(:,i), stds(:,i), 'k.','linewidth',3); hold on
    end
    b = bar(1:ngroups, means); hold on
    for i = 1:nbars
        set(b(i),'edgecolor',bar_face(i,:),'facecolor',bar_face(i,:));
    end
    xticks(1:ngroups); ylim(ylims_); yticks(yticks_);
    set(gca,'fontname','calibri','fontsize',fsize,'xticklabel',X,'linewidth',1.5,'XColor','k','YColor','k'); box off
    ylabel(ylab,'fontname','calibri','fontsize',ysize,'color','k');
end

function add_dataset_sig_brackets(ax, means, stds, posthoc_dataset, ylims_, fsize)
    % posthoc_dataset{r}: table with Group_A, Group_B, p_holm, rows in order (7TS-7TL, 7TS-3T, 7TL-3T)
    if nargin < 6, fsize = 13; end
    axes(ax); %#ok<LAXES>
    ngroups = size(means,1); nbars = 3;
    groupwidth = min(0.8, nbars/(nbars+1.5));
    yrange = ylims_(2) - ylims_(1);
    for r = 1:ngroups
        xpos = (r) - groupwidth/2 + (2*(1:nbars)-1)*groupwidth/(2*nbars);
        tops = means(r,:) + stds(r,:);
        t = posthoc_dataset{r};
        pairs = [1 2; 1 3; 2 3]; % 7TS-7TL, 7TS-3T, 7TL-3T
        stacklevel = 0;
        for row = 1:height(t)
            p = t.p_holm(row);
            if p >= 0.05, continue; end
            mark = '*'; if p < 0.01, mark = '**'; end
            iA = pairs(row,1); iB = pairs(row,2);
            y = max(tops(iA), tops(iB)) + yrange*0.04 + stacklevel*yrange*0.09;
            line([xpos(iA) xpos(iA) xpos(iB) xpos(iB)], [y-yrange*0.015 y y y-yrange*0.015], 'color','k','linewidth',1.2);
            text(mean([xpos(iA) xpos(iB)]), y+yrange*0.01, mark, 'HorizontalAlignment','center','FontName','calibri','FontSize',fsize,'Color','k');
            stacklevel = stacklevel + 1;
        end
    end
end

function [beta, pval] = get_lmm_coef(mdl, predname)
    names = mdl.Coefficients.Name;
    row = find(strcmp(names, predname));
    beta = mdl.Coefficients.Estimate(row);
    pval = mdl.Coefficients.pValue(row);
end

function plot_lmm_scatter(ax, X, Y, idx_7TS, idx_7TL, idx_3T, resol_color, dataset_marker, xlab, ylab, xlims, ylims_, beta, pv, r2m, r2c)
    axes(ax); %#ok<LAXES>
    idxs = {idx_7TS, idx_7TL, idx_3T};
    for kk = 1:3 % resolution
        for dd = 1:3 % dataset
            plot(X(idxs{dd},kk), Y(idxs{dd},kk), dataset_marker{dd}, 'linewidth',1.5, 'color', resol_color(kk,:), 'MarkerSize',6);
            hold on
        end
    end
    xlim(xlims); ylim(ylims_);
    set(gca,'fontname','calibri','fontsize',11,'linewidth',1.2); box off
    xlabel(xlab,'fontname','calibri','fontsize',12);
    ylabel(ylab,'fontname','calibri','fontsize',12);
    if pv < 0.001, pstr = 'p<0.001'; else, pstr = sprintf('p=%.3f', pv); end
    txt = sprintf('\\beta=%.3f, %s\nR^2_m=%.2f, R^2_c=%.2f', beta, pstr, r2m, r2c);
    xr = xlims(2)-xlims(1); yr = ylims_(2)-ylims_(1);
    text(xlims(1)+0.03*xr, ylims_(2)-0.14*yr, txt, 'FontName','calibri','FontSize',9, 'VerticalAlignment','top');
end

function add_scatter_legend(fig, resol_color, dataset_marker, dsnames, resnames, ypos)
    ax = axes(fig, 'Position',[0.4 ypos 0.2 0.03], 'Visible','off'); hold(ax,'on');
    for kk = 1:3
        plot(ax, kk, 1, 's', 'MarkerFaceColor', resol_color(kk,:), 'MarkerEdgeColor', resol_color(kk,:), 'MarkerSize',8, 'HandleVisibility','on');
    end
    for dd = 1:3
        plot(ax, 3+dd, 1, dataset_marker{dd}, 'color','k','MarkerSize',7, 'HandleVisibility','on');
    end
    lg = legend(ax, [resnames, dsnames], 'Orientation','horizontal', 'FontName','calibri','FontSize',9, 'Location','south');
    legend(ax,'boxoff'); lg.AutoUpdate = 'off';
    xlim(ax,[0 7]); ylim(ax,[0.5 1.5]);
end

function save_figure(fig, basepath)
    exportgraphics(fig, [basepath '.tif'], 'Resolution', 600);
    exportgraphics(fig, [basepath '.png'], 'Resolution', 600);
end
