% SNR calculation — corrected version (v2)
%
% Fixes applied relative to the original "SNR calculation 20220510.m":
%   1. Path fix: raw data has moved from E:\dsi_data_7T_20200901_try\ to
%      G:\dsi_data_7T_20200901_try\ since the original script was written.
%      The original script referenced E:\ throughout and was never actually
%      run to completion after the drive letter changed -- confirmed by the
%      complete absence of "SNR_results_20220527.mat" anywhere on disk.
%   2. Statistics: this script now only computes SNR itself (unchanged raw
%      computation logic). Dataset/resolution group comparisons are handled
%      by SNR_PVE_wDICE_dataset_resolution_stats.m, which applies the same
%      corrected Dunn's test / Friedman procedure used for the other metrics.

clearvars; clc; close all
warning('off','all');

cd('G:\dsi_data_7T_20200901_try\SRC_Batch')
noise_15 = single(niftiread('rcT1_pre_brain_roi_noise_1.5mm.nii'));
noise_20 = single(niftiread('rcT1_pre_brain_roi_noise_2.0mm.nii'));
noise_25 = single(niftiread('rcT1_pre_brain_roi_noise_2.5mm.nii'));

noise_15 = noise_15 + fliplr(noise_15);
noise_20 = noise_20 + fliplr(noise_20);
noise_25 = noise_25 + fliplr(noise_25);
noise_15(noise_15 == 0) = nan;
noise_20(noise_20 == 0) = nan;
noise_25(noise_25 == 0) = nan;
    % 7T Longer TE : HY LWK MS OK SY [4 8 10 11 15]
    % 7T Shorter TE : IS MA HM HR OS KW [2 3 5 7 9 12]
    % 3T : CH KHC Peter RL YP YRC [1 6 13 14 16 17]

%% 1p5mm SNR
cd('G:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*1p5mm');
for ii = [1:2:size(folder,1)-1]
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(niftiread('c2rcT1_pre.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    elseif ii == 9
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1;
    else
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    end
    SNR_15(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_15p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    fprintf('Finished SNR computation at %s\n', folder(ii).name);
end

%% 2p0mm SNR
cd('G:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*2p0mm');
for ii = [1:2:size(folder,1)-1]
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(niftiread('c2rcT1_pre.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    elseif ii == 9
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1;
    elseif ii == 23
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1;
    else
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    end
    SNR_20(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_20p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    fprintf('Finished SNR computation at %s\n', folder(ii).name);
end

%% 2p5mm SNR
cd('G:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*2p5mm');
for ii = [1:2:size(folder,1)-1]
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('G:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    elseif ii == 9
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1;
    elseif ii == 23
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1;
    else
        WMmask = double(niftiread('c2rcT1_pre_brain.nii'));
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1;
    end
    SNR_25(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_25p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    fprintf('Finished SNR computation at %s\n', folder(ii).name);
end

SNR_15 = 0.5*(SNR_15+SNR_15p);
SNR_20 = 0.5*(SNR_20+SNR_20p);
SNR_25 = 0.5*(SNR_25+SNR_25p);

SNR = cat(1,SNR_15([1:2:35]),SNR_20([1:2:35]),SNR_25([1:2:35]));
SNR_7TS = SNR(:,[2 3 7 9 12])';
SNR_7TL = SNR(:,[11 10 15 4 8])';
SNR_3T = SNR(:,[1 6 14 16 17 18])';

cd('G:\dsi_data_7T_20200901_try')
save('SNR_results_20220527.mat','SNR','SNR_7TS','SNR_7TL','SNR_3T')
fprintf('\n=== SNR COMPUTATION COMPLETE. Saved to SNR_results_20220527.mat ===\n');
