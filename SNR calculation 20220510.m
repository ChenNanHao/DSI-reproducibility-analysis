clearvars;clc;close all
%%
cd('E:\dsi_data_7T_20200901_try\SRC_Batch')
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
cd('E:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*1p5mm'); 
for ii = [1:2:size(folder,1)-1]
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(load_untouch_nii('c2rcT1_pre.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif ii == 9
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    else
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    end
    SNR_15(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_15p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    disp(append('Finished SNR computation at ',folder(ii).name))
end
figure;
% Firstly, use gray scle on anatomical images
ax1 = axes;
imagesc(imrotate((V1(:,:,35,2)-V1(:,:,35,1)),-90));colorbar
colormap(ax1,'gray'); % choose this figure layer into gray scale
% Secondly, use colormap on onerlay images
ax2 = axes;
imagesc(ax2,imrotate(WMmask(:,:,35),-90)*0.4,'alphadata',0.5);
colormap(ax2,'hot');
caxis(ax2,[0 1]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;


%% 2p0mm SNR
cd('E:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*2p0mm'); 
People{1} = 'NIPS_DSI_7T\HM';People{2} = 'NIPS_DSI_7T\HR';People{3} = 'NIPS_DSI_7T\KW';
People{4} = 'NIPS_DSI_7T\IS';People{5} = 'NIPS_DSI_7T\MA';People{6} = 'NIPS_DSI_7T\OS';
People{7} = 'NIPS_DSI_7T\OK';People{8} = 'NIPS_DSI_7T\MS';People{9} = 'NIPS_DSI_7T\SY';
People{10} = 'NIPS_DSI_7T\HY';People{11} = 'NIPS_DSI_7T\LWK';
People{12} = 'PRISMA3T_NH\CH';People{13} = 'PRISMA3T_NH\KHC';People{14} = 'PRISMA3T_NH\Peter';
People{15} = 'PRISMA3T_NH\RL';People{16} = 'PRISMA3T_NH\YP';People{17} = 'PRISMA3T_NH\YRC';
People{18} = 'PRISMA3T_NH\YR';

for ii = [1:2:size(folder,1)-1]
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(load_untouch_nii('c2rcT1_pre.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif ii == 9
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif ii == 23
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    else
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    end
    SNR_20(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_20p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    disp(append('Finished SNR computation at ',folder(ii).name))
end

%% 2p5mm SNR
cd('E:\dsi_data_7T_20200901_try\SRC_Batch')
folder = dir('*2p5mm'); 
for ii = [1:2:size(folder,1)-1]
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii).name))
    V = niftiread('data.nii.gz');
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',folder(ii+1).name))
    V1 = niftiread('data.nii.gz');
    if ii == 11
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif ii == 9
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif ii == 23
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    else
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    end
    SNR_25(ii) = mean((0.5*(V(:,:,:,1)+V(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V(:,:,:,1)-V(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    SNR_25p(ii) = mean((0.5*(V1(:,:,:,1)+V1(:,:,:,2))).*WMmask(:,:,:),'all','omitnan')./std((V1(:,:,:,1)-V1(:,:,:,2)).*WMmask(:,:,:),[],'all','omitnan');
    disp(append('Finished SNR computation at ',folder(ii).name))
end
SNR_15 = 0.5*(SNR_15+SNR_15p);
SNR_20 = 0.5*(SNR_20+SNR_20p);
SNR_25 = 0.5*(SNR_25+SNR_25p);

SNR = cat(1,SNR_15([1:2:35]),SNR_20([1:2:35]),SNR_25([1:2:35])); 
SNR_7TS = SNR(:,[2 3 7 9 12])';
SNR_7TL = SNR(:,[11 10 15 4 8])';
SNR_3T = SNR(:,[1 6 14 16 17 18])';
SNR_mean = [mean(SNR_7TS)' mean(SNR_7TL)' mean(SNR_3T)'];
SNR_std = [std(SNR_7TS)' std(SNR_7TL)' std(SNR_3T)'];

X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(SNR_mean, 1);
nbars = size(SNR_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
figure
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, SNR_mean(:,i), SNR_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, SNR_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 45]);yticks([0:15:45])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' SNR (a.u.) ','fontname','calibri','fontsize',18)

% Group comparison
p_SNR = [];
c_SNR = [];
for ii = 1:3
x=cat(1,SNR_7TS(:,ii),SNR_7TL(:,ii),SNR_3T(:,ii));
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SNR = cat(2,p_SNR,p);
c_SNR= cat(2,c_SNR,c(:,6));
end
c_SNR = c_SNR';

p_SNR_v = [];
c_SNR_v = [];
for ii = 1:3
    if ii == 1
        CC = SNR_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = SNR_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = SNR_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SNR_v = cat(2,p_SNR_v,p);
c_SNR_v= cat(2,c_SNR_v,c(:,6));
end
c_SNR_v = c_SNR_v';

save('SNR_results_20220527.mat','SNR','SNR_7TS','SNR_7TL','SNR_3T')
