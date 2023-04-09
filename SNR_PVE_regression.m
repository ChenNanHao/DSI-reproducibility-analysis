clearvars;clc;close all
%% Read files 
cd('E:\dsi_data_7T_20200901_try')
load('ODF results')
load('wDICE CST results')
load('SNR_results_20220527.mat')
cd('PVE')
load('PVE results')

People{1} = 'HM';People{2} = 'HR';People{3} = 'KW';
People{4} = 'IS';People{5} = 'MA';People{6} = 'OS';
People{7} = 'OK';People{8} = 'MS';People{9} = 'SY';
People{10} = 'HY';People{11} = 'LWK';
People{12} = 'CH';People{13} = 'KHC';People{14} = 'Peter';
People{15} = 'RL';People{16}= 'YP';People{17} = 'YRC';
People{18} = 'YR';
for ii = 1:18
        a = People{ii};
    for jj = 1:3
        if jj == 1
            slice = 60;
        elseif jj == 2
            slice = 46;
        else 
            slice = 36;
        end
        voxelsize = ['1p5mm';'2p0mm';'2p5mm'];
        cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',a,'_pre_',voxelsize(jj,:)))
        load('CST_nqa.mat');
        Length(jj,:,ii,1) = length;
        cum_length = [0 cumsum(length)];       
        for kk = 1:size(cum_length,2)-1
            tracts_separate{kk} = tracts(:,cum_length(kk)+1:cum_length(kk+1));
        end
        Tracts.pre{jj,ii,1} = tracts_separate;
        
        cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',a,'_post_',voxelsize(jj,:)))
        load('CST_nqa.mat');
        Length(jj,:,ii,2) = length;
        aaa = dir('*ROIs.nii*').name;
        Mask = niftiread(aaa);
        cum_length = [0 cumsum(length)];       
        for kk = 1:size(cum_length,2)-1
            tracts_separate{kk} = tracts(:,cum_length(kk)+1:cum_length(kk+1));
        end
        Tracts.post{jj,ii,1} = tracts_separate;

        if ii == 13 && jj == 1
        [a1,a2] = find(Mask(:,:,slice/2) == 3);
        Midslice(jj,ii) = unique(a1);
        [b1,b2,b3] = ind2sub(size(Mask),find(Mask == 1));
        Bottomslice(jj,ii) = unique(b3);
        elseif ii == 14 && jj == 2
        [a1,a2] = find(Mask(:,:,slice/2) == 1);
        Midslice(jj,ii) = unique(a1);
        [b1,b2,b3] = ind2sub(size(Mask),find(Mask == 3));
        Bottomslice(jj,ii) = unique(b3);
        else
        [a1,a2] = find(Mask(:,:,slice/2) == 2);
        Midslice(jj,ii) = unique(a1);
        [b1,b2,b3] = ind2sub(size(Mask),find(Mask == 1));
        Bottomslice(jj,ii) = unique(b3);
        end
        disp(append('Finished in ', People{ii},' with voxel size ', voxelsize(jj,:)))
    end
end

length_interval = linspace(0,150,101);length_center = movmean(length_interval,2);
Length(1,:,:,:) = Length(1,:,:,:)*0.75;
Length(2,:,:,:) = Length(2,:,:,:)*1;
Length(3,:,:,:) = Length(3,:,:,:)*1.25;

Bottomslice(1,:,2) = Bottomslice(1,:,1) + 40;
Bottomslice(2,:,2) = Bottomslice(2,:,1) + 30;
Bottomslice(3,:,2) = Bottomslice(3,:,1) + 24;

clear length

%% Calculate PVE in tract bundle
People{1} = 'HM';People{2} = 'HR';People{3} = 'KW';
People{4} = 'IS';People{5} = 'MA';People{6} = 'OS';
People{7} = 'OK';People{8} = 'MS';People{9} = 'SY';
People{10} = 'HY';People{11} = 'LWK';
People{12} = 'CH';People{13} = 'KHC';People{14} = 'Peter';
People{15} = 'RL';People{16}= 'YP';People{17} = 'YRC';
People{18} = 'YR';
    % 7T Longer TE : HY LWK MS OK SY [7:11]
    % 7T Shorter TE : IS MA HM HR OS KW [1:6]
    % 3T : CH KHC Peter RL YP YRC [12:17]
% PVE_repro_15 = abs(1-normalize_vf_pre_15(:,:,:,3) - 1+normalize_vf_post_15(:,:,:,3))./(1-normalize_vf_pre_15(:,:,:,3))*100;
% PVE_repro_20 = abs(1-normalize_vf_pre_20(:,:,:,3) - 1+normalize_vf_post_20(:,:,:,3))./(1-normalize_vf_pre_20(:,:,:,3))*100;
% PVE_repro_25 = abs(1-normalize_vf_pre_25(:,:,:,3) - 1+normalize_vf_post_25(:,:,:,3))./(1-normalize_vf_pre_25(:,:,:,3))*100;

for aa = [1:13 15:18]
people = People{aa};
load(append('E:\dsi_data_7T_20200901_try\1p5mm_results_20220411\ODF_info_15_',people),'singlefibermask','crossfibermask');
singlefibermask_15 = singlefibermask; crossfibermask_15 = crossfibermask; 
load(append('E:\dsi_data_7T_20200901_try\2p0mm_results_20220411\ODF_info_20_',people),'singlefibermask','crossfibermask');
singlefibermask_20 = singlefibermask; crossfibermask_20 = crossfibermask; 
load(append('E:\dsi_data_7T_20200901_try\2p5mm_results_20220411\ODF_info_25_',people),'singlefibermask','crossfibermask');
singlefibermask_25 = singlefibermask; crossfibermask_25 = crossfibermask; 

cd(append('E:\dsi_data_7T_20200901_try\PVE\',people))
vf_pre_15 = niftiread('vf_pre_15_norm.nii.gz');
vf_pre_20 = niftiread('vf_pre_20_norm.nii.gz');
vf_pre_25 = niftiread('vf_pre_25_norm.nii.gz');
normalize_vf_pre_15 = vf_pre_15./sum(vf_pre_15,4);
normalize_vf_pre_20 = vf_pre_20./sum(vf_pre_20,4);
normalize_vf_pre_25 = vf_pre_25./sum(vf_pre_25,4);
vf_post_15 = niftiread('vf_post_15_norm.nii.gz');
vf_post_20 = niftiread('vf_post_20_norm.nii.gz');
vf_post_25 = niftiread('vf_post_25_norm.nii.gz');
normalize_vf_post_15 = vf_post_15./sum(vf_post_15,4);
normalize_vf_post_20 = vf_post_20./sum(vf_post_20,4);
normalize_vf_post_25 = vf_post_25./sum(vf_post_25,4);
PVE_SF_15 = 0.5*(vf_pre_15(:,:,:,1)+vf_pre_15(:,:,:,2)+vf_post_15(:,:,:,1)+vf_post_15(:,:,:,2));
PVE_SF_20 = 0.5*(vf_pre_20(:,:,:,1)+vf_pre_20(:,:,:,2)+vf_post_20(:,:,:,1)+vf_post_20(:,:,:,2));
PVE_SF_25 = 0.5*(vf_pre_25(:,:,:,1)+vf_pre_25(:,:,:,2)+vf_post_25(:,:,:,1)+vf_post_25(:,:,:,2));
PVE_15 = 0.5*(1-normalize_vf_pre_15(:,:,:,3) + 1-normalize_vf_post_15(:,:,:,3));
PVE_20 = 0.5*(1-normalize_vf_pre_20(:,:,:,3) + 1-normalize_vf_post_20(:,:,:,3));
PVE_25 = 0.5*(1-normalize_vf_pre_25(:,:,:,3) + 1-normalize_vf_post_25(:,:,:,3));

slice = 60;
slice_interval = 40;
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_pre_1p5mm'))
Pre_ROI = niftiread('CST_ROI_pre.nii'); Pre_ROI(Pre_ROI == 0) = nan; 
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_post_1p5mm'))
Post_ROI = niftiread('CST_ROI_post.nii'); Post_ROI(Post_ROI == 0) = nan; 
PVE_15_CST = PVE_15.*Pre_ROI.*Post_ROI;
PVE_Whole(aa,1) = mean(PVE_15,'all','omitnan');
PVE_homo(aa,1) = mean(PVE_15(:,:,Bottomslice(1,aa):(Bottomslice(1,aa)+slice_interval)),'all','omitnan');
PVE_hetero(aa,1) = mean(PVE_15(:,:,[1:Bottomslice(1,aa)-1 (Bottomslice(1,aa)+slice_interval+1):slice]),'all','omitnan');

slice = 46;
slice_interval = 30;
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_pre_2p0mm'))
Pre_ROI = niftiread('CST_ROI_pre.nii'); Pre_ROI(Pre_ROI == 0) = nan; 
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_post_2p0mm'))
Post_ROI = niftiread('CST_ROI_post.nii'); Post_ROI(Post_ROI == 0) = nan; 
PVE_20_CST = PVE_20.*Pre_ROI.*Post_ROI;
PVE_Whole(aa,2) = mean(PVE_20,'all','omitnan');
PVE_homo(aa,2) = mean(PVE_20(:,:,Bottomslice(2,aa):(Bottomslice(2,aa)+slice_interval)),'all','omitnan');
PVE_hetero(aa,2) = mean(PVE_20(:,:,[1:Bottomslice(2,aa)-1 (Bottomslice(2,aa)+slice_interval+1):slice]),'all','omitnan');

slice = 36;
slice_interval = 24;
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_pre_2p5mm'))
Pre_ROI = niftiread('CST_ROI_pre.nii'); Pre_ROI(Pre_ROI == 0) = nan; 
cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people,'_post_2p5mm'))
Post_ROI = niftiread('CST_ROI_post.nii'); Post_ROI(Post_ROI == 0) = nan; 
PVE_25_CST = PVE_25.*Pre_ROI.*Post_ROI;
PVE_Whole(aa,3) = mean(PVE_25,'all','omitnan');
PVE_homo(aa,3) = mean(PVE_25(:,:,Bottomslice(3,aa):(Bottomslice(3,aa)+slice_interval)),'all','omitnan');
PVE_hetero(aa,3) = mean(PVE_25(:,:,[1:Bottomslice(3,aa)-1 (Bottomslice(3,aa)+slice_interval+1):slice]),'all','omitnan');

disp(append('Finished in ', People{aa}))
end
PVE_Whole = PVE_Whole*100;
PVE_homo = PVE_homo*100;
PVE_hetero = PVE_hetero*100;

%% Regression
GLM_Model = 'Y ~ X';

for kk = 1:3
    tbl = table(SNR(kk,[1:3 5:13 15:18])',PVE_single([1:3 5:13 15:18],kk),'VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    SNR_fitted(kk,:) = mdl.Fitted.Response; 
end
SNR_regressed_single = SNR(:,[1:3 5:13 15:18])-SNR_fitted;
for kk = 1:3
    tbl = table(SNR(kk,[1:3 5:13 15:18])',PVE_crossing([1:3 5:13 15:18],kk),'VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    SNR_fitted(kk,:) = mdl.Fitted.Response; 
end
SNR_regressed_crossing = SNR(:,[1:3 5:13 15:18])-SNR_fitted;
for kk = 1:3
    tbl = table(SNR(kk,[1:3 5:13 15:18])',PVE_Whole([1:3 5:13 15:18],kk),'VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    SNR_fitted(kk,:) = mdl.Fitted.Response; 
end
SNR_regressed_Whole = SNR(:,[1:3 5:13 15:18])-SNR_fitted;
for kk = 1:3
    tbl = table(SNR(kk,[1:3 5:13 15:18])',PVE_homo([1:3 5:13 15:18],kk),'VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    SNR_fitted(kk,:) = mdl.Fitted.Response; 
end
SNR_regressed_homo = SNR(:,[1:3 5:13 15:18])-SNR_fitted;
for kk = 1:3
    tbl = table(SNR(kk,[1:3 5:13 15:18])',PVE_hetero([1:3 5:13 15:18],kk),'VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    SNR_fitted(kk,:) = mdl.Fitted.Response; 
end
SNR_regressed_hetero = SNR(:,[1:3 5:13 15:18])-SNR_fitted;


for kk = 1:3
    tbl = table(PVE_single([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    PVE_fitted(kk,:) = mdl.Fitted.Response; 
end
PVE_regressed_single = PVE_single([1:3 5:13 15:18],:)-PVE_fitted';

for kk = 1:3
    tbl = table(PVE_crossing([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    PVE_fitted(kk,:) = mdl.Fitted.Response; 
end
PVE_regressed_crossing = PVE_crossing([1:3 5:13 15:18],:)-PVE_fitted';
for kk = 1:3
    tbl = table(PVE_Whole([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    PVE_fitted(kk,:) = mdl.Fitted.Response; 
end
PVE_regressed_Whole = PVE_Whole([1:3 5:13 15:18],:)-PVE_fitted';
for kk = 1:3
    tbl = table(PVE_homo([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    PVE_fitted(kk,:) = mdl.Fitted.Response; 
end
PVE_regressed_homo = PVE_homo([1:3 5:13 15:18],:)-PVE_fitted';
for kk = 1:3
    tbl = table(PVE_hetero([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'Y','X'});
    mdl = fitglm(tbl,GLM_Model);
    PVE_fitted(kk,:) = mdl.Fitted.Response; 
end
PVE_regressed_hetero = PVE_hetero([1:3 5:13 15:18],:)-PVE_fitted';

%% Plot regressed values
PVE_7TS = PVE_regressed_single(1:5,:);
PVE_7TL = PVE_regressed_single(6:10,:);
PVE_3T = PVE_regressed_single(11:16,:);
PVE_mean = [mean(PVE_7TS)' mean(PVE_7TL)' mean(PVE_3T)'];
PVE_std = [std(PVE_7TS)' std(PVE_7TL)' std(PVE_3T)'];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(PVE_mean, 1);
nbars = size(PVE_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
figure
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, PVE_mean(:,i), PVE_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, PVE_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 0.08]);yticks([0:0.02:1])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Non-WM signal fraction ','fontname','calibri','fontsize',18)

%% Linear regression (single-fiber)
% PVE vs. repro
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for kk = 1:3
plot(PVE_regressed_single(:,kk),Single_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_regressed_single(:,kk)),max(PVE_regressed_single(:,kk)),100);
[P,S] = polyfit(PVE_regressed_single(:,kk),Single_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 15]);yticks([0:5:40]);
% xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end
x = linspace(min(PVE_regressed_single,[],'all'),max(PVE_regressed_single,[],'all'),100);
[P,S] = polyfit(PVE_regressed_single,Single_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(PVE_regressed_single(:,kk),Single_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_regressed_single(:,kk)),max(PVE_regressed_single(:,kk)),100);
[P,S] = polyfit(PVE_regressed_single(:,kk),Single_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([1 3]);yticks([1:0.5:3]);
% xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
end
x = linspace(min(PVE_regressed_single,[],'all'),max(PVE_regressed_single,[],'all'),100);
[P,S] = polyfit(PVE_regressed_single,Single_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

% SNR vs. repro
for kk = 1:3
plot(SNR_regressed_single(kk,:),Single_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR_regressed_single(kk,:)),max(SNR_regressed_single(kk,:)),100);
[P,S] = polyfit(SNR_regressed_single(kk,:)',Single_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 15]);yticks([0:5:40]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR_regressed_single,[],'all'),max(SNR_regressed_single,[],'all'),100);
[P,S] = polyfit(SNR_regressed_single',Single_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR_regressed_single(kk,:),Single_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR_regressed_single(kk,:)),max(SNR_regressed_single(kk,:)),100);
[P,S] = polyfit(SNR_regressed_single(kk,:)',Single_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([1 3]);yticks([1:0.5:3]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR_regressed_single,[],'all'),max(SNR_regressed_single,[],'all'),100);
[P,S] = polyfit(SNR_regressed_single',Single_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

%% Linear regression (crossing-fiber)
% PVE vs. repro
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for kk = 1:3
plot(PVE_regressed_crossing(:,kk),Crossing_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_regressed_crossing(:,kk)),max(PVE_regressed_crossing(:,kk)),100);
[P,S] = polyfit(PVE_regressed_crossing(:,kk),Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([15 40]);yticks([15:5:40]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end
x = linspace(min(PVE_regressed_crossing,[],'all'),max(PVE_regressed_crossing,[],'all'),100);
[P,S] = polyfit(PVE_regressed_crossing,Crossing_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(PVE_regressed_crossing(:,kk),Crossing_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_regressed_crossing(:,kk)),max(PVE_regressed_crossing(:,kk)),100);
[P,S] = polyfit(PVE_regressed_crossing(:,kk),Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 2]);yticks([0:0.5:2]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
end
x = linspace(min(PVE_regressed_crossing,[],'all'),max(PVE_regressed_crossing,[],'all'),100);
[P,S] = polyfit(PVE_regressed_crossing,Crossing_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

% SNR vs. repro
for kk = 1:3
plot(SNR_regressed_single(kk,:),Crossing_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR_regressed_single(kk,:)),max(SNR_regressed_single(kk,:)),100);
[P,S] = polyfit(SNR_regressed_single(kk,:)',Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([15 40]);yticks([15:5:40]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR_regressed_single,[],'all'),max(SNR_regressed_single,[],'all'),100);
[P,S] = polyfit(SNR_regressed_single',Crossing_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR_regressed_single(kk,:),Crossing_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR_regressed_single(kk,:)),max(SNR_regressed_single(kk,:)),100);
[P,S] = polyfit(SNR_regressed_single(kk,:)',Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 2]);yticks([0:0.5:2]);
% xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR_regressed_single,[],'all'),max(SNR_regressed_single,[],'all'),100);
[P,S] = polyfit(SNR_regressed_single',Crossing_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

%% Linear regression (Tract-based repro)
% PVE vs. repro
for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE_homo(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE_homo(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE_homo(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE_hetero(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE_hetero(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE_hetero(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

% SNR vs. repro
for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE_homo(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE_homo(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE_homo(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),wDICE_hetero(kk,[1:3 5:13 15:18]),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18]),wDICE_hetero(kk,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([50 100]);yticks([0:10:100]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' wDICE(%) ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18]),wDICE_hetero(:,[1:3 5:13 15:18]),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on
