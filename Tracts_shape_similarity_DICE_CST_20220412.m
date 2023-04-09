clearvars;clc;close all
%% Read Tracts file
    % 7T Longer TE : HY LWK MS OK SY [1:5]
    % 7T Shorter TE : IS MA HM HR OS KW [6:11]
    % 3T : CH KHC Peter RL YP YRC [12:18]

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

cd('E:\dsi_data_7T_20200901_try\SRC_Batch\YR_post_1p5mm')
info_15 = niftiinfo('rcT1_pre_brain.nii.gz');
cd('E:\dsi_data_7T_20200901_try\SRC_Batch\YR_post_2p0mm')
info_20 = niftiinfo('rcT1_pre_brain.nii.gz');
cd('E:\dsi_data_7T_20200901_try\SRC_Batch\YR_post_2p5mm')
info_25 = niftiinfo('rcT1_pre_brain.nii.gz');

%% Calculate weighted DICE coefficient (whole tract bundles)
for jj = 1:3
    for ii = [1:13 15:18]
        tpre_resam = Tracts.pre{jj,ii};
        tpost_resam = Tracts.post{jj,ii};
        if jj == 1
            pre_spatial = zeros(144,144,60);
            post_spatial = zeros(144,144,60);
        elseif jj == 2
            pre_spatial = zeros(108,108,46);
            post_spatial = zeros(108,108,46);
        else
            pre_spatial = zeros(96,96,36);
            post_spatial = zeros(96,96,36);
        end
        
for kk = 1:size(tpre_resam,2)
    if find(tpre_resam{kk} <= 0)
        [idx_x,idx_y] = find(tpre_resam{kk} <= 0);
        tpre_resam{kk}(idx_x,idx_y) = tpre_resam{kk}(idx_x,idx_y)+1;
    end
    tpre_resam_ceil{kk} = ceil(tpre_resam{kk});
    for kkk = 1:size(tpre_resam_ceil{kk},2)
        pre_spatial(tpre_resam_ceil{kk}(1,kkk),tpre_resam_ceil{kk}(2,kkk),tpre_resam_ceil{kk}(3,kkk)) = ...
            pre_spatial(tpre_resam_ceil{kk}(1,kkk),tpre_resam_ceil{kk}(2,kkk),tpre_resam_ceil{kk}(3,kkk))+1;
    end
end
for kk = 1:size(tpost_resam,2)
    if find(tpost_resam{kk} <= 0)
        [idx_x,idx_y] = find(tpost_resam{kk} <= 0);
        tpost_resam{kk}(idx_x,idx_y) = tpost_resam{kk}(idx_x,idx_y)+1;
    end
    tpost_resam_ceil{kk} = ceil(tpost_resam{kk});
    for kkk = 1:size(tpost_resam_ceil{kk},2)
        post_spatial(tpost_resam_ceil{kk}(1,kkk),tpost_resam_ceil{kk}(2,kkk),tpost_resam_ceil{kk}(3,kkk)) = ...
            post_spatial(tpost_resam_ceil{kk}(1,kkk),tpost_resam_ceil{kk}(2,kkk),tpost_resam_ceil{kk}(3,kkk))+1;
    end
end

pre_spatial = pre_spatial/10000;
post_spatial = post_spatial/10000;
[xx,yy,zz] = ind2sub(size(pre_spatial),find(pre_spatial > 0 & post_spatial > 0));

for kkkk = 1:size(xx,1)
    mpre_spatial(kkkk) = pre_spatial(xx(kkkk),yy(kkkk),zz(kkkk));
    mpost_spatial(kkkk) = post_spatial(xx(kkkk),yy(kkkk),zz(kkkk));
end
    
wDICE(jj,ii) = (sum(mpre_spatial)+sum(mpost_spatial))/...
    (sum(pre_spatial,'all')+sum(post_spatial,'all'))*100;
    clear mpre_spatial mpost_spatial
    
    % Save ROI
    pre_spatial_roi = pre_spatial; pre_spatial_roi(pre_spatial_roi > 0) = 1; pre_spatial_roi(pre_spatial_roi <= 0) = 0; 
    post_spatial_roi = post_spatial; post_spatial_roi(post_spatial_roi > 0) = 1; post_spatial_roi(post_spatial_roi <= 0) = 0; 
    a = People{ii};
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',a,'_pre_',voxelsize(jj,:)))
    if jj == 1
        info = info_15;
    elseif jj == 2
        info = info_20;
    else
        info = info_25; 
    end
    niftiwrite(single(pre_spatial_roi),'CST_ROI_pre.nii',info)
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',a,'_post_',voxelsize(jj,:)))
    niftiwrite(single(post_spatial_roi),'CST_ROI_post.nii',info)
    end
    clear tpre_resam_ceil tpost_resam_ceil pre_spatial post_spatial pre_spatial_roi post_spatial_roi
end

wDICE_7TS = mean(wDICE(:,[1:3 5:6]),3);
wDICE_7TL = mean(wDICE(:,7:11),3);
wDICE_3T = mean(wDICE(:,[12:13 15:18]),3);
wDICE_mean = [mean(wDICE_7TS,2) mean(wDICE_7TL,2) mean(wDICE_3T,2)];
wDICE_std = [std(wDICE_7TS,[],2) std(wDICE_7TL,[],2) std(wDICE_3T,[],2)];

X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(wDICE_mean, 1);
nbars = size(wDICE_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, wDICE_mean(:,i), wDICE_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, wDICE_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
plot(linspace(0.5,3.5,100),ones(1,100)*70,'r-','linewidth',2)
xticks([1 2 3]);ylim([0 120]);yticks([0:40:120])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' wDICE coefficient (%)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% wDICE
p_wDICE = [];
c_wDICE = [];
for ii = 1:3
x=cat(2,wDICE_7TS(ii,:),wDICE_7TL(ii,:),wDICE_3T(ii,:))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE = cat(2,p_wDICE,p);
c_wDICE= cat(2,c_wDICE,c(:,6));
end
c_wDICE = c_wDICE';

p_wDICE_v = [];
c_wDICE_v = [];
for ii = 1:3
    if ii == 1
        CC = wDICE_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = wDICE_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = wDICE_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE_v = cat(2,p_wDICE_v,p);
c_wDICE_v= cat(2,c_wDICE_v,c(:,6));
end
c_wDICE_v = c_wDICE_v';

%% Calculate weighted DICE coefficient (homo. and hetero.)
for jj = 1:3
    for ii = [1:13 15:18]
        if jj == 1
            slice = 60;
            slice_interval = 40;
        elseif jj == 2
            slice = 46;
            slice_interval = 30;
        else 
            slice = 36;
            slice_interval = 24;
        end
        
        tpre_resam = Tracts.pre{jj,ii};
        tpost_resam = Tracts.post{jj,ii};
        if jj == 1
            pre_spatial = zeros(144,144,60);
            post_spatial = zeros(144,144,60);
        elseif jj == 2
            pre_spatial = zeros(108,108,46);
            post_spatial = zeros(108,108,46);
        else
            pre_spatial = zeros(96,96,36);
            post_spatial = zeros(96,96,36);
        end
        
for kk = 1:size(tpre_resam,2)
    if find(tpre_resam{kk} <= 0)
        [idx_x,idx_y] = find(tpre_resam{kk} <= 0);
        tpre_resam{kk}(idx_x,idx_y) = tpre_resam{kk}(idx_x,idx_y)+1;
    end
    tpre_resam_ceil{kk} = ceil(tpre_resam{kk});
    for kkk = 1:size(tpre_resam_ceil{kk},2)
        pre_spatial(tpre_resam_ceil{kk}(1,kkk),tpre_resam_ceil{kk}(2,kkk),tpre_resam_ceil{kk}(3,kkk)) = ...
            pre_spatial(tpre_resam_ceil{kk}(1,kkk),tpre_resam_ceil{kk}(2,kkk),tpre_resam_ceil{kk}(3,kkk))+1;
    end
end
for kk = 1:size(tpost_resam,2)
    if find(tpost_resam{kk} <= 0)
        [idx_x,idx_y] = find(tpost_resam{kk} <= 0);
        tpost_resam{kk}(idx_x,idx_y) = tpost_resam{kk}(idx_x,idx_y)+1;
    end
    tpost_resam_ceil{kk} = ceil(tpost_resam{kk});
    for kkk = 1:size(tpost_resam_ceil{kk},2)
        post_spatial(tpost_resam_ceil{kk}(1,kkk),tpost_resam_ceil{kk}(2,kkk),tpost_resam_ceil{kk}(3,kkk)) = ...
            post_spatial(tpost_resam_ceil{kk}(1,kkk),tpost_resam_ceil{kk}(2,kkk),tpost_resam_ceil{kk}(3,kkk))+1;
    end
end

% Seperate into homogeneous and heterogeneous regions
% Homogeneous : middle part of bundles with 60 mm 
% Heterogeneous : Total bundle - homogeneous part 
pre_spatial = pre_spatial/10000;
post_spatial = post_spatial/10000;
pre_spatial_homo = pre_spatial(:,:,Bottomslice(jj,ii):(Bottomslice(jj,ii)+slice_interval));
post_spatial_homo = post_spatial(:,:,Bottomslice(jj,ii):(Bottomslice(jj,ii)+slice_interval));
pre_spatial_hetero = pre_spatial(:,:,[1:Bottomslice(jj,ii)-1 (Bottomslice(jj,ii)+slice_interval+1):slice]);
post_spatial_hetero = post_spatial(:,:,[1:Bottomslice(jj,ii)-1 (Bottomslice(jj,ii)+slice_interval+1):slice]);

[xx_homo,yy_homo,zz_homo] = ind2sub(size(pre_spatial_homo),find(pre_spatial_homo > 0 & post_spatial_homo > 0));
[xx_hetero,yy_hetero,zz_hetero] = ind2sub(size(pre_spatial_hetero),find(pre_spatial_hetero > 0 & post_spatial_hetero > 0));

for kkkk = 1:size(xx_homo,1)
    mpre_spatial_homo(kkkk) = pre_spatial_homo(xx_homo(kkkk),yy_homo(kkkk),zz_homo(kkkk));
    mpost_spatial_homo(kkkk) = post_spatial_homo(xx_homo(kkkk),yy_homo(kkkk),zz_homo(kkkk));
end
for kkkk = 1:size(xx_hetero,1)
    mpre_spatial_hetero(kkkk) = pre_spatial_hetero(xx_hetero(kkkk),yy_hetero(kkkk),zz_hetero(kkkk));
    mpost_spatial_hetero(kkkk) = post_spatial_hetero(xx_hetero(kkkk),yy_hetero(kkkk),zz_hetero(kkkk));
end

wDICE_homo(jj,ii) = (sum(mpre_spatial_homo)+sum(mpost_spatial_homo))/...
    (sum(pre_spatial_homo,'all')+sum(post_spatial_homo,'all'))*100;
wDICE_hetero(jj,ii) = (sum(mpre_spatial_hetero)+sum(mpost_spatial_hetero))/...
    (sum(pre_spatial_hetero,'all')+sum(post_spatial_hetero,'all'))*100;
    clear mpre_spatial_homo mpost_spatial_homo mpre_spatial_hetero mpost_spatial_hetero
    
    end
    clear tpre_resam_ceil tpost_resam_ceil pre_spatial post_spatial pre_spatial_homo post_spatial_homo pre_spatial_hetero post_spatial_hetero 
end

wDICE_7TS = mean(wDICE_homo(:,[1:3 5:6]),3);
wDICE_7TL = mean(wDICE_homo(:,7:11),3);
wDICE_3T = mean(wDICE_homo(:,[12:13 15:18]),3);
wDICE_mean = [mean(wDICE_7TS,2) mean(wDICE_7TL,2) mean(wDICE_3T,2)];
wDICE_std = [std(wDICE_7TS,[],2) std(wDICE_7TL,[],2) std(wDICE_3T,[],2)];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(wDICE_mean, 1);
nbars = size(wDICE_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, wDICE_mean(:,i), wDICE_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, wDICE_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
plot(linspace(0.5,3.5,100),ones(1,100)*70,'r-','linewidth',2)
xticks([1 2 3]);ylim([0 120]);yticks([0:40:120])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' wDICE coefficient (%)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% wDICE
p_wDICE = [];
c_wDICE = [];
for ii = 1:3
x=cat(2,wDICE_7TS(ii,:),wDICE_7TL(ii,:),wDICE_3T(ii,:))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE = cat(2,p_wDICE,p);
c_wDICE= cat(2,c_wDICE,c(:,6));
end
c_wDICE = c_wDICE';

p_wDICE_v = [];
c_wDICE_v = [];
for ii = 1:3
    if ii == 1
        CC = wDICE_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = wDICE_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = wDICE_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE_v = cat(2,p_wDICE_v,p);
c_wDICE_v= cat(2,c_wDICE_v,c(:,6));
end
c_wDICE_v = c_wDICE_v';


wDICE_7TS = mean(wDICE_hetero(:,[1:3 5:6]),3);
wDICE_7TL = mean(wDICE_hetero(:,7:11),3);
wDICE_3T = mean(wDICE_hetero(:,[12:13 15:18]),3);
wDICE_mean = [mean(wDICE_7TS,2) mean(wDICE_7TL,2) mean(wDICE_3T,2)];
wDICE_std = [std(wDICE_7TS,[],2) std(wDICE_7TL,[],2) std(wDICE_3T,[],2)];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(wDICE_mean, 1);
nbars = size(wDICE_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, wDICE_mean(:,i), wDICE_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, wDICE_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
plot(linspace(0.5,3.5,100),ones(1,100)*70,'r-','linewidth',2)
xticks([1 2 3]);ylim([0 100]);yticks([0:40:120])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' wDICE coefficient (%)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% wDICE
p_wDICE = [];
c_wDICE = [];
for ii = 1:3
x=cat(2,wDICE_7TS(ii,:),wDICE_7TL(ii,:),wDICE_3T(ii,:))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE = cat(2,p_wDICE,p);
c_wDICE= cat(2,c_wDICE,c(:,6));
end
c_wDICE = c_wDICE';

p_wDICE_v = [];
c_wDICE_v = [];
for ii = 1:3
    if ii == 1
        CC = wDICE_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = wDICE_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = wDICE_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_wDICE_v = cat(2,p_wDICE_v,p);
c_wDICE_v= cat(2,c_wDICE_v,c(:,6));
end
c_wDICE_v = c_wDICE_v';