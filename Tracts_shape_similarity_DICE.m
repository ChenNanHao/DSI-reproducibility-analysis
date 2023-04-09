clearvars;clc;close all
%% CST
cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
load('Tracts_resampling100_pre_outremove.mat');
load('Tracts_resampling100_post_outremove.mat');
    % 7T Longer TE : HY LWK MS OK SY [7:11]
    % 7T Shorter TE : IS MA HM HR OS KW [1:6]
    % 3T : CH KHC Peter RL YP YRC [12:17]
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
        cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',a,'_post_',voxelsize(jj,:)))
        aaa = dir('*ROIs.nii*').name;
        Mask = niftiread(aaa);

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

Midslice(1,:) = 144-Midslice(1,:);
Midslice(2,:) = 108-Midslice(2,:);
Midslice(3,:) = 96-Midslice(3,:);

%% Calculate DICE coefficient
for jj = 1:3
    for ii = [1:13 15:18]
        tpre_resam = Tracts_resam_pre_outremove{jj,ii};
        tpost_resam = Tracts_resam_post_outremove{jj,ii};
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
        
for kk = 1:size(tpre_resam,3)
    tpre_resam_ceil(:,:,kk) = ceil(tpre_resam(:,:,kk));
    for kkk = 1:size(tpre_resam_ceil(:,:,kk),2)
        pre_spatial(tpre_resam_ceil(1,kkk,kk),tpre_resam_ceil(2,kkk,kk),tpre_resam_ceil(3,kkk,kk)) = 1;
    end
end
for kk = 1:size(tpost_resam,3)
    tpost_resam_ceil(:,:,kk) = ceil(tpost_resam(:,:,kk));
    for kkk = 1:size(tpost_resam_ceil(:,:,kk),2)
        post_spatial(tpost_resam_ceil(1,kkk,kk),tpost_resam_ceil(2,kkk,kk),tpost_resam_ceil(3,kkk,kk)) = 1;
    end
end

DICE(jj,ii) = 2*size(find(pre_spatial == 1 & post_spatial == 1),1)/...
    (size(find(pre_spatial == 1),1)+size(find(post_spatial == 1),1))*100;
Overlay_voxel(jj,ii) = (size(find(pre_spatial == 1),1)+size(find(post_spatial == 1),1));
    end
    clear tpre_resam_ceil tpost_resam_ceil pre_spatial post_spatial
end

DICE_7TS = mean(DICE(:,[1:3 5:6]),3);
DICE_7TL = mean(DICE(:,7:11),3);
DICE_3T = mean(DICE(:,[12:13 15:18]),3);
DICE_mean = [mean(DICE_7TS,2) mean(DICE_7TL,2) mean(DICE_3T,2)];
DICE_std = [std(DICE_7TS,[],2) std(DICE_7TL,[],2) std(DICE_3T,[],2)];

X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(DICE_mean, 1);
nbars = size(DICE_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, DICE_mean(:,i), DICE_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, DICE_mean);
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 100]);yticks([0:20:100])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' DICE coefficient (%)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% DICE
p_DICE = [];
c_DICE = [];
for ii = 1:3
x=cat(2,DICE_7TS(ii,:),DICE_7TL(ii,:),DICE_3T(ii,:))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_DICE = cat(2,p_DICE,p);
c_DICE= cat(2,c_DICE,c(:,6));
end
c_DICE = c_DICE';

p_DICE_v = [];
c_DICE_v = [];
for ii = 1:3
    if ii == 1
        CC = DICE_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = DICE_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = DICE_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_DICE_v = cat(2,p_DICE_v,p);
c_DICE_v= cat(2,c_DICE_v,c(:,6));
end
c_DICE_v = c_DICE_v';

%% Calculate weighted DICE coefficient
for jj = 1:3
    for ii = [1:13 15:18]
        tpre_resam = Tracts_resam_pre_outremove{jj,ii};
        tpost_resam = Tracts_resam_post_outremove{jj,ii};
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
        
for kk = 1:size(tpre_resam,3)
    tpre_resam_ceil(:,:,kk) = ceil(tpre_resam(:,:,kk));
    for kkk = 1:size(tpre_resam_ceil(:,:,kk),2)
        pre_spatial(tpre_resam_ceil(1,kkk,kk),tpre_resam_ceil(2,kkk,kk),tpre_resam_ceil(3,kkk,kk)) = ...
            pre_spatial(tpre_resam_ceil(1,kkk,kk),tpre_resam_ceil(2,kkk,kk),tpre_resam_ceil(3,kkk,kk))+1;
    end
end
for kk = 1:size(tpost_resam,3)
    tpost_resam_ceil(:,:,kk) = ceil(tpost_resam(:,:,kk));
    for kkk = 1:size(tpost_resam_ceil(:,:,kk),2)
        post_spatial(tpost_resam_ceil(1,kkk,kk),tpost_resam_ceil(2,kkk,kk),tpost_resam_ceil(3,kkk,kk)) = ...
            post_spatial(tpost_resam_ceil(1,kkk,kk),tpost_resam_ceil(2,kkk,kk),tpost_resam_ceil(3,kkk,kk))+1;
    end
end

pre_spatial = pre_spatial/IDX(jj,1,ii);
post_spatial = post_spatial/IDX(jj,2,ii);
[xx,yy,zz] = ind2sub(size(pre_spatial),find(pre_spatial > 0 & post_spatial > 0));

for kkkk = 1:size(xx,1)
    mpre_spatial(kkkk) = pre_spatial(xx(kkkk),yy(kkkk),zz(kkkk));
    mpost_spatial(kkkk) = post_spatial(xx(kkkk),yy(kkkk),zz(kkkk));
end
    
wDICE(jj,ii) = (sum(mpre_spatial)+sum(mpost_spatial))/...
    (sum(pre_spatial,'all')+sum(post_spatial,'all'))*100;
    clear mpre_spatial mpost_spatial
    end
    clear tpre_resam_ceil tpost_resam_ceil pre_spatial post_spatial
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

%%
cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
load('Tracts_centroid_line_pre_outremove.mat')
load('Tracts_centroid_line_post_outremove.mat')

