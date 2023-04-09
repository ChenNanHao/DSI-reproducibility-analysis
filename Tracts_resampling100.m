clearvars;clc;close all
%% Read Tracts file
    % 7T Longer TE : HY LWK MS OK SY [1:5]
    % 7T Shorter TE : IS MA HM HR OS KW [6:11]
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

%% Truncate certain range of tract bundles
for ii = 1:18
    for jj = 1:3
        tracts_pre = Tracts.pre{jj,ii};
        tracts_post = Tracts.post{jj,ii};
        for kk = 1:10000
            % Pre
            t = tracts_pre{1,kk};
            truncate_index = find(t(3,:) > Bottomslice(jj,ii,1) & t(3,:) < Bottomslice(jj,ii,2));
            tracts_pre1{1,kk} = t(:,truncate_index);
            % Post
            t = tracts_post{1,kk};
            truncate_index1 = find(t(3,:) > Bottomslice(jj,ii,1) & t(3,:) < Bottomslice(jj,ii,2));
            tracts_post1{1,kk} = t(:,truncate_index1);
        end
        Tracts_truncate.pre{jj,ii} = tracts_pre1;
        Tracts_truncate.post{jj,ii} = tracts_post1;
    end
end

jj = 1;
for ii = 1:18
Tractpre = Tracts_truncate.pre{jj,ii};
Tractpost = Tracts_truncate.post{jj,ii};
for kk = 1:10000
    nodes_pre(ii,kk) = size(Tractpre{1,kk},2);
    nodes_post(ii,kk) = size(Tractpost{1,kk},2);
end
end
plot(nodes_pre','b')
hold on
plot(nodes_post','r')
ylim([0 100])

%% Plot tract bundles
for jj = 1:3
figure
tracts_pre = Tracts_truncate.pre{jj,1};
tracts_post = Tracts_truncate.post{jj,1};

for ii = 1:1000
    tpre = tracts_pre{1,ii};
    tpost = tracts_post{1,ii};
    plot3(tpre(1,:),tpre(2,:),tpre(3,:),'b-')
    hold on
    plot3(Tracts_resam.pre{jj,1}(1,:,ii),Tracts_resam.pre{jj,1}(2,:,ii),Tracts_resam.pre{jj,1}(3,:,ii),'r-')
end
end

%% Resample tract nodes into same nodes number 
for ii = 1:18
    for jj = 1:3
        Tractpre = Tracts_truncate.pre{jj,ii};
        Tractpost = Tracts_truncate.post{jj,ii};
        if jj == 1
            slice_interval = 40;
        elseif jj == 2
            slice_interval = 30;
        else 
            slice_interval = 24;
        end
       % Create parallel pool
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolopen = parpool(6);
            poolopenflag = 1;
        else
            poolopenflag = 1;
        end
        N = 10000;
        parfor_progress(N);
        tic
        parfor kk = 1:10000
            tpre = Tractpre{1,kk};
            tpre_cutindex = find((tpre(3,:) > Bottomslice(jj,ii,1)) + (tpre(3,:) < Bottomslice(jj,ii,2)) == 2);
            tpre_cut = tpre(:,tpre_cutindex);
            tpre_resam(:,:,kk) = interparc(slice_interval,tpre_cut(1,:)', tpre_cut(2,:)', tpre_cut(3,:)')';
%             Tractpre_resam{1,1}{1,kk} = tpre_resam;
            
            tpost = Tractpost{1,kk};
            tpost_cutindex = find((tpost(3,:) > Bottomslice(jj,ii,1)) + (tpost(3,:) < Bottomslice(jj,ii,2)) == 2);
            tpost_cut = tpost(:,tpost_cutindex);
            tpost_resam(:,:,kk) = interparc(slice_interval,tpost_cut(1,:)', tpost_cut(2,:)', tpost_cut(3,:)')';
%             Tractpost_resam{1,1}{1,kk} = tpost_resam;
            % Print progress
            parfor_progress;
        end
        parfor_progress(0);
        if poolopenflag == 1
            delete(poolopen)
        end
        toc
        Tracts_resam.pre{jj,ii} = tpre_resam;
        Tracts_resam.post{jj,ii} = tpost_resam;
        clear tpre_resam tpost_resam
    end
end

cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
Tracts_resam_pre = Tracts_resam.pre;
save('Tracts_resampling100_pre.mat','Tracts_resam_pre')
Tracts_resam_post = Tracts_resam.post;
save('Tracts_resampling100_post.mat','Tracts_resam_post')

%% Calculate minimum average direct-flip distance (MDF) and extract centroid path
for ii = 1:18
    for jj = 1:3
        tpre_resam = Tracts_resam.pre{jj,ii};
        tpost_resam = Tracts_resam.post{jj,ii};
        for kk = 1:10000
            Left_cluster_pre(:,:,kk) = tpre_resam(1,:,kk) <= Midslice(1,jj);
            Right_cluster_pre(:,:,kk) = tpre_resam(1,:,kk) > Midslice(1,jj);
            Left_cluster_post(:,:,kk) = tpost_resam(1,:,kk) <= Midslice(1,jj);
            Right_cluster_post(:,:,kk) = tpost_resam(1,:,kk) > Midslice(1,jj);
        end
        Left_cluster_pre = double(Left_cluster_pre);Left_cluster_pre(Left_cluster_pre == 0) = nan;
        Right_cluster_pre = double(Right_cluster_pre);Right_cluster_pre(Right_cluster_pre == 0) = nan;
        mtpre_resam_left(:,:,jj) = mean(tpre_resam.*Left_cluster_pre,3,'omitnan');
        mtpre_resam_right(:,:,jj) = mean(tpre_resam.*Right_cluster_pre,3,'omitnan');
        
        Left_cluster_post = double(Left_cluster_post);Left_cluster_post(Left_cluster_post == 0) = nan;
        Right_cluster_post = double(Right_cluster_post);Right_cluster_post(Right_cluster_post == 0) = nan;
        mtpost_resam_left(:,:,jj) = mean(tpost_resam.*Left_cluster_post,3,'omitnan');
        mtpost_resam_right(:,:,jj) = mean(tpost_resam.*Right_cluster_post,3,'omitnan');
        
        Mean_centroid_pre{jj,1,ii} = mtpre_resam_left;
        Mean_centroid_pre{jj,2,ii} = mtpre_resam_right;
        Mean_centroid_post{jj,1,ii} = mtpost_resam_left;
        Mean_centroid_post{jj,2,ii} = mtpost_resam_right;
    end 
end

mtpre_resam_left15 = mean(tpre_resam.*Left_cluster_pre,3,'omitnan');
mtpre_resam_right15 = mean(tpre_resam.*Right_cluster_pre,3,'omitnan');
plot3(mtpre_resam_left15(1,:),mtpre_resam_left15(2,:),mtpre_resam_left15(3,:),'k-','linewidth',3)
hold on
plot3(mtpre_resam_right15(1,:),mtpre_resam_right15(2,:),mtpre_resam_right15(3,:),'k-','linewidth',3)
for kk = 1:10000
    plot3(tpre_resam(1,:,kk),tpre_resam(2,:,kk),tpre_resam(3,:,kk),'-','color',[0.3 0.3 0.3])
    hold on
end

mpre_left = Mean_centroid_pre{jj,1};mpost_left = Mean_centroid_post{jj,1};
mpre_right = Mean_centroid_pre{jj,2};mpost_right = Mean_centroid_post{jj,2};
plot(mpre_left(1,:,1),mpre_left(3,:,1),'k-','linewidth',3)
hold on
plot(mpost_left(1,:,1),mpost_left(3,:,1),'k-','linewidth',3)
plot(mpre_right(1,:,1),mpre_right(3,:,1),'k-','linewidth',3)
plot(mpost_right(1,:,1),mpost_right(3,:,1),'k-','linewidth',3)


%% Plot Centroid path
% Coronal view
% 1.5 mm
figure(1)
for ii = 1:18
mpre_left = Mean_centroid_pre{1,1,ii};mpost_left = Mean_centroid_post{1,1,ii};
mpre_right = Mean_centroid_pre{1,2,ii};mpost_right = Mean_centroid_post{1,2,ii};
    subplot(3,6,ii)
plot(mpre_left(1,:,1),mpre_left(3,:,1),'b-','linewidth',3)
hold on
plot(mpost_left(1,:,1),mpost_left(3,:,1),'b-','linewidth',2)
plot(mpre_right(1,:,1),mpre_right(3,:,1),'r-','linewidth',2)
plot(mpost_right(1,:,1),mpost_right(3,:,1),'r-','linewidth',2)
if find(ii == 1:6)
t1 = title(append('Subject ',num2str(ii)),'Fontsize',20,'fontname','calibri');
t1.Color =  [0 0.4470 0.7410];
elseif find(ii == 7:11)
t2 = title(append('Subject ',num2str(ii-6)),'Fontsize',20,'fontname','calibri');
t2.Color =  [0.8500 0.3250 0.0980];
else
t3 = title(append('Subject ',num2str(ii-11)),'Fontsize',20,'fontname','calibri');
t3.Color =  [0.4660 0.6740 0.1880];
end
xlim([48 96]);xticks([0:24:144])
ylim([0 50]);yticks([0:10:50])
end
% 2.0 mm
figure(2)
for ii = 1:18
    subplot(3,6,ii)
plot(Left_centroid_pre20(1,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),Left_centroid_pre20(3,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),'b.-','linewidth',2)
hold on
plot(Right_centroid_pre20(1,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),Right_centroid_pre20(3,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),'b.-','linewidth',2)
plot(Left_centroid_post20(1,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),Left_centroid_post20(3,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),'r.-','linewidth',2)
plot(Right_centroid_post20(1,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),Right_centroid_post20(3,Bottomslice(2,jj,1):Bottomslice(2,jj,2),ii),'r.-','linewidth',2)
if find(ii == 1:6)
t1 = title(append('Subject ',num2str(ii)),'Fontsize',20,'fontname','calibri');
t1.Color =  [0 0.4470 0.7410];
elseif find(ii == 7:11)
t2 = title(append('Subject ',num2str(ii-6)),'Fontsize',20,'fontname','calibri');
t2.Color =  [0.8500 0.3250 0.0980];
else
t3 = title(append('Subject ',num2str(ii-11)),'Fontsize',20,'fontname','calibri');
t3.Color =  [0.4660 0.6740 0.1880];
end
xlim([36 72]);xticks([0:12:108])
ylim([0 40]);yticks([0:10:50])
end
% 2.5 mm
figure(3)
for ii = 1:18
    subplot(3,6,ii)
plot(Left_centroid_pre25(1,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),Left_centroid_pre25(3,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),'b.-','linewidth',2)
hold on
plot(Right_centroid_pre25(1,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),Right_centroid_pre25(3,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),'b.-','linewidth',2)
plot(Left_centroid_post25(1,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),Left_centroid_post25(3,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),'r.-','linewidth',2)
plot(Right_centroid_post25(1,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),Right_centroid_post25(3,Bottomslice(3,jj,1):Bottomslice(3,jj,2),ii),'r.-','linewidth',2)
if find(ii == 1:6)
t1 = title(append('Subject ',num2str(ii)),'Fontsize',20,'fontname','calibri');
t1.Color =  [0 0.4470 0.7410];
elseif find(ii == 7:11)
t2 = title(append('Subject ',num2str(ii-6)),'Fontsize',20,'fontname','calibri');
t2.Color =  [0.8500 0.3250 0.0980];
else
t3 = title(append('Subject ',num2str(ii-11)),'Fontsize',20,'fontname','calibri');
t3.Color =  [0.4660 0.6740 0.1880];
end
xlim([36 60]);xticks([0:12:96])
ylim([0 30]);yticks([0:10:50])
end

%% MDF to centroids
for ii = 1:18
    for jj = 1:3
        tpre_resam = Tracts_resam.pre{jj,ii};
        tpost_resam = Tracts_resam.post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};
        for kk = 1:10000
            d_direct = sqrt(sum((tpre_resam(:,:,kk)-mtpre_resam_left(:,:,jj)).^2));
            d_flip = sqrt(sum((flip(tpre_resam(:,:,kk),2)-mtpre_resam_left(:,:,jj)).^2));
            MDF_pre_left(kk,jj,ii) = min(min(d_direct,d_flip));
            
            d_direct = sqrt(sum((tpre_resam(:,:,kk)-mtpre_resam_right(:,:,jj)).^2));
            d_flip = sqrt(sum((flip(tpre_resam(:,:,kk),2)-mtpre_resam_right(:,:,jj)).^2));
            MDF_pre_right(kk,jj,ii) = min(min(d_direct,d_flip));
            
            d_direct = sqrt(sum((tpost_resam(:,:,kk)-mtpost_resam_left(:,:,jj)).^2));
            d_flip = sqrt(sum((flip(tpost_resam(:,:,kk),2)-mtpost_resam_left(:,:,jj)).^2));
            MDF_post_left(kk,jj,ii) = min(min(d_direct,d_flip));
            
            d_direct = sqrt(sum((tpost_resam(:,:,kk)-mtpost_resam_right(:,:,jj)).^2));
            d_flip = sqrt(sum((flip(tpost_resam(:,:,kk),2)-mtpost_resam_right(:,:,jj)).^2));
            MDF_post_right(kk,jj,ii) = min(min(d_direct,d_flip));
        end
    end
end

imagesc(MDF_pre_left(:,:,1));colorbar

subplot(2,2,1)
imagesc(squeeze(mean(MDF_pre_left,1)));colorbar
subplot(2,2,2)
imagesc(squeeze(mean(MDF_pre_right,1)));colorbar
subplot(2,2,3)
imagesc(squeeze(mean(MDF_post_left,1)));colorbar
subplot(2,2,4)
imagesc(squeeze(mean(MDF_post_right,1)));colorbar

%% Plot MDF-encoding cluster and Centroid
for ii = 1:18
    for jj = 1:3
        tpre_resam = Tracts_resam.pre{jj,ii};
        tpost_resam = Tracts_resam.post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};
for kk = 1:1000
    Left_cluster_pre(:,:,kk) = tpre_resam(1,:,kk) <= Midslice(1,1);
    Right_cluster_pre(:,:,kk) = tpre_resam(1,:,kk) > Midslice(1,1);
    Left_cluster_post(:,:,kk) = tpost_resam(1,:,kk) <= Midslice(1,1);
    Right_cluster_post(:,:,kk) = tpost_resam(1,:,kk) > Midslice(1,1);
    Left_cluster_pre = double(Left_cluster_pre);Left_cluster_pre(Left_cluster_pre == 0) = nan;
    Right_cluster_pre = double(Right_cluster_pre);Right_cluster_pre(Right_cluster_pre == 0) = nan;
    
    tpre_resam_left(:,:,kk) = tpre_resam(:,:,kk).*Left_cluster_pre(:,:,kk);
    R = MDF_pre_left(kk,1,1)/2;
    if R > 1
        color = [0 0 0 0];
    else
        color = [1 0 0] + ([0 0 1]-[1 0 0])*R;
        color = [color,0.2];
    end
    plot(tpre_resam_left(2,:,kk),tpre_resam_left(3,:,kk),'-','Color',color)
    hold on
end
plot(mpre_left(2,:,1),mpre_left(3,:,1),'k-','linewidth',3)
    end
end
colorbar

cmap = colormap(jet(10000));
cmap = cmap(randperm(length(cmap)),:);
%Set colororder and plot
ax = axes('colororder',cmap);hold on
