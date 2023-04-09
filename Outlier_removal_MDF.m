clearvars;clc;close all
%% CST
cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
load('Tracts_resampling100_pre.mat');
load('Tracts_resampling100_post.mat');
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

%% Generate centroid lines
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        for kk = 1:10000
            if tpre_resam(1,:,kk) <= Midslice(jj,ii)
                left_cluster_pre(kk) = 1;
                right_cluster_pre(kk) = 0;
            else
                left_cluster_pre(kk) = 0;
                right_cluster_pre(kk) = 1;
            end
            if tpost_resam(1,:,kk) <= Midslice(jj,ii)
                left_cluster_post(kk) = 1;
                right_cluster_post(kk) = 0;
            else
                left_cluster_post(kk) = 0;
                right_cluster_post(kk) = 1;
            end
        end
        Left_cluster_pre{jj,ii} = left_cluster_pre;
        Right_cluster_pre{jj,ii} = right_cluster_pre;
        Left_cluster_post{jj,ii} = left_cluster_post;
        Right_cluster_post{jj,ii} = right_cluster_post;

        mtpre_resam_left = mean(tpre_resam(:,:,find(squeeze(Left_cluster_pre{jj,ii}))),3,'omitnan');
        mtpre_resam_right = mean(tpre_resam(:,:,find(squeeze(Right_cluster_pre{jj,ii}))),3,'omitnan');
        mtpost_resam_left = mean(tpost_resam(:,:,find(squeeze(Left_cluster_post{jj,ii}))),3,'omitnan');
        mtpost_resam_right = mean(tpost_resam(:,:,find(squeeze(Right_cluster_post{jj,ii}))),3,'omitnan');
        
        Mean_centroid_pre{jj,1,ii} = mtpre_resam_left;
        Mean_centroid_pre{jj,2,ii} = mtpre_resam_right;
        Mean_centroid_post{jj,1,ii} = mtpost_resam_left;
        Mean_centroid_post{jj,2,ii} = mtpost_resam_right;
    end 
end

% Generate MDF 
Voxel = [1.5 2.0 2.5];
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};
        for kk = 1:10000
            d_direct = sqrt(sum((tpre_resam(:,:,kk)-mtpre_resam_left).^2));
            d_flip = sqrt(sum((flip(tpre_resam(:,:,kk),2)-mtpre_resam_left).^2));
            MDF_pre_left(kk,jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);
            
            d_direct = sqrt(sum((tpre_resam(:,:,kk)-mtpre_resam_right).^2));
            d_flip = sqrt(sum((flip(tpre_resam(:,:,kk),2)-mtpre_resam_right).^2));
            MDF_pre_right(kk,jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);
            
            d_direct = sqrt(sum((tpost_resam(:,:,kk)-mtpost_resam_left).^2));
            d_flip = sqrt(sum((flip(tpost_resam(:,:,kk),2)-mtpost_resam_left).^2));
            MDF_post_left(kk,jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);
            
            d_direct = sqrt(sum((tpost_resam(:,:,kk)-mtpost_resam_right).^2));
            d_flip = sqrt(sum((flip(tpost_resam(:,:,kk),2)-mtpost_resam_right).^2));
            MDF_post_right(kk,jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);
        end
    end
    disp(append('Finished in ', People{ii}))
end

cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
save('Tracts_centroid_line_pre.mat','Mean_centroid_pre')
save('Tracts_centroid_line_post.mat','Mean_centroid_post')

%% Original Tracts 
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};

        tpre_resam_left = tpre_resam(:,:,find(squeeze(Left_cluster_pre{jj,ii})));
        figure
        for kk = 1:size(tpre_resam_left,3)
            plot(tpre_resam_left(2,:,kk),tpre_resam_left(3,:,kk),'Color',[0.7 0 0.3],'linewidth',0.002)
            hold on
        end
            plot(mtpre_resam_left(2,:),mtpre_resam_left(3,:),'y','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([0 60])
            elseif jj == 2
                xlim([24 84]);ylim([0 46])
            else
                xlim([24 72]);ylim([0 36])
            end
            
        tpre_resam_right = tpre_resam(:,:,find(squeeze(Right_cluster_pre{jj,ii})));
        for kk = 1:size(tpre_resam_right,3)
            plot(tpre_resam_right(2,:,kk),tpre_resam_right(3,:,kk),'Color',[0.3 0 0.7],'linewidth',0.002)
            hold on
        end
            plot(mtpre_resam_right(2,:),mtpre_resam_right,'y','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([0 60])
            elseif jj == 2
                xlim([24 84]);ylim([0 46])
            else
                xlim([24 72]);ylim([0 36])
            end
    end
end
%% Outlier removal
Outlier_threshold = 2.5;
cmap = cat(1,linspace(0,1,1000),zeros(1,1000),linspace(1,0,1000));
MDFcmap = linspace(0,10,1000);
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};

        idx = find(MDF_pre_left(find(squeeze(Left_cluster_pre{jj,ii})),jj,ii) < Outlier_threshold);
        tpre_resam_left = tpre_resam(:,:,find(squeeze(Left_cluster_pre{jj,ii})));
        mtpre_resam_left_outremove = mean(tpre_resam_left(:,:,idx),3,'omitnan');
        figure
        for kk = 1:size(idx,1)
            [~,color_idx] = min(abs(MDF_pre_left(idx(kk),jj,ii)-MDFcmap));
            plot(tpre_resam_left(2,:,idx(kk)),tpre_resam_left(3,:,idx(kk)),'Color',cmap(:,color_idx),'linewidth',0.02)
            hold on
        end
            plot(mtpre_resam_left_outremove(2,:),mtpre_resam_left_outremove(3,:),'y','linewidth',2)
            hold on
            plot(mtpre_resam_left(2,:),mtpre_resam_left(3,:),'g','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([35 50])
            elseif jj == 2
                xlim([24 84]);ylim([26.83 38.33])
            else
                xlim([24 72]);ylim([21 30])
            end
            
        idx = find(MDF_pre_right(find(squeeze(Right_cluster_pre(ii,jj,:))),jj,ii) < Outlier_threshold);
        tpre_resam_right = tpre_resam(:,:,find(squeeze(Right_cluster_pre(ii,jj,:))));
        mtpre_resam_left_outremove(:,:,jj) = mean(tpre_resam(:,:,idx),3,'omitnan');
        figure
        for kk = 1:size(idx,1)
            [~,color_idx] = min(abs(MDF_pre_left(idx(kk),jj,ii)-MDFcmap));
            plot(tpre_resam(1,:,idx(kk)),tpre_resam(2,:,idx(kk)),'Color',cmap(:,color_idx),'linewidth',0.02)
            hold on
        end
            plot(mtpre_resam_left_outremove(1,:,jj),mtpre_resam_left_outremove(2,:,jj),'y','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([35 50])
            elseif jj == 2
                xlim([24 84]);ylim([26.83 38.33])
            else
                xlim([24 72]);ylim([21 30])
            end
    end
end

% Grey tract bundles
% Pre
Outlier_threshold = 2.5;
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};

        idx = find(MDF_pre_left(find(squeeze(Left_cluster_pre{jj,ii})),jj,ii) < Outlier_threshold);
        tpre_resam_left = tpre_resam(:,:,find(squeeze(Left_cluster_pre{jj,ii})));
        mtpre_resam_left_outremove = mean(tpre_resam_left(:,:,idx),3,'omitnan');
        figure
        for kk = 1:size(idx,1)
            [~,color_idx] = min(abs(MDF_pre_left(idx(kk),jj,ii)-MDFcmap));
            plot(tpre_resam_left(2,:,idx(kk)),tpre_resam_left(3,:,idx(kk)),'Color',[0.7 0 0.3],'linewidth',0.002)
            hold on
        end
            plot(mtpre_resam_left_outremove(2,:),mtpre_resam_left_outremove(3,:),'y','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([0 60])
            elseif jj == 2
                xlim([24 84]);ylim([0 46])
            else
                xlim([24 72]);ylim([0 36])
            end
            
        idx = find(MDF_pre_right(find(squeeze(Right_cluster_pre{jj,ii})),jj,ii) < Outlier_threshold);
        tpre_resam_right = tpre_resam(:,:,find(squeeze(Right_cluster_pre{jj,ii})));
        mtpre_resam_right_outremove = mean(tpre_resam_right(:,:,idx),3,'omitnan');
        for kk = 1:size(idx,1)
            [~,color_idx] = min(abs(MDF_pre_right(idx(kk),jj,ii)-MDFcmap));
            plot(tpre_resam_right(2,:,idx(kk)),tpre_resam_right(3,:,idx(kk)),'Color',[0.3 0 0.7],'linewidth',0.002)
            hold on
        end
            plot(mtpre_resam_right_outremove(2,:),mtpre_resam_right_outremove(3,:),'y','linewidth',2)
            if jj == 1
                xlim([32 112]);ylim([0 60])
            elseif jj == 2
                xlim([24 84]);ylim([0 46])
            else
                xlim([24 72]);ylim([0 36])
            end
    end
end
% Post
for ii = [1:13 15:18]
    for jj = 1:3
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};

        idx = find(MDF_post_left(find(squeeze(Left_cluster_post{jj,ii})),jj,ii) < Outlier_threshold);
        tpost_resam_left = tpost_resam(:,:,find(squeeze(Left_cluster_post{jj,ii})));
        mtpost_resam_left_outremove = mean(tpost_resam_left(:,:,idx),3,'omitnan');
%         figure
%         for kk = 1:size(idx,1)
%             [~,color_idx] = min(abs(MDF_post_left(idx(kk),jj,ii)-MDFcmap));
%             plot(tpost_resam_left(1,:,idx(kk)),tpost_resam_left(2,:,idx(kk)),'Color',[0.7 0 0.3],'linewidth',0.002)
%             hold on
%         end
%             plot(mtpost_resam_left_outremove(1,:),mtpost_resam_left_outremove(2,:),'y','linewidth',2)
%             if jj == 1
%                 xlim([52 92]);ylim([32 112])
%             elseif jj == 2
%                 xlim([39 69]);ylim([24 84])
%             else
%                 xlim([36 60]);ylim([24 72])
%             end
            
        idx = find(MDF_post_right(find(squeeze(Right_cluster_post{jj,ii})),jj,ii) < Outlier_threshold);
        tpost_resam_right = tpost_resam(:,:,find(squeeze(Right_cluster_post{jj,ii})));
        mtpost_resam_right_outremove = mean(tpost_resam_right(:,:,idx),3,'omitnan');
%         for kk = 1:size(idx,1)
%             [~,color_idx] = min(abs(MDF_post_right(idx(kk),jj,ii)-MDFcmap));
%             plot(tpost_resam_right(1,:,idx(kk)),tpost_resam_right(2,:,idx(kk)),'Color',[0.3 0 0.7],'linewidth',0.002)
%             hold on
%         end
%             plot(mtpost_resam_right_outremove(1,:),mtpost_resam_right_outremove(2,:),'y','linewidth',2)
%             if jj == 1
%                 xlim([52 92]);ylim([32 112])
%             elseif jj == 2
%                 xlim([39 69]);ylim([24 84])
%             else
%                 xlim([36 60]);ylim([24 72])
%             end
    end
end

%% Outlier removal (calculations)
% Pre & Post mean centroid calculation
Outlier_threshold = 2.5;
idx_pre = 0;idx_post = 0;
for jj = 1:3 
    figure(jj)
    for ii = [1:13 15:18]       
        tpre_resam = Tracts_resam_pre{jj,ii};
        tpost_resam = Tracts_resam_post{jj,ii};
        mtpre_resam_left = Mean_centroid_pre{jj,1,ii};
        mtpre_resam_right = Mean_centroid_pre{jj,2,ii};
        mtpost_resam_left = Mean_centroid_post{jj,1,ii};
        mtpost_resam_right = Mean_centroid_post{jj,2,ii};

        idx = find(MDF_pre_left(find(squeeze(Left_cluster_pre{jj,ii})),jj,ii) < Outlier_threshold);
        tpre_resam_left = tpre_resam(:,:,find(squeeze(Left_cluster_pre{jj,ii})));
        mtpre_resam_left_outremove = mean(tpre_resam_left(:,:,idx),3,'omitnan');
        Tracts_resam_pre_outremove{jj,1,ii} = tpre_resam_left(:,:,idx);
        idx_pre = idx_pre+size(idx,1);
        
        idx = find(MDF_pre_right(find(squeeze(Right_cluster_pre{jj,ii})),jj,ii) < Outlier_threshold);
        tpre_resam_right = tpre_resam(:,:,find(squeeze(Right_cluster_pre{jj,ii})));
        mtpre_resam_right_outremove = mean(tpre_resam_right(:,:,idx),3,'omitnan');
        Tracts_resam_pre_outremove{jj,2,ii} = tpre_resam_right(:,:,idx);
        idx_pre = idx_pre+size(idx,1);

        idx = find(MDF_post_left(find(squeeze(Left_cluster_post{jj,ii})),jj,ii) < Outlier_threshold);
        tpost_resam_left = tpost_resam(:,:,find(squeeze(Left_cluster_post{jj,ii})));
        mtpost_resam_left_outremove = mean(tpost_resam_left(:,:,idx),3,'omitnan');
        Tracts_resam_post_outremove{jj,1,ii} = tpost_resam_left(:,:,idx);
        idx_post = idx_post+size(idx,1);

        idx = find(MDF_post_right(find(squeeze(Right_cluster_post{jj,ii})),jj,ii) < Outlier_threshold);
        tpost_resam_right = tpost_resam(:,:,find(squeeze(Right_cluster_post{jj,ii})));
        mtpost_resam_right_outremove = mean(tpost_resam_right(:,:,idx),3,'omitnan');
        Tracts_resam_post_outremove{jj,2,ii} = tpost_resam_right(:,:,idx);
        idx_post = idx_post+size(idx,1);

        if find(ii == [12 13])
            subplot(3,6,ii+1)
        else
            subplot(3,6,ii)
        end      
        
        plot(mtpre_resam_left_outremove(1,:),mtpre_resam_left_outremove(3,:),'b','linewidth',2)
        hold on
        plot(mtpre_resam_right_outremove(1,:),mtpre_resam_right_outremove(3,:),'b','linewidth',2)
        plot(mtpost_resam_left_outremove(1,:),mtpost_resam_left_outremove(3,:),'r','linewidth',2)
        plot(mtpost_resam_right_outremove(1,:),mtpost_resam_right_outremove(3,:),'r','linewidth',2)

        
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
        
            if jj == 1
                xlim([32 112]);ylim([0 60])
            elseif jj == 2
                xlim([24 84]);ylim([0 46])
            else
                xlim([24 72]);ylim([0 36])
            end        
        Mean_centroid_pre_outremove{jj,1,ii} = mtpre_resam_left_outremove;
        Mean_centroid_pre_outremove{jj,2,ii} = mtpre_resam_right_outremove;
        Mean_centroid_post_outremove{jj,1,ii} = mtpost_resam_left_outremove;
        Mean_centroid_post_outremove{jj,2,ii} = mtpost_resam_right_outremove;
        
        IDX(jj,1,ii) = idx_pre;
        IDX(jj,2,ii) = idx_post;
        idx_pre = 0;idx_post = 0;
    end
end

cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
save('Tracts_resampling100_pre_outremove.mat','Tracts_resam_pre_outremove','IDX')
save('Tracts_resampling100_post_outremove.mat','Tracts_resam_post_outremove','IDX')
save('Tracts_centroid_line_pre_outremove.mat','Mean_centroid_pre_outremove')
save('Tracts_centroid_line_post_outremove.mat','Mean_centroid_post_outremove')

%% Centroid line quantification
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
% Distance (Chamfer)
for jj = 1:3
for ii = [1:13 15:18]
    mpre_left = Mean_centroid_pre_outremove{jj,1,ii};
    mpre_right = Mean_centroid_pre_outremove{jj,2,ii}; 
    mpost_left = Mean_centroid_post_outremove{jj,1,ii};
    mpost_right = Mean_centroid_post_outremove{jj,2,ii}; 

    d_left(jj,ii) = mean(sqrt(sum((mpre_left-mpost_left).^2)))*Voxel(jj);
    d_right(jj,ii) = mean(sqrt(sum((mpre_right-mpost_right).^2)))*Voxel(jj);
end
end
Centroidline_chamferDIST = (d_left+d_right)*0.5;
Center_chamferDIST_7TS = mean(Centroidline_chamferDIST(:,1:6),3);
Center_chamferDIST_7TL = mean(Centroidline_chamferDIST(:,7:11),3);
Center_chamferDIST_3T = mean(Centroidline_chamferDIST(:,[12:13 15:18]),3);
Center_chamferDIST = [mean(Center_chamferDIST_7TS,2) mean(Center_chamferDIST_7TL,2) mean(Center_chamferDIST_3T,2)];
Center_chamferDIST_std = [std(Center_chamferDIST_7TS,[],2) std(Center_chamferDIST_7TL,[],2) std(Center_chamferDIST_3T,[],2)];

ngroups = size(Center_chamferDIST, 1);
nbars = size(Center_chamferDIST, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, Center_chamferDIST(:,i), Center_chamferDIST_std(:,i), '.','color',[0 0 0],'linewidth',2)
    elseif i == 2
        errorbar(x, Center_chamferDIST(:,i), Center_chamferDIST_std(:,i), '.','color',[0.3 0.3 0.3],'linewidth',2)
    else
        errorbar(x, Center_chamferDIST(:,i), Center_chamferDIST_std(:,i), '.','color',[0.6 0.6 0.6],'linewidth',2)
    end
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, Center_chamferDIST);
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.3 0.3 0.3],'facecolor',[0.3 0.3 0.3])
set(b(1,3),'edgecolor',[0.6 0.6 0.6],'facecolor',[0.6 0.6 0.6])
xticks([1 2 3]);ylim([0 6]);yticks([0:2:6])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Chamfer distance (mm)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% Distance (MDF)
for jj = 1:3
for ii = [1:13 15:18]
    mpre_left = Mean_centroid_pre_outremove{jj,1,ii};
    mpre_right = Mean_centroid_pre_outremove{jj,2,ii}; 
    mpost_left = Mean_centroid_post_outremove{jj,1,ii};
    mpost_right = Mean_centroid_post_outremove{jj,2,ii}; 

    d_direct = sqrt(sum((mpre_left-mpost_left).^2));
    d_flip = sqrt(sum((flip(mpre_left,2)-mpost_left).^2));
    d_left(jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);

    d_direct = sqrt(sum((mpre_right-mpost_right).^2));
    d_flip = sqrt(sum((flip(mpre_right,2)-mpost_right).^2));
    d_right(jj,ii) = min(min(d_direct,d_flip))*Voxel(jj);
end
end
Centroidline_MDF = (d_left+d_right)*0.5;
Center_MDF_7TS = mean(Centroidline_MDF(:,1:6),3);
Center_MDF_7TL = mean(Centroidline_MDF(:,7:11),3);
Center_MDF_3T = mean(Centroidline_MDF(:,[12:13 15:18]),3);
Center_MDF = [mean(Center_MDF_7TS,2) mean(Center_MDF_7TL,2) mean(Center_MDF_3T,2)];
Center_MDF_std = [std(Center_MDF_7TS,[],2) std(Center_MDF_7TL,[],2) std(Center_MDF_3T,[],2)];

ngroups = size(Center_MDF, 1);
nbars = size(Center_MDF, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, Center_MDF(:,i), Center_MDF_std(:,i), '.','color',[0 0 0],'linewidth',2)
    elseif i == 2
        errorbar(x, Center_MDF(:,i), Center_MDF_std(:,i), '.','color',[0.3 0.3 0.3],'linewidth',2)
    else
        errorbar(x, Center_MDF(:,i), Center_MDF_std(:,i), '.','color',[0.6 0.6 0.6],'linewidth',2)
    end
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, Center_MDF);
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.3 0.3 0.3],'facecolor',[0.3 0.3 0.3])
set(b(1,3),'edgecolor',[0.6 0.6 0.6],'facecolor',[0.6 0.6 0.6])
xticks([1 2 3]);ylim([0 6]);yticks([0:2:6])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Chamfer distance (mm)','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% Vector Inner product
for jj = 1:3
for ii = [1:13 15:18]
    mpre_left = Mean_centroid_pre_outremove{jj,1,ii};
    mpre_right = Mean_centroid_pre_outremove{jj,2,ii}; 
    mpost_left = Mean_centroid_post_outremove{jj,1,ii};
    mpost_right = Mean_centroid_post_outremove{jj,2,ii}; 
    % Inner product of vector (Left side)
    l = diff(mpre_left,1,2);
    l1 = diff(mpost_left,1,2);
    if jj == 1
        slice_interval = 40;
    elseif jj == 2
        slice_interval = 30;
    else
        slice_interval = 24;
    end
        
    for aa = 1:slice_interval-1
        inner(aa) = dot(l(:,aa),l1(:,aa))/(norm(l(:,aa))*norm(l1(:,aa)));
    end
    Vector_inner_left(jj,ii,1) = mean(abs(inner),'omitnan');
    % Inner product of vector (Left side)
    l = diff(mpre_right,1,2);
    l1 = diff(mpost_right,1,2);
    for aa = 1:slice_interval-1
    inner(aa) = dot(l(:,aa),l1(:,aa))/(norm(l(:,aa))*norm(l1(:,aa)));
    end
    Vector_inner_right(jj,ii,1) = mean(abs(inner),'omitnan');
end
end
Vector_inner_Centroidline = (Vector_inner_left+Vector_inner_right)*0.5;
Center_Vinner_7TS = mean(Vector_inner_Centroidline(:,1:6),3);
Center_Vinner_7TL = mean(Vector_inner_Centroidline(:,7:11),3);
Center_Vinner_3T = mean(Vector_inner_Centroidline(:,[12:13 15:18]),3);
Center_Vinner = [mean(Center_Vinner_7TS,2) mean(Center_Vinner_7TL,2) mean(Center_Vinner_3T,2)];
Center_Vinner_std = [std(Center_Vinner_7TS,[],2) std(Center_Vinner_7TL,[],2) std(Center_Vinner_3T,[],2)];

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, Center_Vinner(:,i), Center_Vinner_std(:,i), '.','color',[0 0 0],'linewidth',2)
    elseif i == 2
        errorbar(x, Center_Vinner(:,i), Center_Vinner_std(:,i), '.','color',[0.3 0.3 0.3],'linewidth',2)
    else
        errorbar(x, Center_Vinner(:,i), Center_Vinner_std(:,i), '.','color',[0.6 0.6 0.6],'linewidth',2)
    end
    hold on
end
b = bar(x1, Center_Vinner);
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.3 0.3 0.3],'facecolor',[0.3 0.3 0.3])
set(b(1,3),'edgecolor',[0.6 0.6 0.6],'facecolor',[0.6 0.6 0.6])
xticks([1 2 3]);ylim([0.95 1.01]);yticks([0.95:0.02:1.01])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Vector inner (a.u.)','fontname','calibri','fontsize',18)
% legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% Length Ratio
for jj = 1:3
for ii = [1:13 15:18]
    voxelsize = [1.5 2.0 2.5];
    mpre_left = Mean_centroid_pre_outremove{jj,1,ii};
    mpre_right = Mean_centroid_pre_outremove{jj,2,ii}; 
    mpost_left = Mean_centroid_post_outremove{jj,1,ii};
    mpost_right = Mean_centroid_post_outremove{jj,2,ii}; 
    l = diff(mpre_left,1,2);
    l1 = diff(mpost_left,1,2);
    Length_pre_left(ii,jj) = sum(sqrt(sum(l.^2)))*voxelsize(jj);
    Length_post_left(ii,jj) = sum(sqrt(sum(l1.^2)))*voxelsize(jj);
    
    l = diff(mpre_right,1,2);
    l1 = diff(mpost_right,1,2);
    Length_pre_right(ii,jj) = sum(sqrt(sum(l.^2)))*voxelsize(jj);
    Length_post_right(ii,jj) = sum(sqrt(sum(l1.^2)))*voxelsize(jj);
end
end
Length_ratio(:,:,1) = abs(1-(Length_pre_left./Length_post_left));
Length_ratio(:,:,2) = abs(1-(Length_pre_right./Length_post_right));
Center_Lengthratio_7TS = mean(Length_ratio(1:6,:,:),3)';
Center_Lengthratio_7TL = mean(Length_ratio(7:11,:,:),3)';
Center_Lengthratio_3T = mean(Length_ratio([12:13 15:18],:,:),3)';
Center_Lengthratio = [mean(Center_Lengthratio_7TS,2) mean(Center_Lengthratio_7TL,2) mean(Center_Lengthratio_3T,2)];
Center_Lengthratio_std = [std(Center_Lengthratio_7TS,[],2) std(Center_Lengthratio_7TL,[],2) std(Center_Lengthratio_3T,[],2)];

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, Center_Lengthratio(:,i), Center_Lengthratio_std(:,i), '.','color',[0 0 0],'linewidth',2)
    elseif i == 2
        errorbar(x, Center_Lengthratio(:,i), Center_Lengthratio_std(:,i), '.','color',[0.3 0.3 0.3],'linewidth',2)
    else
        errorbar(x, Center_Lengthratio(:,i), Center_Lengthratio_std(:,i), '.','color',[0.6 0.6 0.6],'linewidth',2)
    end
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, Center_Lengthratio);
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.3 0.3 0.3],'facecolor',[0.3 0.3 0.3])
set(b(1,3),'edgecolor',[0.6 0.6 0.6],'facecolor',[0.6 0.6 0.6])
xticks([1 2 3]);ylim([0 0.04]);yticks([0:0.01:0.6])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Length ratio (a.u.)','fontname','calibri','fontsize',18)
% legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

cd 'E:\dsi_data_7T_20200901_try\Tracts_results_20220210'
save('Tracts_centroid_line_results.mat','Center_chamferDIST_7TS','Center_chamferDIST_7TL','Center_chamferDIST_3T',...
    'Center_Vinner_7TS','Center_Vinner_7TL','Center_Vinner_3T','Center_Lengthratio_7TS','Center_Lengthratio_7TL','Center_Lengthratio_3T')

%% H-test with posthoc test (different datasets)
% Chamfer distance
p_chamferDIST = [];
c_chamferDIST = [];
for ii = 1:3
x=cat(2,Center_chamferDIST_7TS(ii,:),Center_chamferDIST_7TL(ii,:),Center_chamferDIST_3T(ii,:))';
group=[1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_chamferDIST = cat(2,p_chamferDIST,p);
c_chamferDIST= cat(2,c_chamferDIST,c(:,6));
end
c_chamferDIST = c_chamferDIST';

% Vector inner
p_Vinner = [];
c_Vinner = [];
for ii = 1:3
x=cat(2,Center_Vinner_7TS(ii,:),Center_Vinner_7TL(ii,:),Center_Vinner_3T(ii,:))';
group=[1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Vinner=cat(2,p_Vinner,p);
c_Vinner= cat(2,c_Vinner,c(:,6));
end
c_Vinner = c_Vinner';

% Length ratio
p_Lengthratio = [];
c_Lengthratio = [];
for ii = 1:3
x=cat(2,Center_Lengthratio_7TS(ii,:),Center_Lengthratio_7TL(ii,:),Center_Lengthratio_3T(ii,:))';
group=[1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Lengthratio=cat(2,p_Lengthratio,p);
c_Lengthratio= cat(2,c_Lengthratio,c(:,6));
end
c_Lengthratio = c_Lengthratio';

%% H-test with posthoc test (different voxel sizes)
% Chamfer distance
p_chamferDIST_v = [];
c_chamferDIST_v = [];
for ii = 1:3
    if ii == 1
        CC = Center_chamferDIST_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = Center_chamferDIST_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = Center_chamferDIST_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_chamferDIST_v = cat(2,p_chamferDIST_v,p);
c_chamferDIST_v= cat(2,c_chamferDIST_v,c(:,6));
end
c_chamferDIST_v = c_chamferDIST_v';

% Vector inner
p_Vinner_v = [];
c_Vinner_v = [];
for ii = 1:3
    if ii == 1
        CC = Center_Vinner_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = Center_Vinner_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = Center_Vinner_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Vinner_v=cat(2,p_Vinner_v,p);
c_Vinner_v= cat(2,c_Vinner_v,c(:,6));
end
c_Vinner_v = c_Vinner_v';

% Length ratio
p_Lengthratio_v = [];
c_Lengthratio_v = [];
for ii = 1:3
    if ii == 1
        CC = Center_Lengthratio_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = Center_Lengthratio_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = Center_Lengthratio_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Lengthratio_v=cat(2,p_Lengthratio_v,p);
c_Lengthratio_v= cat(2,c_Lengthratio_v,c(:,6));
end
c_Lengthratio_v = c_Lengthratio_v';

