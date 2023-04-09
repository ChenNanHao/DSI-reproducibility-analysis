clearvars;clc;close all
%% Read files
for kk = 1:3
    if kk == 1
        cd('E:\dsi_data_7T_20200901_try\1p5mm_results_20220411')
        voxelsize = '15';
    elseif kk == 2
        cd('E:\dsi_data_7T_20200901_try\2p0mm_results_20220411')
        voxelsize = '20';
    else
        cd('E:\dsi_data_7T_20200901_try\2p5mm_results_20220411')
        voxelsize = '25';
    end       
for aa =  [1:13 15:18]
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

people = People{aa};
load(append('ODF_info_',voxelsize,'_',people),'singlefibermask','crossfibermask');
load(append('ODF_repro_',voxelsize,'_',people));
Single_SIMODF(aa,kk) = mean(ODFsim.*singlefibermask,'all','omitnan');
Crossing_SIMODF(aa,kk) = mean(ODFsim.*crossfibermask,'all','omitnan');
Single_Angledev(aa,kk) = mean(Pa1,'all','omitnan');
Crossing_Angledev(aa,kk) = mean(Pa4,'all','omitnan');

end
end
Single_SIMODF = atanh(Single_SIMODF);
Crossing_SIMODF = atanh(Crossing_SIMODF);

save('ODF results.mat','Single_SIMODF','Crossing_SIMODF','Single_Angledev','Crossing_Angledev')
%% ODFSIM
SIMODF_7TS = Single_SIMODF([1:3 5:6],:);
SIMODF_7TL = Single_SIMODF(7:11,:);
SIMODF_3T = Single_SIMODF([12:13 15:18],:);
SIMODF_mean = [mean(SIMODF_7TS)' mean(SIMODF_7TL)' mean(SIMODF_3T)'];
SIMODF_std = [std(SIMODF_7TS)' std(SIMODF_7TL)' std(SIMODF_3T)'];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(SIMODF_mean, 1);
nbars = size(SIMODF_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
figure
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, SIMODF_mean(:,i), SIMODF_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, SIMODF_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 3]);yticks([0:1:3])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
% legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff
% Group comparison
p_SIMODF = [];
c_SIMODF = [];
for ii = 1:3
x=cat(1,SIMODF_7TS(:,ii),SIMODF_7TL(:,ii),SIMODF_3T(:,ii));
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SIMODF = cat(2,p_SIMODF,p);
c_SIMODF= cat(2,c_SIMODF,c(:,6));
end
c_SIMODF = c_SIMODF';

p_SIMODF_v = [];
c_SIMODF_v = [];
for ii = 1:3
    if ii == 1
        CC = SIMODF_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = SIMODF_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = SIMODF_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SIMODF_v = cat(2,p_SIMODF_v,p);
c_SIMODF_v= cat(2,c_SIMODF_v,c(:,6));
end
c_SIMODF_v = c_SIMODF_v';


SIMODF_7TS = Crossing_SIMODF([1:3 5:6],:);
SIMODF_7TL = Crossing_SIMODF(7:11,:);
SIMODF_3T = Crossing_SIMODF([12:13 15:18],:);
SIMODF_mean = [mean(SIMODF_7TS)' mean(SIMODF_7TL)' mean(SIMODF_3T)'];
SIMODF_std = [std(SIMODF_7TS)' std(SIMODF_7TL)' std(SIMODF_3T)'];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(SIMODF_mean, 1);
nbars = size(SIMODF_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
figure
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, SIMODF_mean(:,i), SIMODF_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, SIMODF_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 3]);yticks([0:1:3])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff
% Group comparison
p_SIMODF = [];
c_SIMODF = [];
for ii = 1:3
x=cat(1,SIMODF_7TS(:,ii),SIMODF_7TL(:,ii),SIMODF_3T(:,ii));
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SIMODF = cat(2,p_SIMODF,p);
c_SIMODF= cat(2,c_SIMODF,c(:,6));
end
c_SIMODF = c_SIMODF';

p_SIMODF_v = [];
c_SIMODF_v = [];
for ii = 1:3
    if ii == 1
        CC = SIMODF_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = SIMODF_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = SIMODF_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_SIMODF_v = cat(2,p_SIMODF_v,p);
c_SIMODF_v= cat(2,c_SIMODF_v,c(:,6));
end
c_SIMODF_v = c_SIMODF_v';

%% AngleDev
AngleDev_7TS = Single_Angledev([1:3 5:6],:);
AngleDev_7TL = Single_Angledev(7:11,:);
AngleDev_3T = Single_Angledev([12:13 15:18],:);
AngleDev_mean = [mean(AngleDev_7TS)' mean(AngleDev_7TL)' mean(AngleDev_3T)'];
AngleDev_std = [std(AngleDev_7TS)' std(AngleDev_7TL)' std(AngleDev_3T)'];
X = categorical({'1.5 mm', '2.0 mm', '2.5 mm'});
X = reordercats(X,{'1.5 mm', '2.0 mm', '2.5 mm'});
ngroups = size(AngleDev_mean, 1);
nbars = size(AngleDev_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, AngleDev_mean(:,i), AngleDev_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, AngleDev_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 20]);yticks([0:5:40])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff
% Angle
p_Angle = [];
c_Angle = [];
for ii = 1:3
x=cat(1,AngleDev_7TS(:,ii),AngleDev_7TL(:,ii),AngleDev_3T(:,ii))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Angle = cat(2,p_Angle,p);
c_Angle= cat(2,c_Angle,c(:,6));
end
c_Angle = c_Angle';

p_Angle_v = [];
c_Angle_v = [];
for ii = 1:3
    if ii == 1
        CC = AngleDev_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = AngleDev_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = AngleDev_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Angle_v = cat(2,p_Angle_v,p);
c_Angle_v= cat(2,c_Angle_v,c(:,6));
end
c_Angle_v = c_Angle_v';

AngleDev_7TS = Crossing_Angledev([1:3 5:6],:);
AngleDev_7TL = Crossing_Angledev(7:11,:);
AngleDev_3T = Crossing_Angledev([12:13 15:18],:);
AngleDev_mean = [mean(AngleDev_7TS)' mean(AngleDev_7TL)' mean(AngleDev_3T)'];
AngleDev_std = [std(AngleDev_7TS)' std(AngleDev_7TL)' std(AngleDev_3T)'];
ngroups = size(AngleDev_mean, 1);
nbars = size(AngleDev_mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, AngleDev_mean(:,i), AngleDev_std(:,i), 'k.','linewidth',3)
    hold on
end
x1 = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
b = bar(x1, AngleDev_mean);hold on
set(b(1,1),'edgecolor',[0 0 0],'facecolor',[0 0 0])
set(b(1,2),'edgecolor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7])
set(b(1,3),'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
xticks([1 2 3]);ylim([0 60]);yticks([0:20:60])
set(gca,'fontname','calibri','fontsize',16,'xticklabel',X,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
legend([b(1,1) b(1,2) b(1,3)],' 7TS',' 7TL',' 3T','fontname','calibri','numcolumns',3,'fontsize',14);legend boxoff

% Angle
p_Angle = [];
c_Angle = [];
for ii = 1:3
x=cat(1,AngleDev_7TS(:,ii),AngleDev_7TL(:,ii),AngleDev_3T(:,ii))';
group=[1;1;1;1;1;2;2;2;2;2;3;3;3;3;3;3];
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Angle = cat(2,p_Angle,p);
c_Angle= cat(2,c_Angle,c(:,6));
end
c_Angle = c_Angle';

p_Angle_v = [];
c_Angle_v = [];
for ii = 1:3
    if ii == 1
        CC = AngleDev_7TS;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    elseif ii == 2
        CC = AngleDev_7TL;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    else
        CC = AngleDev_3T;
        group = reshape(repmat([1 2 3],size(CC,2),1),[size(CC,2)*3 1]);
    end  
x=cat(2,CC(1,:),CC(2,:),CC(3,:))';
[p,tbl,stats] = kruskalwallis(x',group,'off');
c = multcompare(stats);
p_Angle_v = cat(2,p_Angle_v,p);
c_Angle_v= cat(2,c_Angle_v,c(:,6));
end
c_Angle_v = c_Angle_v';

