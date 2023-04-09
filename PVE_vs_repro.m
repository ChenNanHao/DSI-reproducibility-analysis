clearvars;clc;close all
%%
cd('E:\dsi_data_7T_20200901_try')
load('ODF results')
load('wDICE CST results')
load('SNR_results_20220527.mat')
cd('PVE')
load('PVE results')

PVE_7TS = PVE_single([1:3 5:6],:);
PVE_7TL = PVE_single(7:11,:);
PVE_3T = PVE_single([12:13 15:18],:);
SIMODF_7TS = Single_SIMODF([1:3 5:6],:);
SIMODF_7TL = Single_SIMODF(7:11,:);
SIMODF_3T = Single_SIMODF([12:13 15:18],:);
AngleDev_7TS = Single_Angledev([1:3 5:6],:);
AngleDev_7TL = Single_Angledev(7:11,:);
AngleDev_3T = Single_Angledev([12:13 15:18],:);
wDICE_7TS = mean(wDICE(:,[1:3 5:6]),3);
wDICE_7TL = mean(wDICE(:,7:11),3);
wDICE_3T = mean(wDICE(:,[12:13 15:18]),3);
wDICE_homo_7TS = mean(wDICE_homo(:,[1:3 5:6]),3);
wDICE_homo_7TL = mean(wDICE_homo(:,7:11),3);
wDICE_homo_3T = mean(wDICE_homo(:,[12:13 15:18]),3);
wDICE_hetero_7TS = mean(wDICE_hetero(:,[1:3 5:6]),3);
wDICE_hetero_7TL = mean(wDICE_hetero(:,7:11),3);
wDICE_hetero_3T = mean(wDICE_hetero(:,[12:13 15:18]),3);

for kk = 1:3
    [r,p] = corrcoef(PVE_7TS(:,kk),SIMODF_7TS(:,kk));
    R(1,3*(kk-1)+1) = r(2,1); P(1,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),SIMODF_7TL(:,kk));
    R(1,3*(kk-1)+2) = r(2,1); P(1,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),SIMODF_3T(:,kk));
    R(1,3*(kk-1)+3) = r(2,1); P(1,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),AngleDev_7TS(:,kk));
    R(2,3*(kk-1)+1) = r(2,1); P(2,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),AngleDev_7TL(:,kk));
    R(2,3*(kk-1)+2) = r(2,1); P(2,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),AngleDev_3T(:,kk));
    R(2,3*(kk-1)+3) = r(2,1); P(2,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_7TS(kk,:));
    R(3,3*(kk-1)+1) = r(2,1); P(3,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_7TL(kk,:));
    R(3,3*(kk-1)+2) = r(2,1); P(3,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_3T(kk,:));
    R(3,3*(kk-1)+3) = r(2,1); P(3,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_homo_7TS(kk,:));
    R(4,3*(kk-1)+1) = r(2,1); P(4,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_homo_7TL(kk,:));
    R(4,3*(kk-1)+2) = r(2,1); P(4,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_homo_3T(kk,:));
    R(4,3*(kk-1)+3) = r(2,1); P(4,3*(kk-1)+3) = p(2,1); 

    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_hetero_7TS(kk,:));
    R(5,3*(kk-1)+1) = r(2,1); P(5,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_hetero_7TL(kk,:));
    R(5,3*(kk-1)+2) = r(2,1); P(5,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_hetero_3T(kk,:));
    R(5,3*(kk-1)+3) = r(2,1); P(5,3*(kk-1)+3) = p(2,1); 
end
imagesc(R,[-1 1]);colormap jet;set(gca,'xtick',[],'ytick',[])

PVE_7TS = PVE_crossing([1:3 5:6],:);
PVE_7TL = PVE_crossing(7:11,:);
PVE_3T = PVE_crossing([12:13 15:18],:);
SIMODF_7TS = Crossing_SIMODF([1:3 5:6],:);
SIMODF_7TL = Crossing_SIMODF(7:11,:);
SIMODF_3T = Crossing_SIMODF([12:13 15:18],:);
AngleDev_7TS = Crossing_Angledev([1:3 5:6],:);
AngleDev_7TL = Crossing_Angledev(7:11,:);
AngleDev_3T = Crossing_Angledev([12:13 15:18],:);

for kk = 1:3
    [r,p] = corrcoef(PVE_7TS(:,kk),SIMODF_7TS(:,kk));
    R(1,3*(kk-1)+1) = r(2,1); P(1,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),SIMODF_7TL(:,kk));
    R(1,3*(kk-1)+2) = r(2,1); P(1,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),SIMODF_3T(:,kk));
    R(1,3*(kk-1)+3) = r(2,1); P(1,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),AngleDev_7TS(:,kk));
    R(2,3*(kk-1)+1) = r(2,1); P(2,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),AngleDev_7TL(:,kk));
    R(2,3*(kk-1)+2) = r(2,1); P(2,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),AngleDev_3T(:,kk));
    R(2,3*(kk-1)+3) = r(2,1); P(2,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_7TS(kk,:));
    R(3,3*(kk-1)+1) = r(2,1); P(3,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_7TL(kk,:));
    R(3,3*(kk-1)+2) = r(2,1); P(3,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_3T(kk,:));
    R(3,3*(kk-1)+3) = r(2,1); P(3,3*(kk-1)+3) = p(2,1); 
    
    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_homo_7TS(kk,:));
    R(4,3*(kk-1)+1) = r(2,1); P(4,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_homo_7TL(kk,:));
    R(4,3*(kk-1)+2) = r(2,1); P(4,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_homo_3T(kk,:));
    R(4,3*(kk-1)+3) = r(2,1); P(4,3*(kk-1)+3) = p(2,1); 

    [r,p] = corrcoef(PVE_7TS(:,kk),wDICE_hetero_7TS(kk,:));
    R(5,3*(kk-1)+1) = r(2,1); P(5,3*(kk-1)+1) = p(2,1); 
    [r,p] = corrcoef(PVE_7TL(:,kk),wDICE_hetero_7TL(kk,:));
    R(5,3*(kk-1)+2) = r(2,1); P(5,3*(kk-1)+2) = p(2,1); 
    [r,p] = corrcoef(PVE_3T(:,kk),wDICE_hetero_3T(kk,:));
    R(5,3*(kk-1)+3) = r(2,1); P(5,3*(kk-1)+3) = p(2,1); 
end
imagesc(R,[-1 1]);colormap jet;set(gca,'xtick',[],'ytick',[])

%%
for kk = 1:3
    [r,p] = corrcoef(PVE_crossing([1:3 5:13 15:18],kk),Crossing_SIMODF([1:3 5:13 15:18],kk));
    R(kk) = r(2,1); P(kk) = p(2,1); 
end

for kk = 1:3
    [r,p] = corrcoef(PVE_single([1:3 5:13 15:18],kk),Single_Angledev([1:3 5:13 15:18],kk));
    R(kk) = r(2,1); P(kk) = p(2,1); 
end

for kk = 1:3
    [r,p] = corrcoef(PVE_crossing([1:3 5:13 15:18],kk),Crossing_Angledev([1:3 5:13 15:18],kk));
    R(kk) = r(2,1); P(kk) = p(2,1); 
end

wDICE_7TS = mean(wDICE(:,[1:3 5:6]),3);
wDICE_7TL = mean(wDICE(:,7:11),3);
wDICE_3T = mean(wDICE(:,[12:13 15:18]),3);

for kk = 1:3
    [r,p] = corrcoef(PVE_single([1:3 5:13 15:18],kk),wDICE_hetero(kk,[1:3 5:13 15:18]));
    R(kk) = r(2,1); P(kk) = p(2,1); 
end

for kk = 1:3
    [r,p] = corrcoef(PVE_crossing([1:3 5:13 15:18],kk),wDICE_hetero(kk,[1:3 5:13 15:18]));
    R(kk) = r(2,1); P(kk) = p(2,1); 
end

%% Linear regression (single-fiber)
% PVE vs. repro
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for kk = 1:3
plot(PVE_single([1:3 5:13 15:18],kk),Single_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_single([1:3 5:13 15:18],kk)),max(PVE_single([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],kk),Single_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 15]);yticks([0:5:40]);xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end
x = linspace(min(PVE_single([1:3 5:13 15:18],:),[],'all'),max(PVE_single([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],:),Single_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(PVE_single([1:3 5:13 15:18],kk),Single_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_single([1:3 5:13 15:18],kk)),max(PVE_single([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],kk),Single_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([1 3]);yticks([1:0.5:3]);xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
end
x = linspace(min(PVE_single([1:3 5:13 15:18],:),[],'all'),max(PVE_single([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],:),Single_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

% SNR vs. repro
for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),Single_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Single_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 15]);yticks([0:5:40]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18])',Single_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),Single_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Single_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([1 3]);yticks([1:0.5:3]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18])',Single_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

%% Linear regression (crossing-fiber)
% PVE vs. repro
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for kk = 1:3
plot(PVE_crossing([1:3 5:13 15:18],kk),Crossing_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:13 15:18],kk)),max(PVE_crossing([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],kk),Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([15 40]);yticks([15:5:40]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end
x = linspace(min(PVE_crossing([1:3 5:13 15:18],:),[],'all'),max(PVE_crossing([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],:),Crossing_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(PVE_crossing([1:3 5:13 15:18],kk),Crossing_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:13 15:18],kk)),max(PVE_crossing([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],kk),Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 2]);yticks([0:0.5:2]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
end
x = linspace(min(PVE_crossing([1:3 5:13 15:18],:),[],'all'),max(PVE_crossing([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],:),Crossing_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

% SNR vs. repro
for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),Crossing_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([15 40]);yticks([15:5:40]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18])',Crossing_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(SNR(kk,[1:3 5:13 15:18]),Crossing_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 2]);yticks([0:0.5:2]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' SNR (a.u.)','fontname','calibri','fontsize',18)
end
x = linspace(min(SNR(:,[1:3 5:13 15:18]),[],'all'),max(SNR(:,[1:3 5:13 15:18]),[],'all'),100);
[P,S] = polyfit(SNR(:,[1:3 5:13 15:18])',Crossing_SIMODF([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

%% Linear regression (Tract-based repro)
% PVE vs. repro
resol_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for kk = 1:3
plot(PVE_single([1:3 5:13 15:18],kk),Single_Angledev([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_single([1:3 5:13 15:18],kk)),max(PVE_single([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],kk),Single_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([0 15]);yticks([0:5:40]);xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end
x = linspace(min(PVE_single([1:3 5:13 15:18],:),[],'all'),max(PVE_single([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],:),Single_Angledev([1:3 5:13 15:18],:),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on

for kk = 1:3
plot(PVE_single([1:3 5:13 15:18],kk),Single_SIMODF([1:3 5:13 15:18],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_single([1:3 5:13 15:18],kk)),max(PVE_single([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],kk),Single_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
ylim([1 3]);yticks([1:0.5:3]);xlim([5 20]);xticks([0:5:20])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
end
x = linspace(min(PVE_single([1:3 5:13 15:18],:),[],'all'),max(PVE_single([1:3 5:13 15:18],:),[],'all'),100);
[P,S] = polyfit(PVE_single([1:3 5:13 15:18],:),Single_SIMODF([1:3 5:13 15:18],:),1);
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

%% Linear regression (separate different datasets
% PVE vs. repro
for kk = 1:3
    subplot(2,3,kk)
plot(PVE_crossing([1:3 5:6],kk),Crossing_SIMODF([1:3 5:6],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:6],kk)),max(PVE_crossing([1:3 5:6],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:6],kk),Crossing_SIMODF([1:3 5:6],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(PVE_crossing([7:11],kk),Crossing_SIMODF([7:11],kk),'*','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([7:11],kk)),max(PVE_crossing([7:11],kk)),100);
[P,S] = polyfit(PVE_crossing([7:11],kk),Crossing_SIMODF([7:11],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(PVE_crossing([12:13 15:18],kk),Crossing_SIMODF([12:13 15:18],kk),'s','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([12:13 15:18],kk)),max(PVE_crossing([12:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([12:13 15:18],kk),Crossing_SIMODF([12:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:13 15:18],kk)),max(PVE_crossing([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],kk),Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on
ylim([0 2]);yticks([0:0.5:2]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end

for kk = 1:3
    subplot(2,3,kk+3)
plot(PVE_crossing([1:3 5:6],kk),Crossing_Angledev([1:3 5:6],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:6],kk)),max(PVE_crossing([1:3 5:6],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:6],kk),Crossing_Angledev([1:3 5:6],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(PVE_crossing([7:11],kk),Crossing_Angledev([7:11],kk),'*','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([7:11],kk)),max(PVE_crossing([7:11],kk)),100);
[P,S] = polyfit(PVE_crossing([7:11],kk),Crossing_Angledev([7:11],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(PVE_crossing([12:13 15:18],kk),Crossing_Angledev([12:13 15:18],kk),'s','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([12:13 15:18],kk)),max(PVE_crossing([12:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([12:13 15:18],kk),Crossing_Angledev([12:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
x = linspace(min(PVE_crossing([1:3 5:13 15:18],kk)),max(PVE_crossing([1:3 5:13 15:18],kk)),100);
[P,S] = polyfit(PVE_crossing([1:3 5:13 15:18],kk),Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on
ylim([15 40]);yticks([15:5:40]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end

% SNR vs. repro
for kk = 1:3
    subplot(2,3,kk+3)
plot(SNR(kk,[1:3 5:6]),Crossing_Angledev([1:3 5:6],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:6])),max(SNR(kk,[1:3 5:6])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:6])',Crossing_Angledev([1:3 5:6],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(SNR(kk,[7:11]),Crossing_Angledev([7:11],kk),'*','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[7:11])),max(SNR(kk,[7:11])),100);
[P,S] = polyfit(SNR(kk,[7:11])',Crossing_Angledev([7:11],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(SNR(kk,[12:13 15:18]),Crossing_Angledev([12:13 15:18],kk),'s','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[12:13 15:18])),max(SNR(kk,[12:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[12:13 15:18])',Crossing_Angledev([12:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Crossing_Angledev([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on
ylim([15 40]);yticks([15:5:40]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' Angle_D_e_v ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end

for kk = 1:3
    subplot(2,3,kk)
plot(SNR(kk,[1:3 5:6]),Crossing_SIMODF([1:3 5:6],kk),'o','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:6])),max(SNR(kk,[1:3 5:6])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:6])',Crossing_SIMODF([1:3 5:6],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(SNR(kk,[7:11]),Crossing_SIMODF([7:11],kk),'*','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[7:11])),max(SNR(kk,[7:11])),100);
[P,S] = polyfit(SNR(kk,[7:11])',Crossing_SIMODF([7:11],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
plot(SNR(kk,[12:13 15:18]),Crossing_SIMODF([12:13 15:18],kk),'s','linewidth',2,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[12:13 15:18])),max(SNR(kk,[12:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[12:13 15:18])',Crossing_SIMODF([12:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,':','linewidth',3,'color',resol_color(kk,:))
hold on
x = linspace(min(SNR(kk,[1:3 5:13 15:18])),max(SNR(kk,[1:3 5:13 15:18])),100);
[P,S] = polyfit(SNR(kk,[1:3 5:13 15:18])',Crossing_SIMODF([1:3 5:13 15:18],kk),1);
[f,delta] = polyconf(P,x,S);
plot(x,f,'k-.','linewidth',3)
hold on
ylim([0 2]);yticks([0:0.5:2]);xlim([0 40]);xticks([0:10:40])
set(gca,'fontname','calibri','fontsize',16,'linewidth',2);box off
ylabel(' SIM_O_D_F ','fontname','calibri','fontsize',18)
xlabel(' Non-WM volume fraction (%)','fontname','calibri','fontsize',18)
% legend(' 1.5 mm',' 2.0 mm',' 2.5 mm','fontname','calibri','fontsize',18);legend boxoff
end

%%
GLM_Model = 'Value ~ SNR + PVE';
for kk = 1:3
tbl = table(PVE_single([1:3 5:13 15:18],kk),Single_SIMODF([1:3 5:13 15:18],kk),SNR(kk,[1:3 5:13 15:18])','VariableNames',{'SNR','Value','PVE'});
mdl = fitglm(tbl,GLM_Model); 
SNR_coefficient(:,kk) = table2array(mdl.Coefficients(2,[1 4]));
PVE_coefficient(:,kk) = table2array(mdl.Coefficients(3,[1 4]));
end




