clear all;clc;close all
%% Read multiple files (1p5mm)
% Parent_dir = 'E:\dsi_data_7T_20200901_try\';
% ODF_recon_pre = '\ODF reconstruction\1p5mm comparison\Pre';
% ODF_recon_post = '\ODF reconstruction\1p5mm comparison\Post';
% ODF_register = '\Registration\1p5mm';
for aa = [1:13 15:18]
People{1} = 'NIPS_DSI_7T\HM';People{2} = 'NIPS_DSI_7T\HR';People{3} = 'NIPS_DSI_7T\KW';
People{4} = 'NIPS_DSI_7T\IS';People{5} = 'NIPS_DSI_7T\MA';People{6} = 'NIPS_DSI_7T\OS';
People{7} = 'NIPS_DSI_7T\OK';People{8} = 'NIPS_DSI_7T\MS';People{9} = 'NIPS_DSI_7T\SY';
People{10} = 'NIPS_DSI_7T\HY';People{11} = 'NIPS_DSI_7T\LWK';
People{12} = 'PRISMA3T_NH\CH';People{13} = 'PRISMA3T_NH\KHC';People{14} = 'PRISMA3T_NH\Peter';
People{15} = 'PRISMA3T_NH\RL';People{16} = 'PRISMA3T_NH\YP';People{17} = 'PRISMA3T_NH\YRC';
People{18} = 'PRISMA3T_NH\YR';

people = People{aa};
    disp(append('Start working for ',People{aa}, '_2.0mm at ', datestr(datetime(now,'ConvertFrom','datenum'))))
    if aa > 11
        Institution = 'NIPS_DSI_7T';
    else
        Institution = 'PRISMA3T_NH';
    end
% Pre
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people(13:end),'_pre_2p0mm'))
    load(dir('*17.mat').name)   
    if aa == 13
        WMmask = double(load_untouch_nii('c2rcT1_pre.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif aa == 4
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    elseif aa == 6
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/248;
        WMmask(WMmask <= 0.1) = nan; WMmask(~isnan(WMmask)) = 1; 
    else
        WMmask = double(load_untouch_nii('c2rcT1_pre_brain.nii').img);
        WMmask = WMmask/255;
        WMmask(WMmask <= 0.9) = nan; WMmask(~isnan(WMmask)) = 1; 
    end
    ODFratio = 0.4;

    for bb = 1:length(who)-29
        O{bb} =  eval(sprintf('odf%d',bb-1));
    end
    odfs = cell2mat(O); 
    clear O;
    clear odf0 odf1 odf2 odf3 odf4 odf5 odf6 odf7 odf8 odf9 odf10 odf11 odf12 odf13 odf14 odf15 odf16 odf17 odf18 odf19 odf20 odf21 odf22 odf23 odf24 odf25 odf26 odf27 odf28 odf29 odf30 odf31 odf32 odf33 odf34 odf35

index = [index0;index1;index2;index3;index4];
%the 1st fiber orientation at coordinate (10,10,10) stored as dir0
index0 = reshape(index(1,:),dimension);
index1 = reshape(index(2,:),dimension);
% Produce index0 into mask (fa0 > 0 for calculation)
Fa0 = reshape(fa0,[dimension]);Fa1 = reshape(fa1,[dimension]);Fa2 = reshape(fa2,[dimension]);
Fa3 = reshape(fa3,[dimension]);Fa4 = reshape(fa4,[dimension]);
Fa0(Fa0 == 0) = nan;Fa1(Fa1 == 0) = nan;Fa2(Fa2 == 0) = nan;Fa3(Fa3 == 0) = nan;Fa4(Fa4 == 0) = nan;
if aa == 11
Fa0 = Fa0.*imrotate(WMmask,180);Fa1 = Fa1.*imrotate(WMmask,180);Fa2 = Fa2.*imrotate(WMmask,180);
Fa3 = Fa3.*imrotate(WMmask,180);Fa4 = Fa4.*imrotate(WMmask,180);
elseif aa == 13
Fa0 = Fa0.*imrotate(WMmask,180);Fa1 = Fa1.*imrotate(WMmask,180);Fa2 = Fa2.*imrotate(WMmask,180);
Fa3 = Fa3.*imrotate(WMmask,180);Fa4 = Fa4.*imrotate(WMmask,180);
elseif aa == 17
Fa0 = Fa0.*imrotate(WMmask,180);Fa1 = Fa1.*imrotate(WMmask,180);Fa2 = Fa2.*imrotate(WMmask,180);
Fa3 = Fa3.*imrotate(WMmask,180);Fa4 = Fa4.*imrotate(WMmask,180);
else
Fa0 = Fa0.*WMmask;Fa1 = Fa1.*WMmask;Fa2 = Fa2.*WMmask;
Fa3 = Fa3.*WMmask;Fa4 = Fa4.*WMmask;
end

% figure(1)
% subplot(3,6,aa)
% imshow(imrotate(Fa0p(:,:,22),-90),[0 1]);colormap hot
% figure(2)
% subplot(3,6,aa)
% imshow(imrotate(WMmask(:,:,22),-90),[0 1]);colormap hot

facondition = Fa1./Fa0;
crossfibermask = zeros(dimension(1),dimension(2),dimension(3));
singlefibermask = zeros(dimension(1),dimension(2),dimension(3));
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if facondition(ii,jj,kk) >= ODFratio
                crossfibermask(ii,jj,kk) = 1;
            elseif facondition(ii,jj,kk) < ODFratio
                singlefibermask(ii,jj,kk) = 1;
            end
        end
    end
end
singlefibermask(singlefibermask == 0) = nan;
crossfibermask(crossfibermask == 0) = nan;

A1 = singlefibermask.*index0;
A2 = crossfibermask.*index0;
A3 = crossfibermask.*index1;

zvec = [0;0;1];
dirA1 = cell(dimension);
dirA2 = cell(dimension);
dirA3 = cell(dimension);

% Calculate angle rotated to z-axis and find rotation matrix
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(A1(ii,jj,kk))
                dirA1{ii,jj,kk} = odf_vertices(:,A1(ii,jj,kk)+1);
            else 
                dirA1{ii,jj,kk} = [];
            end
        end
    end
end
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(A2(ii,jj,kk))
                dirA2{ii,jj,kk} = odf_vertices(:,A2(ii,jj,kk)+1);
            else
                dirA2{ii,jj,kk} = [];
            end
        end
    end
end
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(A3(ii,jj,kk))
                dirA3{ii,jj,kk} = odf_vertices(:,A3(ii,jj,kk)+1);
            else
                dirA3{ii,jj,kk} = [];
            end
        end
    end
end
FF = reshape(fa0,[dimension(1) dimension(2) dimension(3)]);
[x1,x2,x3] = ind2sub(size(FF) ,find(FF >0));

% Post
    cd(append('E:\dsi_data_7T_20200901_try\SRC_Batch\',people(13:end),'_post_2p0mm'))
    load(dir('*17.mat').name)
    for cc = 1:bb
        O{cc} =  eval(sprintf('odf%d',cc-1));
    end
    odfsp = cell2mat(O); 
    clear O;
    clear odf0 odf1 odf2 odf3 odf4 odf5 odf6 odf7 odf8 odf9 odf10 odf11 odf12 odf13 odf14 odf15 odf16 odf17 odf18 odf19 odf20 odf21 odf22 odf23 odf24 odf25 odf26 odf27 odf28 odf29 odf30 odf31 odf32 odf33 odf34 odf35
    
% Produce index0 into mask (fa0 > 0 for calculation)
Fa0p = reshape(fa0,[dimension]);Fa1p = reshape(fa1,[dimension]);Fa2p = reshape(fa2,[dimension]);
Fa3p = reshape(fa3,[dimension]);Fa4p = reshape(fa4,[dimension]);
Fa0p(Fa0p == 0) = nan;Fa1p(Fa1p == 0) = nan;Fa2p(Fa2p == 0) = nan;Fa3p(Fa3p == 0) = nan;Fa4p(Fa4p == 0) = nan;
if aa == 11
Fa0p = Fa0p.*imrotate(WMmask,180);Fa1p = Fa1p.*imrotate(WMmask,180);Fa2p = Fa2p.*imrotate(WMmask,180);
Fa3p = Fa3p.*imrotate(WMmask,180);Fa4p = Fa4p.*imrotate(WMmask,180);
elseif aa == 13
Fa0p = Fa0p.*imrotate(WMmask,180);Fa1p = Fa1p.*imrotate(WMmask,180);Fa2p = Fa2p.*imrotate(WMmask,180);
Fa3p = Fa3p.*imrotate(WMmask,180);Fa4p = Fa4p.*imrotate(WMmask,180);
elseif aa == 17
Fa0p = Fa0p.*imrotate(WMmask,180);Fa1p = Fa1p.*imrotate(WMmask,180);Fa2p = Fa2p.*imrotate(WMmask,180);
Fa3p = Fa3p.*imrotate(WMmask,180);Fa4p = Fa4p.*imrotate(WMmask,180);
else
Fa0p = Fa0p.*WMmask;Fa1p = Fa1p.*WMmask;Fa2p = Fa2p.*WMmask;
Fa3p = Fa3p.*WMmask;Fa4p = Fa4p.*WMmask;
end

indexp = [index0;index1;index2;index3;index4];
%the 1st fiber orientation at coordinate (10,10,10) stored as dir0
index0p = reshape(index0,dimension);
index1p = reshape(index1,dimension);
B1 = singlefibermask.*index0p;
B2 = crossfibermask.*index0p;
B3 = crossfibermask.*index1p;
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(B1(ii,jj,kk))
                dirB1{ii,jj,kk} = odf_vertices(:,B1(ii,jj,kk)+1);
            else
                dirB1{ii,jj,kk} = [];
            end
        end
    end
end
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(B2(ii,jj,kk))
                dirB2{ii,jj,kk} = odf_vertices(:,B2(ii,jj,kk)+1);
            else
                dirB2{ii,jj,kk} = [];
            end
        end
    end
end
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isnan(B3(ii,jj,kk))
                dirB3{ii,jj,kk} = odf_vertices(:,B3(ii,jj,kk)+1);
            else
                dirB3{ii,jj,kk} = [];
            end
        end
    end
end
FFp = reshape(fa0,[dimension(1) dimension(2) dimension(3)]);
[x1p,x2p,x3p] = ind2sub(size(FFp) ,find(FFp >0));

%% Deviation angle
% Single-fiber groups
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isempty(dirA1{ii,jj,kk}) && ~isempty(dirB1{ii,jj,kk})
                Pa1(ii,jj,kk) = atan2(norm(cross(dirA1{ii,jj,kk},dirB1{ii,jj,kk})), dot(dirA1{ii,jj,kk},dirB1{ii,jj,kk}))*180/pi;
            else
                Pa1(ii,jj,kk) = nan;
            end
            if Pa1(ii,jj,kk) >= 90
                Pa1(ii,jj,kk) = 180-Pa1(ii,jj,kk);
            end
        end
    end
end
Pa1 = singlefibermask.*Pa1;
% Crossing-fiber groups
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            if ~isempty(dirA2{ii,jj,kk}) && ~isempty(dirB2{ii,jj,kk})
                Px4_1(ii,jj,kk) = atan2(norm(cross(dirA2{ii,jj,kk},dirB2{ii,jj,kk})), dot(dirA2{ii,jj,kk},dirB2{ii,jj,kk}))*180/pi;
                Px4_4(ii,jj,kk) = atan2(norm(cross(dirA3{ii,jj,kk},dirB3{ii,jj,kk})), dot(dirA3{ii,jj,kk},dirB3{ii,jj,kk}))*180/pi;
                if Px4_1(ii,jj,kk) >= 90
                    Px4_1(ii,jj,kk) = 180-Px4_1(ii,jj,kk);
                end
                if Px4_4(ii,jj,kk) >= 90
                    Px4_4(ii,jj,kk) = 180-Px4_4(ii,jj,kk);
                end
                    Pa4(ii,jj,kk) = mean([Px4_1(ii,jj,kk) Px4_4(ii,jj,kk)]);
            else
                Pa4(ii,jj,kk) = nan;
                Px4_1(ii,jj,kk) = nan;
                Px4_4(ii,jj,kk) = nan;
            end
        end
    end
end

%% ODF
Odfs = odfs(:,odfs(1, : ) ~= 0.00);
Odfsp = odfsp(:,odfsp(1, : ) ~= 0.00);
ODF = zeros(dimension(1),dimension(2),dimension(3),321);
ODF = ODF*nan;
for ii = 1:size(Odfs,2)
    ODF(x1(ii),x2(ii),x3(ii),:) = Odfs(:,ii);
end
ODFp = zeros(dimension(1),dimension(2),dimension(3),321);
ODFp = ODFp*nan;
for ii = 1:size(Odfsp,2)
    ODFp(x1p(ii),x2p(ii),x3p(ii),:) = Odfsp(:,ii);
end
% ODF similarity
for ii = 1:dimension(1)
    for jj = 1:dimension(2)
        for kk = 1:dimension(3)
            C = corrcoef(ODF(ii,jj,kk,:),ODFp(ii,jj,kk,:));
            ODFsim(ii,jj,kk) = C(2,1);
        end
    end
end

cd 'E:\dsi_data_7T_20200901_try\2p0mm_results_20220411'
people = People{aa};
save(append('ODF_info_20_',people(13:end)),'Fa0','Fa1','Fa2','Fa3','Fa4','singlefibermask','crossfibermask','facondition','index','dimension','odf_faces','odf_vertices');
save(append('ODF_info_20_',people(13:end),'p'),'Fa0p','Fa1p','Fa2p','Fa3p','Fa4p','indexp','dimension');
save(append('ODF_repro_20_',people(13:end)),'Pa1','Pa4','Px4_1','Px4_4','ODFsim');

disp(append('Finished working for ',People{aa}, '_2.0mm at ', datestr(datetime(now,'ConvertFrom','datenum'))))
clearvars
end

imshow(imrotate(Pa1(:,:,25),-90),[-2 10]);colormap jet;colorbar
imshow(imrotate(Pa4(:,:,25),-90),[-2 30]);colormap jet;colorbar

imshow(imrotate(ODF(:,:,25),-90),[]);colormap jet
imshow(imrotate(ODFp(:,:,25),-90),[]);colormap jet

imshow(imrotate(ODFsim(:,:,25),-90),[]);colormap jet
imshow(imrotate(ODFsim(:,:,25).*WMmask(:,:,25),-90),[]);colormap jet

mean(ODFsim.*singlefibermask,'all','omitnan')
mean(ODFsim.*crossfibermask,'all','omitnan')

%% Results of ODFSIM
cd('E:\dsi_data_7T_20200901_try\2p0mm_results_20220411')
for aa = [1:13 15:17]   
People{1} = 'NIPS_DSI_7T\HM';People{2} = 'NIPS_DSI_7T\HR';People{3} = 'NIPS_DSI_7T\KW';
People{4} = 'NIPS_DSI_7T\IS';People{5} = 'NIPS_DSI_7T\MA';People{6} = 'NIPS_DSI_7T\OS';
People{7} = 'NIPS_DSI_7T\OK';People{8} = 'NIPS_DSI_7T\MS';People{9} = 'NIPS_DSI_7T\SY';
People{10} = 'NIPS_DSI_7T\HY';People{11} = 'NIPS_DSI_7T\LWK';
People{12} = 'PRISMA3T_NH\CH';People{13} = 'PRISMA3T_NH\KHC';People{14} = 'PRISMA3T_NH\Peter';
People{15} = 'PRISMA3T_NH\RL';People{16} = 'PRISMA3T_NH\YP';People{17} = 'PRISMA3T_NH\YRC';
People{18} = 'PRISMA3T_NH\YR';

people = People{aa};
load(append('ODF_info_20_',people(13:end)),'singlefibermask','crossfibermask');
load(append('ODF_repro_20_',people(13:end)));
Single_SIMODF(aa) = mean(ODFsim.*singlefibermask,'all','omitnan');
Crossing_SIMODF(aa) = mean(ODFsim.*crossfibermask,'all','omitnan');
Single_Angledev(aa) = mean(Pa1,'all','omitnan');
Crossing_Angledev(aa) = mean(Pa4,'all','omitnan');
end
plot(Crossing_SIMODF)


V = niftiread('rcT1_pre_brain.nii');
spm_get_bbox('rcT1_pre_brain.nii')