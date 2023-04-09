clearvars;clc
%%
rDSI = niftiread('eddy_corrected_data_pre_15.nii');
rDSI_dwi = rDSI(:,:,:,2:end);
rDSI_dwi_line = reshape(rDSI_dwi,[size(rDSI_dwi,1)*size(rDSI_dwi,2)*size(rDSI_dwi,3) 203]);
% 3-shell normalization (with 3 b-values)
bval = [1000 2000 3000];
% DSI b-vector in theory
rDSI_grad = generate_DSI_vectors(203);
% The radius of DSI sphere
rDSI_grad_radius = sqrt(rDSI_grad(:,1).^2+rDSI_grad(:,2).^2+rDSI_grad(:,3).^2); 
% Generate the shell vector
[icosa_tri icosa_vec] = trisphere(2);
for zz = 1:3
Shell_vector(:,:,zz) = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
end

rDSI_shell = zeros(size(rDSI_dwi_line,1),162,3);
for zz = 1:3
tic
for kk = 1:size(rDSI_dwi_line,1)
    % Normalization based on the ratio of maximal DSI radius
    shell_vector = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
    % For each voxel, create a rDSI_space (original q-space)
    rDSI_space = zeros(11,11,11);
    % Based on the assumed q-space, the shell vector is displaced to the center
    shell_vector = shell_vector+6; 
    for k = 1:size(rDSI_grad,1)
        rDSI_space(6+rDSI_grad(k,1),6+rDSI_grad(k,2),6+rDSI_grad(k,3)) = rDSI_dwi_line(kk,k);
    end
    rDSI_Space{kk} = rDSI_space;
end
toc

tic
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolopen = parpool(6);
    poolopenflag = 1;
else
    poolopenflag = 1;
end
N = size(rDSI_dwi_line,1);
parfor_progress(N);
parfor kk = 1:size(rDSI_dwi_line,1)
    shell_data = interp3(rDSI_Space{kk},shell_vector(1,:),shell_vector(2,:),shell_vector(3,:));
    rDSI_shell(kk,:,zz) = shell_data;
%     parfor_progress
end
parfor_progress(0);
if poolopenflag == 1
    delete(poolopen)
end
toc
end

rDSI_shell_3d = reshape(rDSI_shell,[size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3]);
info1 = niftiinfo('eddy_corrected_data_pre_15.nii');
info1.ImageSize = [size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3];
info1.Transform.T = [1.5,0,0,0;0,1.5,0,0;0,0,1.5,0;108,108,-45,1];
niftiwrite(single(rDSI_shell_3d),'eddy_corrected_data_pre_15_shell.nii',info1,'Compressed',1);

bval_all = reshape(repmat(bval,162,1),[486 1]);
writematrix(bval_all,'bval_162shell.txt')
bvec_all = cat(1,Shell_vector(:,:,1)',Shell_vector(:,:,2)',Shell_vector(:,:,3)');
dlmwrite('Shell_vector_162shell.txt',bvec_all,'delimiter',' ')


%%
clearvars
rDSI = niftiread('eddy_corrected_data_pre_20.nii');
rDSI_dwi = rDSI(:,:,:,2:end);
rDSI_dwi_line = reshape(rDSI_dwi,[size(rDSI_dwi,1)*size(rDSI_dwi,2)*size(rDSI_dwi,3) 203]);
% 3-shell normalization (with 3 b-values)
bval = [1000 2000 3000];
% DSI b-vector in theory
rDSI_grad = generate_DSI_vectors(203);
% The radius of DSI sphere
rDSI_grad_radius = sqrt(rDSI_grad(:,1).^2+rDSI_grad(:,2).^2+rDSI_grad(:,3).^2); 
% Generate the shell vector
[icosa_tri icosa_vec] = trisphere(2);
for zz = 1:3
Shell_vector(:,:,zz) = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
end

rDSI_shell = zeros(size(rDSI_dwi_line,1),162,3);
for zz = 1:3
tic
for kk = 1:size(rDSI_dwi_line,1)
    % Normalization based on the ratio of maximal DSI radius
    shell_vector = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
    % For each voxel, create a rDSI_space (original q-space)
    rDSI_space = zeros(11,11,11);
    % Based on the assumed q-space, the shell vector is displaced to the center
    shell_vector = shell_vector+6; 
    for k = 1:size(rDSI_grad,1)
        rDSI_space(6+rDSI_grad(k,1),6+rDSI_grad(k,2),6+rDSI_grad(k,3)) = rDSI_dwi_line(kk,k);
    end
    rDSI_Space{kk} = rDSI_space;
end
toc

tic
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolopen = parpool(6);
    poolopenflag = 1;
else
    poolopenflag = 1;
end
N = size(rDSI_dwi_line,1);
parfor_progress(N);
parfor kk = 1:size(rDSI_dwi_line,1)
    shell_data = interp3(rDSI_Space{kk},shell_vector(1,:),shell_vector(2,:),shell_vector(3,:));
    rDSI_shell(kk,:,zz) = shell_data;
%     parfor_progress
end
parfor_progress(0);
if poolopenflag == 1
    delete(poolopen)
end
toc
end

rDSI_shell_3d = reshape(rDSI_shell,[size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3]);
info1 = niftiinfo('eddy_corrected_data_pre_20.nii');
info1.ImageSize = [size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3];
info1.Transform.T = [1.5,0,0,0;0,1.5,0,0;0,0,1.5,0;108,108,-45,1];
niftiwrite(single(rDSI_shell_3d),'eddy_corrected_data_pre_20_shell.nii',info1,'Compressed',1);

%%
clearvars
rDSI = niftiread('eddy_corrected_data_pre_25.nii');
rDSI_dwi = rDSI(:,:,:,2:end);
rDSI_dwi_line = reshape(rDSI_dwi,[size(rDSI_dwi,1)*size(rDSI_dwi,2)*size(rDSI_dwi,3) 203]);
% 3-shell normalization (with 3 b-values)
bval = [1000 2000 3000];
% DSI b-vector in theory
rDSI_grad = generate_DSI_vectors(203);
% The radius of DSI sphere
rDSI_grad_radius = sqrt(rDSI_grad(:,1).^2+rDSI_grad(:,2).^2+rDSI_grad(:,3).^2); 
% Generate the shell vector
[icosa_tri icosa_vec] = trisphere(2);
for zz = 1:3
Shell_vector(:,:,zz) = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
end

rDSI_shell = zeros(size(rDSI_dwi_line,1),162,3);
for zz = 1:3
tic
for kk = 1:size(rDSI_dwi_line,1)
    % Normalization based on the ratio of maximal DSI radius
    shell_vector = 3.6056.*icosa_vec*sqrt(bval(zz)/4500);
    % For each voxel, create a rDSI_space (original q-space)
    rDSI_space = zeros(11,11,11);
    % Based on the assumed q-space, the shell vector is displaced to the center
    shell_vector = shell_vector+6; 
    for k = 1:size(rDSI_grad,1)
        rDSI_space(6+rDSI_grad(k,1),6+rDSI_grad(k,2),6+rDSI_grad(k,3)) = rDSI_dwi_line(kk,k);
    end
    rDSI_Space{kk} = rDSI_space;
end
toc

tic
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolopen = parpool(6);
    poolopenflag = 1;
else
    poolopenflag = 1;
end
N = size(rDSI_dwi_line,1);
parfor_progress(N);
parfor kk = 1:size(rDSI_dwi_line,1)
    shell_data = interp3(rDSI_Space{kk},shell_vector(1,:),shell_vector(2,:),shell_vector(3,:));
    rDSI_shell(kk,:,zz) = shell_data;
%     parfor_progress;
end
parfor_progress(0);
if poolopenflag == 1
    delete(poolopen)
end
toc
end

rDSI_shell_3d = reshape(rDSI_shell,[size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3]);
info1 = niftiinfo('eddy_corrected_data_pre_25.nii');
info1.ImageSize = [size(rDSI_dwi,1) size(rDSI_dwi,2) size(rDSI_dwi,3) 162*3];
info1.Transform.T = [1.5,0,0,0;0,1.5,0,0;0,0,1.5,0;108,108,-45,1];
niftiwrite(single(rDSI_shell_3d),'eddy_corrected_data_pre_25_shell.nii',info1,'Compressed',1);





