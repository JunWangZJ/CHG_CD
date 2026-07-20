clc;
close all;
clear all;

%% data loading 
addpath(genpath(pwd));

dataName = 'IMAGE_NAME';   % 
X = imread(dataName+ "" + '1.png');
Y = imread(dataName+ "" + '2.png');
Ref_gt = imread(dataName+ "" + '3.png');
ref=Ref_gt(:,:,1);

X_or = double(X(:,:,1))+0.01;
Y_or = double(Y(:,:,1))+0.01;

X_nor = image_normlized(X_or,'sar')+0.001;
Y_nor = image_normlized(Y_or,'sar')+0.001;
[m,n] = size(X_nor);
tic;

%% Parameter Settings
ps = 1; % patch size = ps*2+1
pt = 1; % step
adj_rate = 0.01;  % Ar
Kp =round(m*n*adj_rate); 

Sn = round(m*n/((ps*2+1)^2)); 
Ksp = round(adj_rate*Sn); 

%% patch division
X_p = imageTodata(X_nor,ps,pt); 
Y_p = imageTodata(Y_nor,ps,pt);  
X_p_t = X_p';
Y_p_t = Y_p';

%% superpixel segmentation
Compactness = 1;
[sup_img] = SLIC_Cosegmentation_v2(X_nor,Y_nor,Sn,Compactness);
[t1_feature,t2_feature] = MSMfeature_extraction(sup_img,X_nor,Y_nor) ;
t1_f = t1_feature';
t2_f = t2_feature';

G=fspecial('average',ps*2+1);
Xm_nor=conv2(double(X_nor(:,:,1)),G,'same');
Ym_nor=conv2(double(Y_nor(:,:,1)),G,'same');
[t1_f2,t2_f2] = MSMfeature_extraction(sup_img,Xm_nor,Ym_nor) ; 
%% Information transmission
X_p_m = mean(X_p_t,2)';
Y_p_m = mean(Y_p_t,2)';

[X_agg_p2, id_pspX3] = CHG(X_p_t,t1_f,t1_f2,sup_img,ps,pt,Kp,Ksp);

[Y_agg_p2, id_pspY3] = CHG(Y_p_t,t2_f,t2_f2,sup_img,ps,pt,Kp,Ksp);

%%  change measure

X_fagg = X_p_m + X_agg_p2 + id_pspX3;    
Y_fagg = Y_p_m + Y_agg_p2 + id_pspY3;

f_dp = abs(log(X_fagg)-log(Y_fagg));
dif_sp  = dataToimage(f_dp,ps,pt,X_nor); 

%% feedback changes
[m,n] = size(dif_sp);
DI_sp = dif_sp(ps+1:m-ps, ps+1:n-ps);   
DI = DI_sp./max(DI_sp(:));
figure,imshow(DI);

indexedImage = gray2ind(uint8(DI*255), 256); 
rgbImage = ind2rgb(indexedImage, jet(256)); 
figure,imshow(rgbImage);

%% Otsu
level=graythresh(DI);
CM_OTSU = im2bw(DI, level);

%% visualization
TimeCT = toc;
ref = ref/max(ref(:));

[PRE, REC] = PR_plot(DI,ref,500);
[TPR, FPR]= Roc_plot(DI,ref, 500);
[AUC,~] = AUC_Diagdistance(TPR, FPR);
[AUP, ~] = AUP_Diagdistance(PRE, REC);

ROC_ASPG = [FPR; TPR];
PRC_ASPG = [REC, PRE];
% figure; plot(FPR,TPR);title('ROC curves'); 
% figure; plot(REC,PRE);title('PR curves');

[tp,fp,tn,fn,fplv,fnlv,~,~,pcc,kappa,imw]=performance(CM_OTSU,1*ref);
F1 = 2*tp/(2*tp + fp + fn);

[FP_x, FP_y] = find(imw==0);
[FN_x, FN_y] = find(imw==255);
CM_map_OTSU(:,:,1) = uint8(CM_OTSU)*255;
CM_map_OTSU(:,:,2) = uint8(CM_OTSU)*255;
CM_map_OTSU(:,:,3) = uint8(CM_OTSU)*255;
for i = 1 : max(size(FP_x))
    CM_map_OTSU(FP_x(i), FP_y(i), 1:3) = [255 0 0];
end
for i = 1 : max(size(FN_x))
    CM_map_OTSU(FN_x(i), FN_y(i), 1:3) = [0 255 0];
end

filename_OTSU = sprintf(dataName + "" + '_OTSU_ASPG_FN is %d; FP is %d; AUC is %4.4f; AUP is %4.4f;  PCC is %4.4f; F1 is %4.4f; KC is %4.4f.png',   fn, fp, AUC, AUP, pcc, F1, kappa)
