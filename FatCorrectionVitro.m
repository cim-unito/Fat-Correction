clear all
close all
clc
% We investigated the standard method (asym method) and four post-processing 
% approaches to overcome fat signal influences that enable correct CEST contrast 
% calculation and tumor pH measurements for in VITRO studies:
% - standard method (asym method)
% - positive method (#1)
% - linear method (#2)
% - lorentzian method (#3)
% - interpolation method (#4)
% Lorentzian CEST curve fitting was implemented with the open-source Matlab-based code
% (https://github.com/cest-sources/CEST_EVAL)
%% Uncomment to add the repo to the path CEST_EVAL
% addpath('utils', 'CEST_EVAL-master/levmar_fit')  
addpath('utils')  
%% Select Data
directory = uigetdir('','Select the patient folder');
cd(directory);

directory_morph=uigetdir('','Select the folder(number) where the anatomical reference is stored');
[img_morph , hdr_morph ] = readBruker(directory_morph);
img_morph=img_morph./hdr_morph.slope;
image=double(img_morph);
[size_x, size_y, size_z]=size(img_morph);

cd(directory)
directory_scan=uigetdir('','Select the folder(number) where the CEST pre image is stored');
[img_pre , hdr_pre ] = readBruker(directory_scan);
img_pre=img_pre./hdr_pre.slope;
image_zpre=double(img_pre);
[size_zx, size_zy, size_f, size_zz]= size(img_pre);
x=hdr_pre.freq_list';
num_ppm=hdr_pre.num_ppm;
[path, init1]=fileparts(directory_scan);
%% User selects parameters
R2factor=str2double(char(inputdlg('Which R2 factor to suppress Z-spectra noisy?' ,'',1, {'0.9'} )));
reg_weight=str2double(char(inputdlg('Which reg factor?' ,'',1, {'0.99'} )));

% Contrast agent params (Iopamidol)
list_iodu={'iopamidol'};
[composto]=listdlg('PromptString','Which contrast agent?','SelectionMode','single','ListString',list_iodu);
if strcmp(list_iodu(composto),'iopamidol')
    st_ppm=4.2;
    st_ppm2=5.5;
end

%% Segmentation Morphological Image
clear image_seg
clear image_morph
clear image_seg_hires 
clear level

for z=1:size_z
    [level,mask1(:,:,z)]=thresh_tool(mat2gray(squeeze(image(:,:,z,1))),'gray');
    image_seg(:,:,z)=immultiply(mask1(:,:,z),image(:,:,z,1));
    image_seg_hires(:,:,z)=image_seg(:,:,z);
    h_seg=figure; imshow(image_seg(:,:,z),[]);
    if size_x ~= size_y
        axis square
    end  
    title(['SEGMENTED MORPHOLOGICAL IMAGE: slice:' num2str(z)]);
end

for z=1:size_z
    mask_back(:,:,z)=imresize(squeeze(mask1(:,:,z)), [size_zx size_zy]);
end

mask_back(mask_back>=0.5)=1;
mask_back(mask_back<0.5)=0;

%% CEST Calculation
size_zz=size_z;
% Parameters for Lorentzian CEST curve fitting 
P.EVAL.w_fit=(-10:0.01:10)';
P.SEQ.w=sort(x');
P.EVAL.w_interp=P.SEQ.w;
P.EVAL.lowerlim_slices=1;
P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P.FIT.nIter     = 100;
P.FIT.modelnum  = 5;

% 1 pool
%  lb = [ 0  0  0  -1            0.0   0   -3];
%  ub = [ 1    1     10   +1            1   100  -1 ];
%  p0 = [ 1    0.9   1.4   0            0.1   25   -2  ];

% 5 pool
lb = [ 0.5  0.02  0.3  -1       0.0    0.4   +3     0.0    1    -4.5    0.0   10   -3     0.0   1    1     ];
ub = [ 1    1     10   +1       0.2    3     +4     0.4    5    -2        1   100  -1     0.2   3.5  2.5   ];
p0 = [ 1    0.9   1.4   0       0.025  0.5   3.5     0.02   3    -3.5    0.1   25   -2     0.01  1.5  2.2  ];

P.FIT.lower_limit_fit = lb; P.FIT.upper_limit_fit = ub; P.FIT.start_fit = p0;
P.SEQ.tp=3;
P.SEQ.B1=3;
P.SEQ.FREQ=300;
% 1 pool
% P.FIT.fitfunc='lorentzfit1pool'
% 5 pool
P.FIT.fitfunc='lorentzfit5pool';

PV= [P.SEQ.FREQ P.SEQ.tp P.SEQ.B1];

%-----------------------------------------------------------------%
% Image analysis - Z SPECTRUM calculation voxel x voxel           %
%-----------------------------------------------------------------%
% Normalization of the spectrum
for z=1:size_zz
    for j=1:size_zy
        for i=1:size_zx
            
            image_zpre(i,j,:)=image_zpre(i,j,:)./max(image_zpre(i,j,:));
           
        end
    end
end

% Matrix pre-allocation
image_zero_1=-300*ones(size_zx,size_zy,size_zz);
image_zero_2=-300*ones(size_zx,size_zy,size_zz);
image_r2pre=zeros(size_zx,size_zy,size_zz);
image_r2post=zeros(size_zx,size_zy,size_zz);

% Standard method (asym method) matrix pre-allocation
ST1_pre1=zeros(size_zx,size_zy,size_zz);
ST2_pre1=zeros(size_zx,size_zy,size_zz);
ST1_1=zeros(size_zx,size_zy,size_zz);
ST2_1=zeros(size_zx,size_zy,size_zz);

% Positive method (method #1) matrix pre-allocation
ST1_pre2=zeros(size_zx,size_zy,size_zz);
ST2_pre2=zeros(size_zx,size_zy,size_zz);
ST1_2=zeros(size_zx,size_zy,size_zz);
ST2_2=zeros(size_zx,size_zy,size_zz);

% Linear method (method #2) matrix pre-allocation
ST1_pre3=zeros(size_zx,size_zy,size_zz);
ST2_pre3=zeros(size_zx,size_zy,size_zz);

% Lorentzian method (method #3) matrix pre-allocation
ST1_pre4=zeros(size_zx,size_zy,size_zz);
ST2_pre4=zeros(size_zx,size_zy,size_zz);

% Fat fraction matrix pre-allocation
FF1=zeros(size_zx,size_zy,size_zz);
FF2=zeros(size_zx,size_zy,size_zz);
FF3=zeros(size_zx,size_zy,size_zz);

for cont_z=1:size_zz
    for cont_f=1:size_f
        for cont_i=1:size_zx
            for cont_j=1:size_zy
                Dpre(cont_f,cont_i,cont_j,cont_z)=image_zpre(cont_i,cont_j,cont_f,cont_z); 
            end
        end
    end
end

size_D=size(Dpre);
for cont_z=1:size_zz
    id_read = waitbar(cont_z/size_zz);
    for cont_i=1:size_zx
        for cont_j=1:size_zy
            clear Cpre
            clear Cpost        
            if mask_back(cont_i,cont_j,cont_z)>0
                for h=1:size_D(1)
                    Cpre(h,1)=x(h);
                    Cpre(h,2)=Dpre(h,cont_i,cont_j,cont_z);
                end
                
                cs1 = csaps(double(Cpre(:,1)), double(Cpre(:,2)), reg_weight, [],[]);
                % Water range for min [-2.49 2.5]
                [minval1,min_absc1]=fnmin(cs1,[-2.49 2.5]);
                % Fat range for min [-5 -2.5]
                [minval2,min_absc2]=fnmin(cs1,[-5 -2.5]);
                image_zero_1(cont_i,cont_j,cont_z)=min_absc1;
                % Calculate fat fraction FF with 3 different metrics
                FF1(cont_i,cont_j,cont_z)=minval1*100;
                FF2(cont_i,cont_j,cont_z)=((1-minval2)./(1-minval2+1-minval1))*100;
                FF3(cont_i,cont_j,cont_z)=(1-((minval2)./(minval2+minval1)))*100;

                % Correct the x values ​​(RF in ppm) for the shift found
                xc1=Cpre(:,1)-min_absc1;

                % Calculate R2: goodness of the interpolated curve to pass through
                % experimental points 
                ss_reg1=(Cpre(:,2)-fnval(cs1,Cpre(:,1))).^2;
                ss_tot1=(Cpre(:,2)-mean(Cpre(:,2))).^2;
                SS_reg1=sum(ss_reg1(:));
                SS_tot1=sum(ss_tot1(:));
                R2_pixelpre=1-(SS_reg1/SS_tot1);
                image_r2pre(cont_i,cont_j,cont_z)=R2_pixelpre;
                             
                % Calculate ST1 (st_ppm=4.2) standard method (asym method) 
                ST1_pre1(cont_i,cont_j,cont_z)=1-((fnval(cs1,st_ppm+min_absc1)/fnval(cs1,-st_ppm+min_absc1)));
                
                % Calculate ST1 (st_ppm=4.2) positive method (method #1) 
                ST1_pre2(cont_i,cont_j,cont_z)=1-(fnval(cs1,st_ppm+min_absc1));
                
                %----------------------------------------------------------
                % Start linear method (method #2) 
                delta=0.1;
                v=Cpre(:,1);
                v1=v(v<0);
                v2=v(v>0);
                y_abscissa_pos_pre3=(fnval(cs1,v2+min_absc1));
                abscissa_pos3=v2;
                y_abscissa_neg_pre3=flipud(fnval(cs1,v1+min_absc1));
                abscissa_neg3=v1;
                abscissa_neg3_lr=fliplr(abscissa_neg3');
                a=find(abs(abscissa_neg3_lr-(-7.5002))<0.001); 
                b=find(abs(abscissa_neg3_lr-(-2.0001))<0.001);
                y3_pre=interp1([abscissa_neg3_lr(a),abscissa_neg3_lr(b)],[y_abscissa_neg_pre3(a),y_abscissa_neg_pre3(b)],abscissa_neg3_lr(a:b),'linear','extrap');
                y_abscissa_neg_pre3(a:b)=y3_pre;
                y_pre3=[(y_abscissa_neg_pre3);fnval(cs1,min_absc1);y_abscissa_pos_pre3];
                ind1_pos=find(abs((v+min_absc1)-(st_ppm+min_absc1))<0.001);
                ind2_pos=find(abs((v+min_absc1)-(st_ppm2+min_absc1))<0.001);
                ind1_neg=find(abs((v+min_absc1)-(-st_ppm+min_absc1))<0.001);
                ind2_neg=find(abs((v+min_absc1)-(-st_ppm2+min_absc1))<0.001);
                % Calculate ST1 (st_ppm=4.2) linear method (method #2)
                ST1_pre3(cont_i,cont_j,cont_z)=1-(y_pre3(ind1_pos)./y_pre3(ind1_neg));
                % End linear method (method #2)
                %----------------------------------------------------------
                
                %----------------------------------------------------------
                % Start lorentzian method (method #3) 
                v=sort(Cpre(:,1));
                v1=v(v<0);
                v2=v(v>0);
                y_abscissa_pos_pre4=(fnval(cs1,v2+min_absc1));
                abscissa_pos4=v2;
                y_abscissa_neg_pre4=(fnval(cs1,v1+min_absc1));
                abscissa_neg4=v1;
                Zspectra14=Cpre;
                Zspectra4=Zspectra14(:,2);
                P.SEQ.stack_dim=size(Zspectra4);
                tmpzspec(:,1) = Zspectra4;
                [ret1, popt1, info1, covar1] = levmar(P.FIT.fitfunc, p0, tmpzspec, P.FIT.nIter, P.FIT.options, 'bc', lb, ub, sort(xc1), PV);
                [f, fZi, f1, f2, f3, f4, f5, f6, g, g2, g3, g4, g5]= fitmodelfunc_ANA(popt1,P);
                x_inter=P.SEQ.w;
                Zref.onlyWater=fZi(x_inter)-f1(x_inter);
                y_fit_pre4=Zref.onlyWater;
                y_pre4=[y_abscissa_neg_pre4;fnval(cs1,min_absc1);y_abscissa_pos_pre4];
                a=find(y_pre4==y_fit_pre4);
                if(sum(a)>0)
                    c=find_nearest(v(a),0);
                    v(c)
                    if (abs(v(c)+min_absc1)>2) 
                        c=find(abs(v+min_absc1-min_absc1) <0.0001);
                    end
                else
                    c=find(abs(v+min_absc1-min_absc1)< 0.0001);
                end
                y_pre4(1:c)=y_fit_pre4(1:c); 
                e=find(v==0);
                y_abscissa_neg_pre4=y_pre4(1:e-1);
                y_abscissa_pos_pre4=y_pre4(e+1:end);
                v=double(v);
                v=round(v,2);
                ind1_pos=find(abs((v+min_absc1)-(st_ppm+min_absc1))<0.0001);
                ind2_pos=find(abs((v+min_absc1)-(st_ppm2+min_absc1)) <0.0001);
                ind1_neg=find(abs((v+min_absc1)-(-st_ppm+min_absc1))<0.0001);
                ind2_neg=find(abs((v+min_absc1)-(-st_ppm2+min_absc1))<0.0001);
                % Calculate ST1 (st_ppm=4.2) lorentzian method (method #3)
                ST1_pre4(cont_i,cont_j,cont_z)=1-(y_pre4(ind1_pos)./y_pre4(ind1_neg));
                % End lorentzian method (method #3)
                %----------------------------------------------------------
                
                % Calculate ST2 (st_ppm2=5.5)
                % Calculate ST2 (st_ppm2=5.5) standard method (asym method) 
                ST2_pre1(cont_i,cont_j,cont_z)=1-((fnval(cs1,st_ppm2+min_absc1)/fnval(cs1,-st_ppm2+min_absc1)));
                % Calculate ST2 (st_ppm2=5.5) positive method (method #1)
                ST2_pre2(cont_i,cont_j,cont_z)=1-(fnval(cs1,st_ppm2+min_absc1));
                % Calculate ST2 (st_ppm2=5.5) linear method (method #2)
                ST2_pre3(cont_i,cont_j,cont_z)=1-(y_pre3(ind2_pos)./y_pre3(ind2_neg));
                % Calculate ST2 (st_ppm2=5.5) lorentzian method (method #3)
                ST2_pre4(cont_i,cont_j,cont_z)=1-(y_pre4(ind2_pos)./y_pre4(ind2_neg));
            end
        end
    end
end

%% Create the results folder called FF and save the fat fraction maps for the 3 metrics
cd(directory)
ResultFolder = 'FF';
if ~exist(ResultFolder, 'dir')
  mkdir(ResultFolder);
end
cd(ResultFolder)

% FF1
h4=figure('Position', get(0, 'Screensize'));
imshow(FF1,[0 100]);
colormap(jet(10))
colorbar
title("FF1")
saveas(h4,'FF1.png')
% FF2
h5=figure('Position', get(0, 'Screensize'));
imshow(FF2,[0 100]);
colormap(jet(10))
colorbar
title("FF2")
saveas(h5,'FF2.png')
% FF3
h6=figure('Position', get(0, 'Screensize'));
imshow(FF3,[0 100]);
colormap(jet(10))
colorbar
title("FF3")
saveas(h6,'FF3.png')

%----------------------------------------------%
% ROI selection in the morphological image     %
%----------------------------------------------%
num_roi=str2double(char(inputdlg('How many ROIs do you want to select?' ,'',1 )));
flag_roi = questdlg('You have a mask with ROIs saved?','','Yes','No','No');

if strcmp(flag_roi,'Yes')
    [file_mask,path_mask]=uigetfile('*.*','File with ROI mask');
    cd(path_mask);
    [pathstr,name,ext] = fileparts(file_mask);
    if strcmp(ext,'.nii')
        for r=1:num_roi

           if r>1
               [file_mask,path_mask]=uigetfile('*.*','File with ROI mask');
               cd(path_mask);
           end
        
           mask_tum=load_nii(file_mask);
           size_x_mask=mask_tum.hdr.dime.dim(2);
           size_y_mask=mask_tum.hdr.dime.dim(3);
           size_z_mask=mask_tum.hdr.dime.dim(4);
           mask_tum_image=double(mask_tum.img);
        
           for h=1:size_zz
            mask_slice_roi(:,:,h,r)=fliplr(rot90(mask_tum_image(:,:,h)));
           end

        end

      else
        load(file_mask) 
        num_roi=size(mask_slice_roi,4);
    end
else
    image_roi=zeros(size_zx, size_zy, size_zz);
    mask_slice_roi=zeros(size_zx, size_zy, size_zz,num_roi);
    for z=1:size_zz
        h_roi=figure; imshow(FF1(:,:,z),[]);
        colormap(jet(10))
        title('Select ROIs manually')
        for i=1:num_roi
            title(['Select ROIs manually - roi ', num2str(i) 'slice' num2str(z)])
            [BW2,x_p,y_p]=roipoly;
            for n=1:(size(x_p,1)-1)
                h_line=line([x_p(n) x_p(n+1)],[y_p(n) y_p(n+1)]);
                set(h_line,'color','g');
                set(h_line,'LineWidth',2);
            end
            cx=mean(x_p);
            cy=mean(y_p);
            text(cx, cy, num2str(i), 'FontSize', 16, 'Color','red')
            mask_slice_roi(:,:,z,i)=BW2;
            clear image_pro BW2 x_p y_p cx cy;
        end
    end
    cd(directory);
    cd(directory_morph);
    save('manual_roi_selected','mask_slice_roi','num_roi');
end

for z=1:size_zz
    for i=1:num_roi
        mask(:,:,i,z)=imresize(squeeze(mask_slice_roi(:,:,z,i)), [size_zx size_zy]);
        mask_back(:,:,z)=imresize(squeeze(mask1(:,:,z)), [size_zx size_zy]);
        image_r(:,:,z)=imresize(squeeze(image_seg(:,:,z)), [size_zx size_zy]);
    end
end
mask_back(mask_back>=0.5)=1;
mask_back(mask_back<0.5)=0;
mask(mask>=0.5)=1;
mask(mask<0.5)=0;
mapmask=zeros(size_zx,size_zy,size_zz);
for r=1:num_roi
    mapmask=mapmask+squeeze(mask(:,:,r,:));
end
image_roi=zeros(size_zx, size_zy, size_zz);
for i=1:num_roi
    image_roi=image_roi+squeeze(mask(:,:,i,:));
end
      
h0=figure('Position', get(0, 'Screensize'));
for z=1:size_zz
    BWoutline = bwperim(image_roi(:,:,z));
    Segout = image_r(:,:,z); 
    Segout(BWoutline) = max(image_r(:));
    subplot(2,ceil(size_zz/2),z)
    imshow(Segout,[]);
    title(['Slice: ' num2str(z)])
    clear BWoutline Segout
end
%%
%-------------------------------------------------------------------%
% image analysis - calculation of selected ROI average value        %
%-------------------------------------------------------------------%
cd([directory,'/',ResultFolder])

%--------------------------------------------------------------------------
% Start standard method (asym method) and positive method (#1)
clear Zmean_pre1
clear Zmean_post1
clear Zmean_pre2
clear Zmean_post2
zeroshift={};
Bpre_g={};
Mz_corr_pre_g={};
diff_ST_pre_g={};
size_z1=[];
for cont_r=1:num_roi
    clear Zmean_pre1
    clear Zmean_post1
    clear Zmean_pre2
    clear Zmean_post2
    clear k
    clear c
    for cont_z=1:size_zz   
        for cont_f=1:size_f
            Zmean_pre(cont_r,cont_z,cont_f)=sum(sum(immultiply(squeeze(image_zpre(:,:,cont_f,cont_z)),mask(:,:,cont_r,cont_z))))./sum(sum(mask(:,:,cont_r,cont_z)));
        end 
    end
    c=1;
    d=1;
    for cont_z=1:size_zz
        if(sum(isnan(Zmean_pre(cont_r,cont_z,:)))==0)
            Zmean_pre1(c,:)=Zmean_pre(cont_r,cont_z,:);
            k(c)=cont_z;
            c=c+1;
        else 
            Zmean_pre2(d,:)=Zmean_pre(cont_r,cont_z,:);
            d=d+1;
        end
    end
    size_z1=size(Zmean_pre1,1);
    Bpre=[x;squeeze(Zmean_pre1)];
    Bpre=Bpre';
    Bpre=sortrows(Bpre,1);
    for z=1:size_z1
        Bpre(:,z+1)=Bpre(:,z+1)/max(Bpre(:,z+1));
    end
    % Calculate the B0 shift and shift the data with linear interpolation
    for cont_z=1:size_z1
        cs_pre = csaps(Bpre(:,1),Bpre(:,cont_z+1), reg_weight,[] , []);
        [minval_pre,min_absc1]=fnmin(cs_pre,[-2,2]);
        x_c_pre=Bpre(:,1)-min_absc1;
        zeroshift(cont_r,k(cont_z))={min_absc1};
        delta=0.1;
        v=delta:delta:abs(max(x));
        y_abscissa_pos_pre=(fnval(cs_pre,v+min_absc1));
        abscissa_pos=v;
        y_abscissa_neg_pre=(fnval(cs_pre,-v+min_absc1));
        abscissa_neg=-v;
        diff_ST_pre(cont_z,:)=(1.-y_abscissa_pos_pre);
        Mz_corr_pre(cont_z,:) =([fliplr(y_abscissa_neg_pre),fnval(cs_pre,min_absc1),y_abscissa_pos_pre]);
        xcorr(cont_z,:)=([fliplr(abscissa_neg),0.0,abscissa_pos]);
    end
    Bpre_g(cont_r)={Bpre};
    Mz_corr_pre_g(cont_r)={Mz_corr_pre};
    diff_ST_pre_g(cont_r)={diff_ST_pre};
    size_zz1(cont_r)=size_z1;
    k1(cont_r)={k};
end

colori=[1 1 0; %3 giallo
    1 0 1; %4 magenta
    0 1 1;%5 ciano
    1 0 0;%6 rosso
    0 0.502 0;  %verde scuro
    0 1 0;%7 verde
    0 0 1;%8 blu
    0 0 0;%9 nero
    0.502 0 0;%marrone
    0.753 0.753 0.753;%grigio
    0 0.502 0.502; %verde acqua
    0.502 0 0.502; %viola
    0 0 0.502; %blu scuro
    0.502 0.502 0] ;%verde oliva  
h1=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(C(:,1),C(:,i+1),'o','Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('Normalized Intensity Values, a.u.','Fontname','axial','Fontsize',16);
        title('Z-spectra')
        hold on, plot(xcorr(i,:),D,'LineWidth',2,'Color',[colori(cont_r,:)])
        ylim([0 1.05])
        xlim([-10 10])
    end
end
legend('0% fat','0% fat', '20% fat','20% fat', '40% fat','40% fat','60% fat','60% fat','80% fat','80% fat')
legend('Location','SouthEast')
 
h2=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(abscissa_pos, F(i,:)*100,'LineWidth',2,'Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('ST%','Fontname','axial','Fontsize',16);
        title('ST%')
    end
end
legend('0% fat', '20% fat', '40% fat','60% fat', '80% fat')
legend('Location','NorthEast')

saveas(h1,[num2str(init1) '_zeta-spectra pre. method1' '.fig']);
saveas(h1,[num2str(init1) '_zeta-spectra pre. method1' '.jpeg']);
saveas(h2,[num2str(init1) '_ST pre. method1' '.fig']);
saveas(h2,[num2str(init1) '_ST pre. method1' '.jpeg']);
% End standard method (asym method) and positive method (#1)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Start linear method (#2)
clear Zmean_pre1
clear Zmean_post1
clear Zmean_pre2
clear Zmean_post2
zeroshift={};
Bpre_g={};
Bpost_g={};
Mz_corr_pre_g={};
Mz_corr_post_g={};
diff_ST_pre_g={};
diff_ST_post_g={};
size_z1=[];
for cont_r=1:num_roi
    clear Zmean_pre1
    clear Zmean_post1
    clear Zmean_pre2
    clear Zmean_post2
    clear k
    clear c
    for cont_z=1:size_zz
        for cont_f=1:size_f
            Zmean_pre(cont_r,cont_z,cont_f)=sum(sum(immultiply(squeeze(image_zpre(:,:,cont_f,cont_z)),mask(:,:,cont_r,cont_z))))./sum(sum(mask(:,:,cont_r,cont_z)));
        end 
    end 
    c=1;
    d=1;
    for cont_z=1:size_zz
        if(sum(isnan(Zmean_pre(cont_r,cont_z,:)))==0)
            Zmean_pre1(c,:)=Zmean_pre(cont_r,cont_z,:);
            k(c)=cont_z;
            c=c+1;   
        else 
            Zmean_pre2(d,:)=Zmean_pre(cont_r,cont_z,:);
            d=d+1;
        end
    end
    size_z1=size(Zmean_pre1,1); 
    Bpre=[x;squeeze(Zmean_pre1)];
    Bpre=Bpre';
    Bpre=sortrows(Bpre,1);
    for z=1:size_z1
        Bpre(:,z+1)=Bpre(:,z+1)/max(Bpre(:,z+1));
    end
    % Calculate the B0 shift and shift the data with linear interpolation
    for cont_z=1:size_z1
        clear min_p min_po min_absc_p min_absc_po
        cs_pre = csaps(Bpre(:,1),Bpre(:,cont_z+1), reg_weight,[] , []);
        [minval_pre,min_absc1]=fnmin(cs_pre,[-2,2]);
        x_c_pre=Bpre(:,1)-min_absc1;  
        zeroshift(cont_r,k(cont_z))={min_absc1};
        delta=0.1;
        v=delta:delta:abs(max(x));
        y_abscissa_pos_pre=(fnval(cs_pre,v+min_absc1));
        abscissa_pos=v;
        y_abscissa_neg_pre=fliplr(fnval(cs_pre,-v+min_absc1));
        abscissa_neg_lr=fliplr(abscissa_neg);
        a=find(abs(abscissa_neg_lr-(-7.5000))<0.00001);
        b=find(abs(abscissa_neg_lr- (-2.00000))<0.00001,1);
        y3_pre=interp1([abscissa_neg_lr(a),abscissa_neg_lr(b)],[y_abscissa_neg_pre(a),y_abscissa_neg_pre(b)],abscissa_neg_lr(a:b),'linear','extrap');
        y_abscissa_neg_pre(a:b)=(y3_pre);  
        diff_ST_pre(cont_z,:)=(1.-((y_abscissa_pos_pre)./fliplr(y_abscissa_neg_pre)));
        Mz_corr_pre(cont_z,:) =([((y_abscissa_neg_pre)),fnval(cs_pre,min_absc1),y_abscissa_pos_pre]);
        xcorr(cont_z,:)=([(abscissa_neg_lr),0.0,abscissa_pos]);
    end
    Bpre_g(cont_r)={Bpre};
    Mz_corr_pre_g(cont_r)={Mz_corr_pre};
    diff_ST_pre_g(cont_r)={diff_ST_pre};
    size_zz1(cont_r)=size_z1;
    k1(cont_r)={k};
end

h1=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(C(:,1),C(:,i+1),'o','Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('Normalized Intensity Values, a.u.','Fontname','axial','Fontsize',16);
        title('Z-spectra')
        hold on, plot(xcorr(i,:),D,'LineWidth',2,'Color',[colori(cont_r,:)])
        ylim([0 1.05])
        xlim([-10 10])
    end
end
legend('0% fat','0% fat', '20% fat','20% fat', '40% fat','40% fat','60% fat','60% fat','80% fat','80% fat')
legend('Location','SouthEast')
 
h2=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(abscissa_pos, F(i,:)*100,'LineWidth',2,'Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('ST%','Fontname','axial','Fontsize',16);
        title('ST%')
    end
end
legend('0% fat', '20% fat','40% fat','60% fat','80% fat')
legend('Location','NorthEast')

saveas(h1,[num2str(init1) '_zeta-spectra pre. method2' '.fig']);
saveas(h1,[num2str(init1) '_zeta-spectra pre. method2' '.jpeg']);
saveas(h2,[num2str(init1) '_ST pre. method2' '.fig']);
saveas(h2,[num2str(init1) '_ST pre. method2' '.jpeg']);

% End linear method (#2)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Start Lorentzian method (#3)
P.EVAL.w_fit=(-10:0.01:10)';
P.SEQ.w=sort(x');
P.EVAL.w_interp=P.SEQ.w;
P.EVAL.lowerlim_slices=1;
P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P.FIT.nIter     = 100;
P.FIT.modelnum  = 5;  

% 1 pool 
%  lb = [ 0  0  0  -1            0.0   0   -3];
%  ub = [ 1    1     10   +1            1   100  -1 ];
%  p0 = [ 1    0.9   1.4   0            0.1   25   -2  ];
% 5 pool
lb = [ 0.5  0.02  0.3  -1       0.0    0.4   +3     0.0    1    -4.5    0.0   10   -3     0.0   1    1     ];
ub = [ 1    1     10   +1       0.2    3     +4     0.4    5    -2        1   100  -1     0.2   3.5  2.5   ];
p0 = [ 1    0.9   1.4   0       0.025  0.5   3.5     0.02   3    -3.5    0.1   25   -2     0.01  1.5  2.2  ];

P.FIT.lower_limit_fit = lb; P.FIT.upper_limit_fit = ub; P.FIT.start_fit = p0;
P.SEQ.tp=3;
P.SEQ.B1=3;
P.SEQ.FREQ=300;
% 1 pool 
% P.FIT.fitfunc='lorentzfit1pool';
% 5 pool 
P.FIT.fitfunc='lorentzfit5pool';
PV= [P.SEQ.FREQ P.SEQ.tp P.SEQ.B1];

clear Zmean_pre1
clear Zmean_pre2
zeroshift={};
Bpre_g={};
Mz_corr_pre_g={};
diff_ST_pre_g={};
filename='analisys.xlsx';
size_z1=[];
for cont_r=1:num_roi
    clear Zmean_pre1
    clear Zmean_post1
    clear Zmean_pre2
    clear Zmean_post2
    clear k
    clear c
    for cont_z=1:size_zz
        for cont_f=1:size_f
            Zmean_pre(cont_r,cont_z,cont_f)=sum(sum(immultiply(squeeze(image_zpre(:,:,cont_f,cont_z)),mask(:,:,cont_r,cont_z))))./sum(sum(mask(:,:,cont_r,cont_z)));
        end
    end
    c=1;
    d=1;
    for cont_z=1:size_zz
        if(sum(isnan(Zmean_pre(cont_r,cont_z,:)))==0)
                Zmean_pre1(c,:)=Zmean_pre(cont_r,cont_z,:);
                k(c)=cont_z;
                c=c+1;
            else 
                Zmean_pre2(d,:)=Zmean_pre(cont_r,cont_z,:);
                d=d+1;
        end
    end
    size_z1=size(Zmean_pre1,1);
    Bpre=[x;squeeze(Zmean_pre1)];
    Bpre=Bpre';
    Bpre=sortrows(Bpre,1);
    for z=1:size_z1
        Bpre(:,z+1)=Bpre(:,z+1)/max(Bpre(:,z+1));
    end
    % Calculate the B0 shift and shift the data with linear interpolation
    for cont_z=1:size_z1
        cs_pre = csaps(Bpre(:,1),Bpre(:,cont_z+1), reg_weight,[] , []);
        [minval_pre,min_absc1]=fnmin(cs_pre,[-2,2]);
        x_c_pre=Bpre(:,1)-min_absc1;
        zeroshift(cont_r,k(cont_z))={min_absc1};
        delta=0.1;
        v=delta:delta:abs(max(x));
        v1=[fliplr(-v),0,v];
        Zspectra1=Bpre;
        Zspectra=Zspectra1(:,2);
        P.SEQ.stack_dim=size(Zspectra);
        tmpzspec(:,1) = Zspectra;
        [ret1, popt1, info1, covar1] = levmar(P.FIT.fitfunc, p0, tmpzspec, P.FIT.nIter, P.FIT.options, 'bc', lb, ub, x_c_pre, PV);
        [f, fZi, f1, f2, f3, f4, f5, f6, g, g2, g3, g4, g5]= fitmodelfunc_ANA(popt1,P);
        x_inter=P.SEQ.w;
        Zref.onlyWater=fZi(x_inter)-f1(x_inter);
        yfit_int_pre=interp1(Bpre(:,1),Zref.onlyWater,v1);
        y_abscissa_pos_pre=(fnval(cs_pre,v+min_absc1));
        abscissa_pos=v;
        y_abscissa_neg_pre=(fnval(cs_pre,-v+min_absc1));
        abscissa_neg=-v;
        ypre=[fliplr(y_abscissa_neg_pre),fnval(cs_pre,min_absc1),y_abscissa_pos_pre];
        ypre1=ypre;
        a=find(ypre==yfit_int_pre);
        if(sum(a)>0)
            c=find_nearest(v1(a),0);
            v1(c)
            if (abs(v1(c)+min_absc1)>2)
                c=find(abs(v1+min_absc1-min_absc1)<0.0001);
            end
        else
            c=find(abs(v1+min_absc1-min_absc1)<0.0001);
        end
        ypre(1:c)=yfit_int_pre(1:c);  
        e=find(v1==0);
        y_abscissa_neg_pre=ypre(1:e-1);
        y_abscissa_pos_pre=ypre(e+1:end);
        diff_ST_pre(cont_z,:)=(1.-((y_abscissa_pos_pre)./y_abscissa_neg_pre));     
        st1_pre(cont_z)=diff_ST_pre(cont_z,st_ppm*10)*100;        
        st2_pre(cont_z)=diff_ST_pre(cont_z,st_ppm2*10)*100;
        ratio_pre(cont_z)=st1_pre(cont_z)/st2_pre(cont_z);
        Mz_corr_pre(cont_z,:) =([(y_abscissa_neg_pre),fnval(cs_pre,min_absc1),y_abscissa_pos_pre]);
        xcorr(cont_z,:)=([fliplr(abscissa_neg),0.0,abscissa_pos]);
    end
    Bpre_g(cont_r)={Bpre};
    Mz_corr_pre_g(cont_r)={Mz_corr_pre};
    diff_ST_pre_g(cont_r)={diff_ST_pre};
    size_zz1(cont_r)=size_z1;
    k1(cont_r)={k};
end
 
h1=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(C(:,1),C(:,i+1),'o','Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('Normalized Intensity Values, a.u.','Fontname','axial','Fontsize',16);
        title('Z-spectra')
        hold on, plot(xcorr(i,:),D,'LineWidth',2,'Color',[colori(cont_r,:)])
        ylim([0 1.05])
        xlim([-10 10])
    end
end
legend('0% fat','0% fat', '20% fat','20% fat', '40% fat','40% fat','60% fat','60% fat','80% fat','80% fat')
legend('Location','SouthEast')
 
h2=figure;
for cont_r=1:num_roi
    clear C 
    clear D E F
    C=cell2mat(Bpre_g(cont_r));
    D=cell2mat(Mz_corr_pre_g(cont_r));
    F=cell2mat( diff_ST_pre_g(cont_r));
    E=cell2mat(k1(cont_r));
    for i=1:size_zz1(cont_r)
        hold on
        plot(abscissa_pos, F(i,:)*100,'LineWidth',2,'Color',[colori(cont_r,:)])
        xlabel('Sat. offset, ppm','Fontname','axial','Fontsize',16);
        ylabel('ST%','Fontname','axial','Fontsize',16);
        title('ST%')
    end
end
legend('0% fat', '20% fat','40% fat','60% fat','80% fat')
legend('Location','NorthEast')

H1=getframe(h1);
H2=getframe(h2);
saveas(h1,[num2str(init1) '_zeta-spectra pre. method3' '.fig']);
saveas(h1,[num2str(init1) '_zeta-spectra pre. method3' '.jpeg']);
saveas(h2,[num2str(init1) '_ST pre. method3' '.fig']);
saveas(h2,[num2str(init1) '_ST pre. method3' '.jpeg']);
% End lorentzian method (#3)
%--------------------------------------------------------------------------
%% ST Ratio
ST1_1=ST1_pre1;
ST1_2=ST1_pre2;
ST1_3=ST1_pre3;
ST1_4=ST1_pre4;
ST2_1=ST2_pre1;
ST2_2=ST2_pre2;
ST2_3=ST2_pre3;
ST2_4=ST2_pre4;
STratio1 = zeros(size_zx,size_zy,size_zz);
STratio2 = zeros(size_zx,size_zy,size_zz);
STratio3 = zeros(size_zx,size_zy,size_zz);
STratio4 = zeros(size_zx,size_zy,size_zz);

for cont_z=1:size_zz
    for cont_j=1:size_zy
        for cont_i=1:size_zx
                if (mask_back(cont_i,cont_j,cont_z)>0)
                    STratio1(cont_i,cont_j,cont_z)=ST1_1(cont_i,cont_j,cont_z)/ST2_1(cont_i,cont_j,cont_z);
                    STratio2(cont_i,cont_j,cont_z)=ST1_2(cont_i,cont_j,cont_z)/ST2_2(cont_i,cont_j,cont_z);
                    STratio3(cont_i,cont_j,cont_z)=ST1_3(cont_i,cont_j,cont_z)/ST2_3(cont_i,cont_j,cont_z);
                    STratio4(cont_i,cont_j,cont_z)=ST1_4(cont_i,cont_j,cont_z)/ST2_4(cont_i,cont_j,cont_z);
                else
                    STratio1(cont_i,cont_j,cont_z)=0;
                    STratio2(cont_i,cont_j,cont_z)=0;
                    STratio3(cont_i,cont_j,cont_z)=0;
                    STratio4(cont_i,cont_j,cont_z)=0;
                end
        end
    end
end
STratio1(STratio1>10)=10;
STratio2(STratio2>10)=10;
STratio3(STratio3>10)=10;
STratio4(STratio4>10)=10;
STratio1(image_r2pre < R2factor)=0;
STratio2(image_r2pre < R2factor)=0;
STratio3(image_r2pre < R2factor)=0;
STratio4(image_r2pre < R2factor)=0;

%--------------------------------------------------------------------------
% Start interpolation method (#4)
% Cubic smoothing spline interpolation 
load pp_spline.mat
FF1_x = 0:100;
y_7 = fnval(pp7,FF1_x);
y_65 = fnval(pp65,FF1_x);
y_7 = y_7 + 0.4671;
y_65 = y_65-0.3328;

% Interpolated curves
n = 200;
map_matrix = zeros(n, length(FF1_x));
for k = 1:length(FF1_x)
  map_matrix(:,k) = linspace(y_65(k),y_7(k),n);
end
STratio5 = zeros(size_zx,size_zy,size_zz);
for cont_z=1:size_zz
    for cont_j=1:size_zy
        for cont_i=1:size_zx
            [c, index] = min(abs(FF1_x-(FF1(cont_i,cont_j,cont_z)-3.38)));
            [c1, index1] = min(abs(map_matrix(:,index)-STratio1(cont_i,cont_j,cont_z)));
            STratio5(cont_i,cont_j,cont_z) = map_matrix(index1,1);
        end
    end
end
% STratio1 = STratio5;
% End interpolation method (#4)
%--------------------------------------------------------------------------

%% Calcolo pH
pH1 = zeros(size_zx,size_zy,size_zz);
pH2 = zeros(size_zx,size_zy,size_zz);
pH3 = zeros(size_zx,size_zy,size_zz);
pH4 = zeros(size_zx,size_zy,size_zz);
pH5 = zeros(size_zx,size_zy,size_zz);
for cont_z=1:size_zz
    for cont_j=1:size_zy
        for cont_i=1:size_zx
            pH1(cont_i,cont_j,cont_z)=1.2401*(STratio1(cont_i,cont_j,cont_z))+5.3205; % retta di calibrazione
            pH2(cont_i,cont_j,cont_z)=1.6884*(STratio2(cont_i,cont_j,cont_z))+4.7764; % retta di calibrazione
            pH3(cont_i,cont_j,cont_z)=0.63*(STratio3(cont_i,cont_j,cont_z))+5.6604; % retta di calibrazione
            pH4(cont_i,cont_j,cont_z)=1.2286*(STratio4(cont_i,cont_j,cont_z))+5.3339; % retta di calibrazione
            pH5(cont_i,cont_j,cont_z)=1.2401*(STratio5(cont_i,cont_j,cont_z))+5.3205; % retta di calibrazione

        end
    end
end

pH1(STratio1>1.7 & STratio1<2.5)=7.4;
pH1(STratio1>0 & STratio1<0.6)=6.0;
pH1(STratio1<=0 | STratio1>=2.5)=NaN;
pH2(STratio2>1.7 & STratio2<2.5)=7.4;
pH2(STratio2>0 & STratio2<0.6)=6.0;
pH2(STratio2<=0 | STratio2>=2.5)=NaN;
pH3(STratio3>2.6 & STratio3<3.5)=7.4;
pH3(STratio3>0 & STratio3<0.6)=6.0;
pH3(STratio3<=0 | STratio3>=3.5)=NaN;
pH4(STratio4>1.7 & STratio4<2.5)=7.4;
pH4(STratio4>0 & STratio4<0.6)=6.0;
pH4(STratio4<=0 | STratio4>=2.5)=NaN;
pH5(STratio5>1.7 & STratio5<2.5)=7.4;
pH5(STratio5>0 & STratio5<0.6)=6.0;
pH5(STratio5<=0 | STratio5>=2.5)=NaN;

pH1 = pH1.*mask_back;
pH2 = pH2.*mask_back;
pH3 = pH3.*mask_back;
pH4 = pH4.*mask_back;
pH5 = pH5.*mask_back;

close(id_read)

%% Save images
cd([directory,'/',ResultFolder])
load('cest_colormap.mat')
load('iopa_pH_map.mat')

h3=figure('Position', get(0, 'Screensize'));
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(medfilt2(image_zero_1(:,:,i)), [-1 1]);
    colormap hot
    colorbar
    title(['"ZERO" SHIFT MAP PRE. Slice:' num2str(i)])
end

h5=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r2pre(:,:,i), [0.97 1]);
    colormap(mycmap)
    colorbar
    title(['R2 map PRE . Slice:' num2str(i)])
end

H3=getframe(h3);
H5=getframe(h5);
imwrite(H3.cdata,[num2str(init1) '_zero_shift_map_pre' '.jpeg'] , 'jpeg');
imwrite(H5.cdata,[num2str(init1) '_R2 map pre' '.jpeg'] , 'jpeg');
saveas(h3,[num2str(init1) '_zero_shift_map_pre' '.fig']);
saveas(h5,[num2str(init1) '_R2 mappre' '.fig']);
%%
for cont_r=1:num_roi
    for cont_z=1:size_zz
        % ST1
        ST1_roi1(:,:,cont_r,cont_z)=immultiply(squeeze(ST1_1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST1_roi2(:,:,cont_r,cont_z)=immultiply(squeeze(ST1_2(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST1_roi3(:,:,cont_r,cont_z)=immultiply(squeeze(ST1_3(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST1_roi4(:,:,cont_r,cont_z)=immultiply(squeeze(ST1_4(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        image_zero_1_roi(:,:,cont_r,cont_z)=immultiply(squeeze(image_zero_1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        % ST2
        ST2_roi1(:,:,cont_r,cont_z)=immultiply(squeeze(ST2_1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST2_roi2(:,:,cont_r,cont_z)=immultiply(squeeze(ST2_2(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST2_roi3(:,:,cont_r,cont_z)=immultiply(squeeze(ST2_3(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        ST2_roi4(:,:,cont_r,cont_z)=immultiply(squeeze(ST2_4(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        % STRatio
        STratio_roi1(:,:,cont_r,cont_z)=immultiply(squeeze(STratio1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        STratio_roi2(:,:,cont_r,cont_z)=immultiply(squeeze(STratio2(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        STratio_roi3(:,:,cont_r,cont_z)=immultiply(squeeze(STratio3(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        STratio_roi4(:,:,cont_r,cont_z)=immultiply(squeeze(STratio4(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        STratio_roi5(:,:,cont_r,cont_z)=immultiply(squeeze(STratio5(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        % pH
        pH_roi1(:,:,cont_r,cont_z)=immultiply(squeeze(pH1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        pH_roi2(:,:,cont_r,cont_z)=immultiply(squeeze(pH2(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        pH_roi3(:,:,cont_r,cont_z)=immultiply(squeeze(pH3(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        pH_roi4(:,:,cont_r,cont_z)=immultiply(squeeze(pH4(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        pH_roi5(:,:,cont_r,cont_z)=immultiply(squeeze(pH5(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        % Fat Fraction FF
        FF1_roi(:,:,cont_r,cont_z)=immultiply(squeeze(FF1(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        FF2_roi(:,:,cont_r,cont_z)=immultiply(squeeze(FF2(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
        FF3_roi(:,:,cont_r,cont_z)=immultiply(squeeze(FF3(:,:,cont_z)),squeeze(mask(:,:,cont_r,cont_z)));
    end
end
%% Calculate params
for i=1:num_roi
    a0 = mask(:,:,i);
    pix_tot_roi(i) = sum(a0(:));

    % standard method (asym method)
    a=pH_roi1(:,:,i,1);
    b=a(a~=0 & ~isnan(a));
    pix_enh_roi1(i) = length(b);
    mpH1(i)=mean(b);
    devpH1(i)=std(b);
   
    c=ST1_roi1(:,:,i,1);
    d=c(c~=0);
    mST11(i)=mean(d);
    devST11(i)=std(d);
    
    clear c d
    
    c=ST2_roi1(:,:,i,1);
    d=c(c~=0);
    mS21(i)=mean(d);
    devST21(i)=std(d);
    
    clear c d

    c=STratio_roi1(:,:,i,1);
    d=c(c~=0);
    mSTR1(i)=mean(d);
    devSTR1(i)=std(d);

    clear a0 a b c d
    
    % positive method (#1)
    a=pH_roi2(:,:,i,1);
    b=a(a~=0 & ~isnan(a));
    pix_enh_roi2(i) = length(b);
    mpH2(i)=mean(b);
    devpH2(i)=std(b);
    
    c=ST1_roi2(:,:,i,1);
    d=c(c~=0);
    mST12(i)=mean(d);
    devST12(i)=std(d);
    
    clear c d
    
    c=ST2_roi2(:,:,i,1);
    d=c(c~=0);
    mS22(i)=mean(d);
    devST22(i)=std(d);
    
    clear c d

    c=STratio_roi2(:,:,i,1);
    d=c(c~=0);
    mSTR2(i)=mean(d);
    devSTR2(i)=std(d);

    clear a b c d

    % linear method (#2)
    a=pH_roi3(:,:,i,1);
    b=a(a~=0 & ~isnan(a));
    pix_enh_roi3(i) = length(b);
    mpH3(i)=mean(b);
    devpH3(i)=std(b);
    
    c=ST1_roi3(:,:,i,1);
    d=c(c~=0);
    mST13(i)=mean(d);
    devST13(i)=std(d);
    
    clear c d
    
    c=ST2_roi3(:,:,i,1);
    d=c(c~=0);
    mS23(i)=mean(d);
    devST23(i)=std(d);
    
    clear c d

    c=STratio_roi3(:,:,i,1);
    d=c(c~=0);
    mSTR3(i)=mean(d);
    devSTR3(i)=std(d);

    clear a b c d

    % Lorentzian method (#3)
    a=pH_roi4(:,:,i,1);
    b=a(a~=0 & ~isnan(a));
    pix_enh_roi4(i) = length(b);
    mpH4(i)=mean(b);
    devpH4(i)=std(b);
    
    c=ST1_roi4(:,:,i,1);
    d=c(c~=0);
    mST14(i)=mean(d);
    devST14(i)=std(d);
    
    clear c d
    
    c=ST2_roi4(:,:,i,1);
    d=c(c~=0);
    mS24(i)=mean(d);
    devST24(i)=std(d);
    
    clear c d

    c=STratio_roi4(:,:,i,1);
    d=c(c~=0);
    mSTR4(i)=mean(d);
    devSTR4(i)=std(d);

    clear a b c d

    % interpolation method (#4)
    a=pH_roi5(:,:,i,1);
    b=a(a~=0 & ~isnan(a));
    pix_enh_roi5(i) = length(b);
    mpH5(i)=mean(b);
    devpH5(i)=std(b);
    
    c=STratio_roi5(:,:,i,1);
    d=c(c~=0);
    mSTR5(i)=mean(d);
    devSTR5(i)=std(d);

    clear a b c d

    % Fat fraction
    c = FF1_roi(:,:,i,1);
    d=c(c~=0);
    mFF1(i)=mean(d);
    devFF1(i)=std(d);

    clear c d

    c = FF2_roi(:,:,i,1);
    d=c(c~=0);
    mFF2(i)=mean(d);
    devFF2(i)=std(d);

    clear c d

    c = FF3_roi(:,:,i,1);
    d=c(c~=0);
    mFF3(i)=mean(d);
    devFF3(i)=std(d);

    clear c d
      
end
% Fraction pixel
FP_1 = (pix_enh_roi1./pix_tot_roi).*100;
FP_2 = (pix_enh_roi2./pix_tot_roi).*100;
FP_3 = (pix_enh_roi3./pix_tot_roi).*100;
FP_4 = (pix_enh_roi4./pix_tot_roi).*100;
FP_5 = (pix_enh_roi5./pix_tot_roi).*100;
%% construction of parametric maps for each slice
% standard method (asym method)
mapST1_1=zeros(size_zx,size_zy,size_zz);
mapST2_1=zeros(size_zx,size_zy,size_zz);
mapSTratio1=zeros(size_zx,size_zy,size_zz);
mapzero1=zeros(size_zx,size_zy,size_zz);
% positive method (#1)
mapST1_2=zeros(size_zx,size_zy,size_zz);
mapST2_2=zeros(size_zx,size_zy,size_zz);
mapSTratio2=zeros(size_zx,size_zy,size_zz);
% linear method (#2)
mapST1_3=zeros(size_zx,size_zy,size_zz);
mapST2_3=zeros(size_zx,size_zy,size_zz);
mapSTratio3=zeros(size_zx,size_zy,size_zz);
% lorentzian method (#3)
mapST1_4=zeros(size_zx,size_zy,size_zz);
mapST2_4=zeros(size_zx,size_zy,size_zz);
mapSTratio4=zeros(size_zx,size_zy,size_zz);
% interpolation method (#4)
mapSTratio5=zeros(size_zx,size_zy,size_zz);
% pH
map_pH1=zeros(size_zx,size_zy,size_zz);
map_pH2=zeros(size_zx,size_zy,size_zz);
map_pH3=zeros(size_zx,size_zy,size_zz);
map_pH4=zeros(size_zx,size_zy,size_zz);
map_pH5=zeros(size_zx,size_zy,size_zz);


for r=1:num_roi
    % ST1
    mapST1_1=mapST1_1+squeeze(ST1_roi1(:,:,r,:));
    mapST1_2=mapST1_2+squeeze(ST1_roi2(:,:,r,:));
    mapST1_3=mapST1_3+squeeze(ST1_roi3(:,:,r,:));
    mapST1_4=mapST1_4+squeeze(ST1_roi4(:,:,r,:));
    
    mapzero1=mapzero1+squeeze(image_zero_1_roi(:,:,r,:));
    % ST2
    mapST2_1=mapST2_1+squeeze(ST2_roi1(:,:,r,:));
    mapST2_2=mapST2_2+squeeze(ST2_roi2(:,:,r,:));
    mapST2_3=mapST2_3+squeeze(ST2_roi3(:,:,r,:));
    mapST2_4=mapST2_4+squeeze(ST2_roi4(:,:,r,:));
    % ST ratio
    mapSTratio1=mapSTratio1+squeeze(STratio_roi1(:,:,r,:));
    mapSTratio2=mapSTratio2+squeeze(STratio_roi2(:,:,r,:));
    mapSTratio3=mapSTratio3+squeeze(STratio_roi3(:,:,r,:));
    mapSTratio4=mapSTratio4+squeeze(STratio_roi4(:,:,r,:));
    mapSTratio5=mapSTratio5+squeeze(STratio_roi5(:,:,r,:));
    % pH
    map_pH1=map_pH1+squeeze(pH_roi1(:,:,r,:));
    map_pH2=map_pH2+squeeze(pH_roi2(:,:,r,:));
    map_pH3=map_pH3+squeeze(pH_roi3(:,:,r,:));
    map_pH4=map_pH4+squeeze(pH_roi4(:,:,r,:));
    map_pH5=map_pH4+squeeze(pH_roi4(:,:,r,:));
end
mapzero1(mapzero1==0)=-300;
%% Save map figures
%--------------------------------------------------------------------------
% standard method (asym method)
h7=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2)  ,i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST1_1(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST1_1(:,:,i)>0).*0.7); 
    title(['deltaST at ' num2str(st_ppm) ' ppm. Slice:' num2str(i)])
end

h8=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST2_1(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST2_1(:,:,i)>0).*0.7);
    title(['deltaST at ' num2str(st_ppm2) ' ppm. Slice:' num2str(i)])
end

h9=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapSTratio1(:,:,i), [0 3]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapSTratio1(:,:,i)>0).*0.7);
    title( ['Delta ratio. Slice:' num2str(i)])
end

h11=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(pH1(:,:,i), [5 8]);
    set(gcf(),'Colormap',pHmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (pH1(:,:,i)>0).*0.7);
    title(['Image pH. Slice:' num2str(i)])
end

H7=getframe(h7);
H8=getframe(h8);
H9=getframe(h9);
imwrite(H8.cdata,[num2str(init1) '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method_asym' '.jpg'] , 'jpeg');
imwrite(H9.cdata,[num2str(init1) '_' char(list_iodu(composto)) '_div_map_' 'method_asym' '.jpg'] , 'jpeg');
saveas(h8,[num2str(init1)  '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method_asym' '.fig']);
saveas(h9,[num2str(init1)  '_' char(list_iodu(composto)) '_div_map_' 'method_asym''.fig']);
H11=getframe(h11);
imwrite(H11.cdata,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method_asym' '.jpg'] , 'jpeg');    
saveas(h11,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method_asym' '.fig']);
imwrite(H7.cdata,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method_asym' '.jpg'] , 'jpeg'); 
saveas(h7,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method_asym' '.fig']);

%--------------------------------------------------------------------------
% positive method (#1)
h7=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2)  ,i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST1_2(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST1_2(:,:,i)>0).*0.7); 
    title(['deltaST at ' num2str(st_ppm) ' ppm. Slice:' num2str(i)])
end

h8=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST2_2(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST2_2(:,:,i)>0).*0.7);
    title(['deltaST at ' num2str(st_ppm2) ' ppm. Slice:' num2str(i)])
end

h9=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapSTratio2(:,:,i), [0 3]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapSTratio2(:,:,i)>0).*0.7);
    title( ['Delta ratio. Slice:' num2str(i)])
end

h11=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(pH2(:,:,i), [5 8]);
    set(gcf(),'Colormap',pHmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (pH2(:,:,i)>0).*0.7);
    title(['Image pH. Slice:' num2str(i)])
end

H7=getframe(h7);
H8=getframe(h8);
H9=getframe(h9);
imwrite(H8.cdata,[num2str(init1) '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method1' '.jpg'] , 'jpeg');
imwrite(H9.cdata,[num2str(init1) '_' char(list_iodu(composto)) '_div_map_' 'method1' '.jpg'] , 'jpeg');
saveas(h8,[num2str(init1)  '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method1' '.fig']);
saveas(h9,[num2str(init1)  '_' char(list_iodu(composto)) '_div_map_' 'method1''.fig']);
H11=getframe(h11);
imwrite(H11.cdata,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method1' '.jpg'] , 'jpeg');    
saveas(h11,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method1' '.fig']);
imwrite(H7.cdata,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method1' '.jpg'] , 'jpeg'); 
saveas(h7,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method1' '.fig']);

%--------------------------------------------------------------------------
% linear method (#2)
h7=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2)  ,i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST1_3(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST1_3(:,:,i)>0).*0.7); 
    title(['deltaST at ' num2str(st_ppm) ' ppm. Slice:' num2str(i)])
end

h8=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST2_3(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST2_3(:,:,i)>0).*0.7);
    title(['deltaST at ' num2str(st_ppm2) ' ppm. Slice:' num2str(i)])
end

h9=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapSTratio3(:,:,i), [0 3]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapSTratio3(:,:,i)>0).*0.7);
    title( ['Delta ratio. Slice:' num2str(i)])
end

h11=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(pH3(:,:,i), [5 8]);
    set(gcf(),'Colormap',pHmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (pH3(:,:,i)>0).*0.7);
    title(['Image pH. Slice:' num2str(i)])
end

H7=getframe(h7);
H8=getframe(h8);
H9=getframe(h9);
imwrite(H8.cdata,[num2str(init1) '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method2' '.jpg'] , 'jpeg');
imwrite(H9.cdata,[num2str(init1) '_' char(list_iodu(composto)) '_div_map_' 'method2' '.jpg'] , 'jpeg');
saveas(h8,[num2str(init1)  '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method2' '.fig']);
saveas(h9,[num2str(init1)  '_' char(list_iodu(composto)) '_div_map_' 'method2''.fig']);
H11=getframe(h11);
imwrite(H11.cdata,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method2' '.jpg'] , 'jpeg');    
saveas(h11,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method2' '.fig']);
imwrite(H7.cdata,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method2' '.jpg'] , 'jpeg'); 
saveas(h7,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method2' '.fig']);

%--------------------------------------------------------------------------
% lorentzian method (#3)
h7=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2)  ,i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST1_4(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST1_4(:,:,i)>0).*0.7); 
    title(['deltaST at ' num2str(st_ppm) ' ppm. Slice:' num2str(i)])
end

h8=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapST2_4(:,:,i), [0 0.2]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapST2_4(:,:,i)>0).*0.7);
    title(['deltaST at ' num2str(st_ppm2) ' ppm. Slice:' num2str(i)])
end

h9=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapSTratio4(:,:,i), [0 3]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapSTratio4(:,:,i)>0).*0.7);
    title( ['Delta ratio. Slice:' num2str(i)])
end

h11=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(pH4(:,:,i), [5 8]);
    set(gcf(),'Colormap',pHmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (pH4(:,:,i)>0).*0.7);
    title(['Image pH. Slice:' num2str(i)])
end

H7=getframe(h7);
H8=getframe(h8);
H9=getframe(h9);
imwrite(H8.cdata,[num2str(init1) '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method3' '.jpg'] , 'jpeg');
imwrite(H9.cdata,[num2str(init1) '_' char(list_iodu(composto)) '_div_map_' 'method3' '.jpg'] , 'jpeg');
saveas(h8,[num2str(init1)  '_deltaST_map_' num2str(st_ppm2) '_ppm' '_method3' '.fig']);
saveas(h9,[num2str(init1)  '_' char(list_iodu(composto)) '_div_map_' 'method3''.fig']);
H11=getframe(h11);
imwrite(H11.cdata,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method3' '.jpg'] , 'jpeg');    
saveas(h11,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method3' '.fig']);
imwrite(H7.cdata,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method3' '.jpg'] , 'jpeg'); 
saveas(h7,[num2str(init1)  '_deltaST_map_' num2str(st_ppm) '_ppm' '_method3' '.fig']);

%--------------------------------------------------------------------------
% interpolation method (#4)
h9=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(mapSTratio5(:,:,i), [0 3]);
    set(gcf(),'Colormap',mycmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (mapSTratio5(:,:,i)>0).*0.7);
    title( ['Delta ratio. Slice:' num2str(i)])
end

h11=figure('Position', get(0, 'Screensize')); 
for i=1:size_zz
    subplot(2,ceil(size_zz/2),i)
    imshow(image_r(:,:,i), []);
    title(['slice' num2str(i)])
    colormap gray
    freezeColors
    hold on
    h=imshow(pH5(:,:,i), [5 8]);
    set(gcf(),'Colormap',pHmap) 
    c=colorbar();
    set(c,'FontSize',18, 'FontWeight', 'bold')
    set(h, 'AlphaData', (pH5(:,:,i)>0).*0.7);
    title(['Image pH. Slice:' num2str(i)])
end

H9=getframe(h9);
imwrite(H9.cdata,[num2str(init1) '_' char(list_iodu(composto)) '_div_map_' 'method4' '.jpg'] , 'jpeg');
saveas(h9,[num2str(init1)  '_' char(list_iodu(composto)) '_div_map_' 'method4''.fig']);
H11=getframe(h11);
imwrite(H11.cdata,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method4' '.jpg'] , 'jpeg');    
saveas(h11,[ num2str(init1)  '_' char(list_iodu(composto)) '_pH_map_' 'method4' '.fig']);
%% Params
tot_pixel_slice = sum(mask_back(:));

pH_slice1 = pH1(pH1>0 & ~isnan(pH1));
enh_pixel_slice1 = length(pH_slice1);
mean_pH_slice_1 = round(mean(pH_slice1),2);
std_pH_slice_1 = round(std(pH_slice1),2);
FP_slice1 = (enh_pixel_slice1/tot_pixel_slice)*100;

pH_slice2 = pH2(pH2>0 & ~isnan(pH2));
enh_pixel_slice2 = length(pH_slice2);
mean_pH_slice_2 = round(mean(pH_slice2),2);
std_pH_slice_2 = round(std(pH_slice2),2);
FP_slice2 = (enh_pixel_slice2/tot_pixel_slice)*100;

pH_slice3 = pH3(pH3>0 & ~isnan(pH3));
enh_pixel_slice3 = length(pH_slice3);
mean_pH_slice_3 = round(mean(pH_slice3),2);
std_pH_slice_3 = round(std(pH_slice3),2);
FP_slice3 = (enh_pixel_slice3/tot_pixel_slice)*100;

pH_slice4 = pH4(pH4>0 & ~isnan(pH4));
enh_pixel_slice4 = length(pH_slice4);
mean_pH_slice_4 = round(mean(pH_slice4),2);
std_pH_slice_4 = round(std(pH_slice4),2);
FP_slice4 = (enh_pixel_slice4/tot_pixel_slice)*100;

pH_slice5 = pH5(pH5>0 & ~isnan(pH5));
enh_pixel_slice5 = length(pH_slice5);
mean_pH_slice_5 = round(mean(pH_slice5),2);
std_pH_slice_5 = round(std(pH_slice5),2);
FP_slice5 = (enh_pixel_slice5/tot_pixel_slice)*100;

close all
save('FFvariables.mat')