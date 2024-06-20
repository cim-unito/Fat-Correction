function [ img , hdr ] = readBruker(pathScan)



% function [ img , hdr ] = readBrucker(pathScan)
%
% This function reads image data in the Bruker format and return it into
% an array.
%
% Daisy Villano
%  Matlab 2015b

%+ pV 5.1
%+ pV 6
%+ pV 360
%+ RARE morfo
%+ CEST 1slice
%+ CEST n-slices
%+ VFA mapping
%+ VTR mapping
%- DCE
%- DWI
%- pV 5 3T Ivrea
rep=1;
freq_size=1;
freq_list=1;
hdr.TR_list=1;
hdr.TRExps=1;
hdr.VFAExps=1;
hdr.NRepetitions=1;
if nargin < 1
    pathScan  = uigetdir('Select the directory containing the Scan');
end

% tic
fprintf('Reading from %s',pathScan)

if exist([pathScan,'/pdata/1/2dseq'],'file')
    imagefile = [pathScan,'/pdata/1/2dseq'];
else
    error('2dseq does not exist!')
end

if exist([pathScan,'/method'],'file')
    methodfile = [pathScan,'/method'];
else
    error('method does not exist!')
end

if exist([pathScan,'/acqp'],'file')
    acqpfile = [pathScan,'/acqp'];
else
    error('acqp does not exist!')
end

if exist([pathScan,'/pdata/1/reco'],'file')
    recofile = [pathScan,'/pdata/1/reco'];
else
    error('reco does not exist!')
end

if exist([pathScan,'/pdata/1/visu_pars'],'file')
    visu_parsfile = [pathScan,'/pdata/1/visu_pars'];
else
    error('visu_pars does not exist!')
end

fprintf('.')
fp = fopen(methodfile,'r');
method = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(method,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch method(idx_sta(idx)+3:idx_end(idx)-1)
        case 'Method' %string: Measuring method
            hdr.Method = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'PVM_EchoTime' %(ms)
            hdr.EchoTime = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_NAverages' %integer: Number of times the signal is accumulated prior to the storage on disk and the reconstruction.
            hdr.NAverages = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_NRepetitions' %integer: Number of repetitions (executions) of the experiment.
            hdr.NRepetitions = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_RepetitionTime' %(ms)
            hdr.RepetitionTime = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_RareFactor' %integer: Number of times the signal is accumulated prior to the storage on disk and the reconstruction.
            hdr.RareFactor = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SliceThick' %(mm)
            hdr.SliceThick = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_FatSupOnOff'
            hdr.FatSupOnOff = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'PVM_ObjOrderScheme' %string: Selects the order in which the slices are excited in a multislice experiment. (Interlaced, Reverse_sequential, Sequential)
            hdr.ObjOrderScheme = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'PVM_Fov' %(ms)
            hdr.Fov =  sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+8),'%f');
        case 'PVM_SPackArrSliceOrient' %array of strings: General orientation of each slice package (axial, sagittal, coronal)
            hdr.SPackArrSliceOrient = sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+8),'%s');
        case 'PVM_ScanTimeStr' %integer
            method1=method(idx_end(idx)+1:end);
            idx1=find(method1=='<',1);
            idx2=find(method1=='>',1);
            hdr.ScanTime = char(method1(idx1:idx2));
            clear method1
            
            %% CEST variables
            % pV5.1
        case 'PVM_MagTransPower' % (uT)
            hdr.CESTPulsePower = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_MagTransPulse' % (ms)
            hdr.CESTPulseLength = str2double(method(idx_end(idx)+2:idx_end(idx)+5));%+find_newline(method,idx_end(idx))-1));
        case 'PVM_MtP0' % (us)
            hdr.MTPulseLength = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_MTPulseLength2' %(ms)
            hdr.CESTPulseLength2 = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
            % pV6
        case 'PVM_MagTransPulse1' % (ms)
            hdr.MagTransPulseLength = str2double(method(idx_end(idx)+2:idx_end(idx)+5));%+find_newline(method,idx_end(idx))-1));
            % pV360
        case 'PVM_SatTransPulseLength2' %(ms)
            hdr.CESTPulseLength2 = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SatTransPulseAmpl_uT' % (uT)
            hdr.CESTPulsePower = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SatTransPulse' % (ms)
            hdr.CESTPulseLength = str2double(method(idx_end(idx)+2:idx_end(idx)+5));
        case 'PVM_SatTransRepetitions'
            freq_size = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SatTransFL'
            method1=method(idx_end(idx)+1+find_newline(method,idx_end(idx)):end);
            idx1=find(method1=='#',1);
            freq_list = sscanf(method(idx_end(idx)+1+find_newline(method,idx_end(idx)):idx_end(idx)+find_newline(method,idx_end(idx))+idx1-1),'%f');
        case 'PVM_MagTransFL'
            hdr.freq=sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+8),'%f');
            
            %% VFA mapping
            %% VTR mapping
        case 'NumT1Exps' %integer: Number of repetitions (executions) of the experiment.
            hdr.TRExps = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
    end
end

fprintf('.')
fp = fopen( acqpfile, 'r' );
acqp = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(acqp,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch acqp(idx_sta(idx)+3:idx_end(idx)-1)
        
        case 'NSLICES' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            hdr.nslices = str2double(acqp(idx_end(idx)+1:idx_end(idx)+find_newline(acqp,idx_end(idx))-1));
        case 'NI' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            hdr.NI = str2double(acqp(idx_end(idx)+1:idx_end(idx)+find_newline(acqp,idx_end(idx))-1));
        case 'BF3' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            hdr.B0 = str2double(acqp(idx_end(idx)+1:idx_end(idx)+find_newline(acqp,idx_end(idx))-1));
            %         case 'ACQ_protocol_name' %integer
            %             acqp1=acqp(idx_end(idx):end);
            %             idx1=find(acqp1=='<',1);
            %             idx2=find(acqp1=='>',1);
            %             hdr.protocol_name = char(acqp1(idx1:idx2));
            %             clear acqp1
        case 'ACQ_time' %integer
            acqp1=acqp(idx_end(idx):end);
            idx1=find(acqp1=='<',1);
            idx2=find(acqp1=='>',1);
            hdr.ACQTime = char(acqp1(idx1:idx2));
            clear acqp1
        case 'ACQ_institution' %integer
            acqp1=acqp(idx_end(idx):end);
            idx1=find(acqp1=='<',1);
            idx2=find(acqp1=='>',1);
            hdr.institution = char(acqp1(idx1:idx2));
            clear acqp1
        case 'ACQ_station' %integer
            acqp1=acqp(idx_end(idx):end);
            idx1=find(acqp1=='<',1);
            idx2=find(acqp1=='>',1);
            hdr.scanner = char(acqp1(idx1:idx2));
            clear acqp1
        case 'ACQ_sw_version' %integer
            acqp1=acqp(idx_end(idx):end);
            idx1=find(acqp1=='<',1);
            idx2=find(acqp1=='>',1);
            hdr.pv = char(acqp1(idx1:idx2));
            clear acqp1
            %% CEST variables
            % pV5.1 pV6
        case 'ACQ_O2_list' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            acqp1=acqp(idx_end(idx)+1+find_newline(acqp,idx_end(idx)):end);
            idx1=find(acqp1=='#',1);
            o2_list = sscanf(acqp(idx_end(idx)+1+find_newline(acqp,idx_end(idx)):idx_end(idx)+find_newline(acqp,idx_end(idx))+idx1-1),'%f');
            
        case 'ACQ_O2_list_size' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            o2_list_size = str2double(acqp(idx_end(idx)+1:idx_end(idx)+find_newline(acqp,idx_end(idx))-1));
            
            %% VFA mapping
            %% VTR mapping
        case 'ACQ_repetition_time' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            acqp1=acqp(idx_end(idx)+1+find_newline(acqp,idx_end(idx)):end);
            idx1=find(acqp1=='#',1);
            hdr.TR_list = sscanf(acqp(idx_end(idx)+1+find_newline(acqp,idx_end(idx)):idx_end(idx)+find_newline(acqp,idx_end(idx))+idx1-1),'%f');
            
    end
end

fprintf('.')
fp = fopen(recofile,'r');
reco = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(reco,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch reco(idx_sta(idx)+3:idx_end(idx)-1)
        
        case 'RECO_ft_size' %string. either SumOfSquares, AddImages, or ShuffleImages
            i1=find_newline(reco,idx_end(idx)):idx_end(idx);
            hdr.size_xy = sscanf(reco(idx_end(idx)+find_newline(reco,idx_end(idx)):idx_end(idx)+find_newline(reco,idx_end(idx))+10),'%f');
        case 'RECO_wordtype' %the parameter describes the structure of the input data. When RecoNumInputChan > 1, reconstruction assumes, that the raw data file consists of RecoNumInputChan blocks of size RECO_inp_size[0] forming the first dimension of the data file.
            hdr.wordtype = char(reco(idx_end(idx)+1:idx_end(idx)+find_newline(reco,idx_end(idx))-1));
        case 'RECO_map_slope' % array of numbers
            hdr.slope = sscanf(reco(idx_end(idx)+1+find_newline(reco,idx_end(idx)):idx_end(idx)+find_newline(reco,idx_end(idx))+20),'%f', [1 1]);
            if isempty(hdr.slope)
                reco1=reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:end);
                idx1=find(reco1=='(');
                idx2=find(reco1==')');
                hdr.slope=str2double(reco1(idx1+1:idx2-1));
            end
    end
end

fprintf('.')
fp = fopen(visu_parsfile,'r');
visu = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(visu,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch visu(idx_sta(idx)+3:idx_end(idx)-1)
        case 'VisuAcqFlipAngle'
            hdr.VFA = str2double(visu(idx_end(idx)+1:idx_end(idx)+find_newline(visu,idx_end(idx))-1));
        case 'VisuAcquisitionProtocol'
            visu1=visu(idx_end(idx):end);
            idx1=find(visu1=='<',1);
            idx2=find(visu1=='>',1);
            hdr.protocol_name=char(visu1(idx1:idx2));
            
    end
end
num_ppm=max(o2_list_size,freq_size);
if num_ppm==o2_list_size
    freq_list=o2_list;
end
hdr.freq_list = freq_list./hdr.B0;
hdr.num_ppm=num_ppm;
switch hdr.wordtype
    case '_32BIT_SGN_INT'
        precision = 'int32';
    case '_16BIT_SGN_INT'
        precision = 'int16';
    case '_8BIT_UNSGN_INT'
        precision = 'uint8';
    case '_32BIT_FLOAT'
        precision = 'single';
end
if hdr.TRExps<length(hdr.TR_list)
    hdr.TRExps=length(hdr.TR_list);
end
if numel(hdr.size_xy)>2
    rep=hdr.size_xy(3);
end
ni=hdr.NI/hdr.nslices;
hdr.ntime=max([hdr.NRepetitions,hdr.num_ppm,hdr.TRExps,hdr.VFAExps,rep,ni]);

fid=fopen(imagefile);

for t=1:hdr.ntime
    
    for z=1:hdr.nslices
        for i=1:hdr.size_xy(2)
            for j=1:hdr.size_xy(1)
                image(i,j,z,t)=fread(fid,1,precision);%,'ieee-le');
            end
        end
    end
end

img=squeeze(image);


fprintf('%s\n','done!')
end

