clear;
close all;
tStart=tic;

%% Define variables

pathname='/Users/christo/Desktop/test/raw';

makeWF=1; % save WF image
makeWR=1; % process Wiener reconstruction
makeHF=1; % save HiFi reconstruction

splitSlices=0; % split slices for multiplane files
paramFirst=0; % only estimate parameters for first image
paramSlice=1; % only estimate parameters for first slice of multiplane files
normImages=0; % normalize image intensities from 0 to 255 at max intensity before saving
savePRS=0; % save in subfolders  

im3D=0; % is 3D-SIM
lambda=525; % emission wavelenght in nm
pixelsize=65; % image pixel size in nm
NA=1.49; % objective numerical aperture

attStrength=0.9;
a=1; % damping factor：β
attFWHM=1.0;

if im3D==0
    nrDirs=3; % number of angles
    nrPhases=3; % number of phases
else
    nrDirs=3; % number of angles
    nrPhases=5; % number of phases  
end

w1=1.2;  % Initial optimization Wiener constant：[0.9-2.5]
w2=0.1;
ApoFWHM=0; % for 0 automatic setting

%% Preparation

if im3D==0
    param.phaOff=0;   
    param.fac=ones(1,2);    
    param.nrBands=2;
else
    param.phaOff=0;   
    param.fac=ones(1,3);    
    param.nrBands=3;
end

param.pathname=pathname;
param.nrDirs=nrDirs;
param.nrPhases=nrPhases;

% create output folders
ipr=strfind(pathname,filesep);
parentpath=pathname(1:ipr(end)-1);
dirname=pathname(ipr(end)+1:strlength(pathname));

if makeWF==1
    savedirWF=[parentpath,filesep,dirname,'_WF'];
    param.savedirWF=savedirWF;
    if ~exist(savedirWF,'dir')
        mkdir(savedirWF);
    end
end

if makeWR==1
    savedirWR=[parentpath,filesep,dirname,'_WR'];
    param.savedirWR=savedirWR;
    if ~exist(savedirWR,'dir')
        mkdir(savedirWR);
    end
end

if makeHF==1
    savedirHF=[parentpath,filesep,dirname,'_HF'];
    param.savedirHF=savedirHF;
    if ~exist(savedirHF,'dir')
        mkdir(savedirHF);
    end
end


%% Loop on files in folder

filelist=dir([pathname,filesep,'*.tif']);
filecount=length(filelist);
disp(['Batch HiFi-SIM started in folder:',pathname]);

for fileid=1:filecount
    
    %% Load raw SIM data
    filename=filelist(fileid).name;
    param.filename=filename;    
    disp(['processing file ',num2str(fileid),'/',num2str(filecount),': ',param.filename]);

    N=param.nrDirs*param.nrPhases;
    info=imfinfo(fullfile(param.pathname, param.filename));

    inSlices=numel(info);
    inWidth=info(1).Width;
    inHeight=info(1).Height;

    if inWidth > 2049       % Case of image stacks (order: phases, angles)
        isMosaic=1;
        inImages=inSlices;
    else                    % Case of mosaic image (lines:phases, columns:angles)
        isMosaic=0; 
        inImages=inSlices/N;
    end

    % Loop on slices in input image
    for sliceid=1:inImages 
    
        disp([' processing slice ',num2str(sliceid),'/',num2str(inImages)]);
        
        % test if output image will be created (deleting existing file) or appended to existing stack
        if sliceid==1
            app=0;
        else
            app=1;
        end
     
        % Case of image stacks (order: phases, angles)
        if isMosaic==0
            NPixel=max(inWidth,inHeight);
            Iraw0=zeros(inHeight,inWidth,N);
            Iraw=zeros(NPixel,NPixel,N); 
            for j=0:N-1
                Iraw0(:,:,j)=double(imread(fullfile(param.pathname, param.filename),sliceid+j));
                if inWidth == inHeight || inWidth > inHeight
                    Iraw(1:inHeight,:,j)=Iraw0(1:inHeight,:,j);
                else
                    Iraw(:,1:inWidth,j)=Iraw0(:,1:inWidth,j);
                end 
            end
        % Case of mosaic image (lines:phases, columns:angles)
        else
            Iraw0=zeros(inHeight,inWidth,1);
            Iraw0(:,:,1)=double(imread(fullfile(param.pathname, param.filename),sliceid));  
            ImHeight=inHeight/param.nrDirs;
            ImWidth=inWidth/param.nrPhases;
            NPixel=max(ImWidth,ImHeight);
            Iraw=zeros(ImHeight,ImWidth,N);        
            for j=1:param.nrDirs
                ShiftY=ImHeight*(j-1);
                for k=1:param.nrPhases
                    ShiftX=ImWidth*(k-1);
                    Slicepos=k+(j-1)*param.nrDirs;
                    Iraw(1:ImHeight,1:ImWidth,Slicepos)=Iraw0(ShiftY+1:ShiftY+ImHeight,ShiftX+1:ShiftX+ImWidth,1); 
                end
            end
        end
        
        Format=info.Format;
        param.Format=Format;
        param.imgSize=NPixel;
        if numel(info)==N
            param.Size1 = info(1).Height; 
            param.Size2 = info(1).Width;
        else
            param.Size1 = ImHeight; 
            param.Size2 = ImWidth; 
        end
        param.Iraw = Iraw;
        
        
        %% Processing part
        
    
        %% Generate approximate OTF/PSF
        Iraw=param.Iraw;
        NPixel=size(Iraw,1);  
        param.micronsPerPixel=pixelsize*10^(-3);
        param.cyclesPerMicron=1/(NPixel*param.micronsPerPixel);
        param.NA=NA;
        param.lambda=lambda;
        param.attStrength=0;
        param.OtfProvider=SimOtfProvider(param,param.NA,param.lambda,1);
        
        psf=abs(otf2psf((param.OtfProvider.otf)));
       
    
        %% Preprocessing
        Temp=importImages(Iraw);  
        IIraw=deconvlucy(Temp,psf,5);
        for I=1:N
            IIrawFFT(:,:,I)=FFT2D(IIraw(:,:,I),false);
        end
        param.IIrawFFT=IIrawFFT;
    
        WF=zeros(NPixel,NPixel);
        Tdeconv=zeros(NPixel,NPixel);
        WFdeconv=zeros(NPixel,NPixel,param.nrDirs);
        WFdeconvFFT=zeros(NPixel,NPixel,param.nrDirs);
        for i=1:param.nrDirs
            for j=1:param.nrPhases
                Tdeconv=Tdeconv+IIraw(:,:,(i-1)*param.nrDirs+j);
                WF(:,:)=WF(:,:)+Iraw(:,:,(i-1)*param.nrDirs+j);
            end
            WFdeconv(:,:,i)=Tdeconv/param.nrPhases;
            WFdeconvFFT(:,:,i)=FFT2D(WFdeconv(:,:,i),false);
        end
        WF=WF/N;
        WF=importImages(WF); 
      
    
        %% Wide-field
        Size1=2*param.Size1;
        Size2=2*param.Size2;
        WF2=zeros(Size1,Size2);
        Temp=zeros(2*size(WF,1),2*size(WF,2));
        fftWF=fftshift(fft2(WF));
        Temp(size(WF,1)/2+1:size(WF,1)/2+size(WF,1),size(WF,2)/2+1:size(WF,2)/2+size(WF,2))=fftWF;
        Temp=abs(ifft2(Temp));
        if normImages==1
            Temp=255*Temp/max(max(Temp));
        end
        WF2(1:Size1,1:Size2)=Temp(1:Size1,1:Size2);
        WF2=importImages(WF2);
        param.WF2=WF2;
     
    
        %% Parameter Estimation (with conditions to be done only on first file or first slice)
        if (paramFirst==1 && fileid == 1) || paramFirst==0
            if (paramSlice == 1 && sliceid == 1) || paramSlice==0
        
            % Kc region MASK
            cnt=[NPixel/2+1,NPixel/2+1];
            param.cutoff=1000/(0.5*param.lambda/param.NA);
            [x,y]=meshgrid(1:NPixel,1:NPixel);
            rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
            Mask=double(rad<=1.0*(param.cutoff/param.cyclesPerMicron+1));
            NotchFilter0=getotfAtt(NPixel,param.OtfProvider.cyclesPerMicron,0.5*param.cutoff,0,0);
            NotchFilter=NotchFilter0.*Mask;
            Mask2=double(rad<=1.10*(param.cutoff/param.cyclesPerMicron+1));
            NotchFilter2=NotchFilter0.*Mask2;
            
            CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),param.nrDirs);
            k0=zeros(1,param.nrDirs);
            
            for I=1:param.nrDirs
                lb=2;
                if param.nrBands==2
                    hb=2;
                    fb=lb;
                elseif param.nrBands==3
                    hb=4;
                    fb=hb;
                end
                 
                param.phaOff=0;    
                param.fac=ones(1,param.nrBands); 
                separateII=separateBands(IIrawFFT(:,:,(I-1)*param.nrPhases+1:I*param.nrPhases),param.phaOff,param.nrBands,param.fac);
                
                SeparateII{1,I}=separateII;
                
                if param.nrBands==2
                    c0=separateII(:,:,1);
                    c2=separateII(:,:,lb);
                   
                    c0=(c0./(max(max(abs(c0)))));
                    c2=(c2./(max(max(abs(c2)))));
                  
                    c0=c0.*NotchFilter;
                    c2=c2.*NotchFilter;
                    
                    c0=FFT2D(c0,false);
                    c2=FFT2D(c2,false);
                    c2=c2.*conj(c0); 
                    c2=c2./max(max(c2));              
            
                    vec=fftshift(FFT2D(c2,true));
                elseif param.nrBands==3
                    c0=separateII(:,:,1);
                    c3=separateII(:,:,hb);
                    
                    c0=c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter;
                    c3=c3./(max(max(abs(separateII(:,:,hb))))).*NotchFilter;
            
                    c0=FFT2D(c0,false);
                    c3=FFT2D(c3,false);
                    c3=c3.*conj(c0);
                    c3=c3./max(max(c3)); 
            
                    vec=fftshift(FFT2D(c3,true));
                end
                CrossCorrelation(:,:,I)=vec;
                %% 
            %       temp=vec.*NotchFilter; 
                temp=vec.*NotchFilter2;
                temp=log(1+abs(temp));
                temp=temp./max(max(temp));
            
                [yPos,xPos]=find(temp==max(max(temp)));
                peak.xPos=xPos(1);
                peak.yPos=yPos(1);
                k0(I)=sqrt((peak.xPos-cnt(1))^2+(peak.yPos-cnt(2))^2);
            end
            
            Flag=0;
            if param.nrDirs>2           % For very few special cases
                if max(k0)-min(k0)>8
                    Flag=1;
                    Kobject=min(k0); 
            %             Kobject=208; 
                    Mask1=rad>=(Kobject+1);
                    Mask2=rad<=(Kobject-1);
                end
            end
            
            for I=1:param.nrDirs
                vec=CrossCorrelation(:,:,I);
                if Flag==1
                    vec(Mask1)=0;
                    vec(Mask2)=0;
                end
                temp=vec.*NotchFilter2;
                temp=log(1+abs(temp));
                temp=temp./max(max(temp));
            
                [yPos,xPos]=find(temp==max(max(temp)));
                peak.xPos=xPos(1);
                peak.yPos=yPos(1);
            
                cntrl=zeros(10,30);
                overlap=0.15;
                step=2.5;
                bn1=(param.nrBands-1)*2;
                kx=(peak.xPos-cnt(2));
                ky=(peak.yPos-cnt(1));
                
                separateII=SeparateII{1,I};
                [peak,cntrl]=fitPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,fb)/(max(max(abs(separateII(:,:,fb))))),1,bn1,param.OtfProvider,-kx,-ky,overlap,step,cntrl);
                
                if lb~=hb
                    if param.nrBands==2
                        peak.kx=peak.kx*2;
                        peak.ky=peak.ky*2;
                    end
                    
                    p1=getPeak(separateII(:,:,1),separateII(:,:,lb),0,1,param.OtfProvider,peak.kx/2,peak.ky/2,overlap);
                    p2=getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,hb)/(max(max(abs(separateII(:,:,1))))),0,2,param.OtfProvider,peak.kx,peak.ky,overlap);
                    
                    param.Dir(I).px=-peak.kx/2;    
                    param.Dir(I).py=-peak.ky/2;
                    param.Dir(I).phaOff=-phase(p1);
                    Temp_m1=abs(p1);
                    Temp_m2=abs(p2);
                    
                    Temp_m1(Temp_m1>1.0)=1;
                    Temp_m2(Temp_m2>1.0)=1.0;
                    param.Dir(I).modul(1)=Temp_m1;
                    param.Dir(I).modul(2)=Temp_m2;
                end
                if lb==hb
                    p1=getPeak(separateII(:,:,1),separateII(:,:,lb),1,lb,param.OtfProvider,peak.kx,peak.ky,overlap);
                    param.Dir(I).px=-peak.kx;
                    param.Dir(I).py=-peak.ky;
                    param.Dir(I).phaOff=-phase(p1);
                    Temp_m=abs(p1);
                    Temp_m(Temp_m>1.0)=1.0;
                    param.Dir(I).modul=Temp_m;
                end
                K0(I)=sqrt((param.Dir(I).px)^2+(param.Dir(I).py)^2);
                %% 
            %             fittedPeak=vec;
            %             for x=1:30;
            %                 for y=1:10
            %                     fittedPeak(y,x)=cntrl(y,x);
            %                 end
            %             end
                    
            %             figure,imshow(temp,[]);
            %             colormap('hot');
            %             hold on;
            %             
            %             if lb~=hb
            %                 f=(param.nrBands-1)/2;
            %             else
            %                 f=1;
            %             end
            %             t=0:0.1:2.1*pi;
            %             x=cnt(2)-peak.kx+10*sin(t);
            %             y=cnt(1)-peak.ky+10*cos(t);
            %             plot(x,y,'-w','LineWidth',2);
            % %             title(['Orientation',num2str(I)]);
            %             %         plotbrowser('on');
            % %         end
            %         %     plotbrowser('off');
            end
            
            %%
            SIMparam=zeros(3,5);
            if param.nrPhases==3
                for i=1:param.nrDirs
                    SIMparam(i,1)=atan(param.Dir(i).py/param.Dir(i).px)*180/pi;
                    SIMparam(i,2)=sqrt((param.Dir(i).px)^2+(param.Dir(i).py)^2);
                    SIMparam(i,3)=param.Dir(i).phaOff;
                    SIMparam(i,4)=param.Dir(i).modul;
                    if param.Dir(i).modul<0.35
                        param.Dir(i).modul=0.7;
                    end
                end
            elseif param.nrPhases==5
                for i=1:param.nrDirs
                    SIMparam(i,1)=atan(param.Dir(i).py/param.Dir(i).px)*180/pi;
                    SIMparam(i,2)=sqrt((param.Dir(i).px)^2+(param.Dir(i).py)^2)*2;
                    SIMparam(i,3)=param.Dir(i).phaOff;
                    SIMparam(i,4)=param.Dir(i).modul(1);
                    SIMparam(i,5)=param.Dir(i).modul(2);
                    if param.Dir(i).modul(1)<0.35
                        param.Dir(i).modul(1)=0.7;
                    end
                    if param.Dir(i).modul(2)<0.35
                        param.Dir(i).modul(2)=0.7;
                    end
                end
            end        
            param.K0=K0;
            disp('  illumination parameters estimation done');
      
            end
        end
        
    
        %% Reconstruction
        
        IIrawFFT=param.IIrawFFT;
        Iraw=param.Iraw;
        siz=size(Iraw(:,:,1));
        w=siz(2);
        h=siz(1);
        
        param.attStrength=attStrength;
        param.a=a;
        param.attFWHM=attFWHM;
        param.OtfProvider=SimOtfProvider(param,param.NA,param.lambda,param.a);
        
               
        %% HiFi-SIM：Spectrum optimization
        
        fftDirectlyCombined=zeros(h*2,w*2);
        for I=1:param.nrDirs
            par=param.Dir(I);
            param.fac(2:param.nrBands)=param.Dir(I).modul(1:param.nrBands-1);   
            param.fac(2:param.nrBands)=param.Dir(I).modul(1:param.nrBands-1);
            separate=separateBands(IIrawFFT(:,:,(I-1)*param.nrPhases+1:I*param.nrPhases),par.phaOff,param.nrBands,param.fac);
            
            shifted=zeros(2*h,2*w,param.nrPhases);
            shifted(:,:,1)=placeFreq(separate(:,:,1));
            
            for b=2:param.nrBands
                pos=b*2-2;
                neg=b*2-1;
                shifted(:,:,pos)=placeFreq(separate(:,:,pos));
                shifted(:,:,neg)=placeFreq(separate(:,:,neg));
                
                shifted(:,:,pos)=NfourierShift(shifted(:,:,pos),-(b-1)*par.px,-(b-1)*par.py);
                shifted(:,:,neg)=NfourierShift(shifted(:,:,neg),(b-1)*par.px,(b-1)*par.py);
            end
                shifted(:,:,1)=applyOtf(shifted(:,:,1),param.OtfProvider,1,0,0,1,0);
                for b=2:param.nrBands
                    pos=b*2-2;
                    neg=b*2-1;
                    shifted(:,:,pos)=applyOtf(shifted(:,:,pos),param.OtfProvider,b,-(b-1)*par.px,-(b-1)*par.py,1,0);
                    shifted(:,:,neg)=applyOtf(shifted(:,:,neg),param.OtfProvider,b,(b-1)*par.px,(b-1)*par.py,1,0);
                end
            for J=1:param.nrBands*2-1
                fftDirectlyCombined=fftDirectlyCombined+shifted(:,:,J);
            end
        end
       
        % Temp1=real(ifft2(fftshift((fftDirectlyCombined))));
        % Temp1(Temp1<0)=0;
        % MIJ.createImage(Temp1);
        disp('  spectrum optimization done');
    
        % w1=1.2;  % Initial optimization Wiener constant：[0.9-2.5]
        % w2=0.1;    
        param.cutoff=1000/(0.5*param.lambda/param.NA);                        
        param.sampleLateral=ceil(param.cutoff/param.cyclesPerMicron)+1;  
        K0=param.K0;
        K=max([ceil(K0)]);
        if param.nrBands==2
            cutoff=floor(1*K)/param.sampleLateral+1.0;
            R=K;
        elseif	param.nrBands==3
            cutoff=floor(2*K)/param.sampleLateral+1.0;
            R=2*K;
        end 
        otfHiFi=zeros(2*h,2*w);
        otfHiFi=writeApoVector(otfHiFi,param.OtfProvider,cutoff);      % Ideal OTF
        Mask=zeros(2*h,2*w);
        Mask(otfHiFi~=0)=1;
        
    
        %% Traditional Wiener-SIM
        if makeWR==1
          
            if size(Iraw,3)==9
                wFilter0=WienerFilterWiener_2D(param);
            else
                wFilter0=WienerFilterWiener_3D(param);
            end
            Wk0=otfHiFi./(wFilter0.wDenom+w2^2);
            fftWiener=real(ifft2(fftshift((fftDirectlyCombined.*Wk0.*Mask))));
            Temp=fftWiener;
            Temp(Temp<0)=0;
            if normImages==1
                Wiener=255*Temp/max(max(Temp));
            end
            Wiener=Temp;
            % MIJ.createImage(Wiener);                  % The reconstruction results are displayed in the imageJ window
            % figure, imshow(Wiener,[]);    % The reconstruction results are displayed in the matlab window
            % colormap('hot');
            param.WR=Wiener;
            disp('  Wiener reconstruction done');
            
            if savePRS==0
                if splitSlices==0
                    outnameWR=strrep(param.filename, '.tif', '_WR.tif');
                else
                    outnameWR=strrep(param.filename, '.tif', ['_s',num2str(sliceid,'%03d'),'_WR.tif']);
                end
                imgsave32seq(Wiener, fullfile(param.savedirWR,outnameWR), app);

            else
                outnameWR=param.filename;
                in = strfind(outnameWR,'_view');
                if isempty(in)
                    outdirnameWR = outnameWR(1:strlength(outnameWR)-4);
                    outshortnameWR = 'WR.tif';
                else
                    outdirnameWR = outnameWR(1:in-1);
                    outshortnameWR = outnameWR(in+1:strlength(outnameWR));
                end
                savedirfullWR = [param.savedirWR,filesep,outdirnameWR];
                if ~exist(savedirfullWR,'dir')
                    mkdir(savedirfullWR);
                end
                imgsave32seq(Wiener, fullfile(savedirfullWR,outshortnameWR), app);          
            end   
        
        end
        
    
        %% HiFi-SIM
       
        if makeHF==1
            
            % Step 1
            % if size(Iraw,3)==9
            if param.nrBands==2
                wFilter1=WienerFilterW1_2D(param);
            else
                wFilter1=WienerFilterW1_3D(param);
            end
            
            Wk1=otfHiFi./(wFilter1.wDenom+w1^2);
            
            fftInitialHiFi=fftDirectlyCombined.*Wk1.*Mask;
            % Temp3=real(ifft2(fftshift((fftInitialHiFi))));
            % Temp3(Temp3<0)=0;
            % MIJ.createImage(Temp3);
            
            % Step 2
            if size(Iraw,3)==9
                wFilter2=WienerFilterW2_2D(param);
            else
                wFilter2=WienerFilterW2_3D(param);
            end
            
            % 
            if ApoFWHM==0
                ApoFWHM=0.5*(cutoff-1);
                ApoFWHM=min(0.5,round(ApoFWHM*100)/100);
            end
            apo=apodize_gauss([2*h,2*w], struct('rad',ApoFWHM));  
            Wk2=apo./(wFilter2.wDenom+w2^2);
            fftHiFi=real(ifft2(fftshift((fftInitialHiFi.*Wk2.*Mask))));
        
        end
        
    
        %% Results of wide field
        if makeWF==1
          
            WF=param.WF2;
            if size(WF,1)~=size(WF,2)
                WF=importImages2(WF);
            end    
        %   MIJ.createImage(WF);      % The reconstruction results are displayed in the imageJ window
        %     figure, imshow(WF,[]);    % The reconstruction results are displayed in the matlab window
        %     colormap('hot');
            param.WF2=WF;
            disp('  widefield reconstruction done');
    
            if savePRS==0
                if splitSlices==0
                    outnameWF=strrep(param.filename, '.tif', '_WF.tif');
                else
                    outnameWF=strrep(param.filename, '.tif', ['_s',num2str(sliceid,'%03d'),'_WF.tif']);
                end
                imgsave32seq(WF, fullfile(param.savedirWF,outnameWF), app);
            else
                outnameWF=param.filename;
                in = strfind(outnameWF,'_view');
                if isempty(in)
                    outdirnameWF = outnameWF(1:strlength(outnameWF)-4);
                    outshortnameWF = 'WF.tif';
                else
                    outdirnameWF = outnameWF(1:in-1);
                    outshortnameWF = outnameWF(in+1:strlength(outnameWF));
                end
                savedirfullWF = [param.savedirWF,filesep,outdirnameWF];
                if ~exist(savedirfullWF,'dir')
                    mkdir(savedirfullWF);
                end
                imgsave32seq(WF, fullfile(savedirfullWF,outshortnameWF), app);          
            end   
    
        end
    
    
        %% Results of HiFi-SIM
        if makeHF==1
            
            Size1=2*param.Size1;
            Size2=2*param.Size2;
            HiFi=zeros(Size1,Size2);
            Temp=fftHiFi;
            Temp(Temp<0)=0;
            if normImages==1
                Temp=255*Temp/max(max(Temp));
            end
            HiFi(1:Size1,1:Size2)=Temp(1:Size1,1:Size2);
            HiFi=importImages2(HiFi);
            param.HiFi=HiFi;
            %    MIJ.createImage(HiFi);  % The reconstruction results are displayed in the imageJ window
            %    figure,imshow(HiFi,[]);   % The reconstruction results are displayed in the matlab window
            %   colormap('hot');
            disp('  HiFi reconstruction done');
                
            if savePRS==0
                if splitSlices==0
                    outnameHF=strrep(param.filename, '.tif', '_HF.tif');
                else
                    outnameHF=strrep(param.filename, '.tif', ['_s',num2str(sliceid,'%03d'),'_HF.tif']);
                end
                imgsave32seq(HiFi, fullfile(savedirHF,outnameHF), app);
            else
                outnameHF=param.filename;
                in = strfind(outnameHF,'_view');
                if isempty(in)
                    outdirnameHF = outnameHF(1:strlength(outnameHF)-4);
                    outshortnameHF = 'HF.tif';
                else
                    outdirnameHF = outnameHF(1:in-1);
                    outshortnameHF = outnameHF(in+1:strlength(outnameHF));
                end
                savedirfullHF = [param.savedirHF,filesep,outdirnameHF];
                if ~exist(savedirfullHF,'dir')
                    mkdir(savedirfullHF);
                end
                imgsave32seq(HiFi, fullfile(savedirfullHF,outshortnameHF), app);          
            end        
       
        end % end if HR

    end % end for loop on slices

end % end for loop on files

tEnd = toc(tStart);
endMessage = 'Batch HiFi-SIM finished after ' + string(tEnd) + ' seconds';
disp(endMessage);