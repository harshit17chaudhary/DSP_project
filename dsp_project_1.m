clc;
close all;
clear all;

% reading the audio file

file='Windows XP Startup.wav'; %this file is used throught the code
[y,samp_f] = audioread(file);
len=length(y);
T=1/samp_f;
x=0:T:(length(y)-1)*T;
plot(x,y);
title('orignal Audio');
xlabel('time In Seconds');
ylabel('Amplitude');

figure;


%% Initiallizing the enable variables

wavelet='haar';
level=5;
framewidth=2048;
psychoacoustic='on '; %if it is off it uses 8 bits/frame as default
w_comp_status = 'on ';
compander='on ';
quantization ='on ';
heavy_compression='off';

%% ENCODER 

step=framewidth;
N=ceil(len/step);

% computational variables

% initializing variables
Cchunks=0;
Lchunks=0;
Csize=0;
PERF0mean=0;
PERFL2mean=0;
n_avg=0;
n_max=0;
n_0=0;
n_vector=[];

for i=1:1:N
    if (i~=N);
        frame=y([(step*(i-1)+1):step*i]);
    else
        frame=y([(step*(i-1)+1):length(x)]);
    end

    % wavelet decomposition of the frame
    [C,L] = wavedec(frame,level,wavelet);

    % wavelet compression scheme
    % denoising the signal
    if w_comp_status =='on '
        [thr,sorh,keepapp] = ddencmp('cmp','wv',frame);
        % ddencmp used to find the threshold whcih is use dto denoise the
        % signal
        if heavy_compression == 'on '
            thr=thr*10^6;
        end
        [XC,CXC,LXC,PERF0,PERFL2] = wdencmp('gbl',C, L, wavelet,level,thr,sorh,keepapp);
        C=CXC;
        L=LXC;
        PERF0mean=+ PERF0;
        PERFL2mean=+PERFL2;
    end

    %Psychoacoustic model
    if psychoacoustic =='on '
            % Taking the DB represntation of the absolute value of the frequency
            % response of he frame
            P=10.*log10((abs(fft(frame,length(frame)))).^2);
        % Empty array of zeros of the length of the DB value    
        Ptm=zeros(1,length(P));

        % Inspect spectrum and find tones maskers
        for k=1:1:length(P)
            % setting a boolean value to decide the value of ptm for each
            % itteration
            if (((k<=1)||(k>=250))||(P(k)<P(k-1))||(P(k)<P(k+1)))
                bool = 0;
            elseif ((k>2) & (k<63)),
                bool = ((P(k)>(P(k+2)+7))& (P(k)>(P(k-2)+7)));
            elseif ((k>=63) & (k<127)),
                bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-3)+7)) & (P(k)>(P(k+3)+7)));
            elseif ((k>=127) & (k<=256)),
                bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-3)+7)) & (P(k)>(P(k+3)+7)) & (P(k)>(P(k-4)+7)) & (P(k)>(P(k+4)+7)) &(P(k)>(P(k-5)+7)) & (P(k)>(P(k+5)+7)) & (P(k)>(P(k-6)+7)) &(P(k)>(P(k+6)+7)));
            else
                bool = 0;
            end

            %  If the boolean value is 1 we update the ith value in the array
            %  ptm it else remains 0 
            if bool==1
                Ptm(k)=10*log10(10.^(0.1.*(P(k-1)))+10.^(0.1.*(P(k)))+10.^(0.1.*P(k+1)));
            end
        end

        % initializing total energy
        totalenergy=0;

        % calculating the total energy
        for k=1:1:length(Ptm)
            totalenergy=+ 0.^(0.1.*(Ptm(k)));
        end

        % calculating SNR and n
        E=10*log10(totalenergy/(length(Ptm))); 
        % SNR = 6.02*n + 1.76 db
        SNR=max(P)-E;
        n=ceil(SNR/6.02);

        if n<=3
            n=4;
            n_0=+1;
        end
        % updating the max value
        if n>=n_max
            n_max=n;
        end

    n_avg= + n;

    n_vector=[n_vector n];
    
    end
    
    %Compander(compressor)
    % compressign based on U-law
    if compander=='on '
        Mu=255;
        C = compand(C,Mu,max(C),'mu/compressor'); 
        % orignalvalue is replaced by the compressed value
    end
    
    %Quantization
    if quantization=='on '
        if psychoacoustic=='off'
            n=8;
        end
        partition = [min(C):((max(C)-min(C))/2^n):max(C)];
        codebook = [1 min(C):((max(C)-min(C))/2^n):max(C)];
        % quantizing the variable
        [index,quant,distor] = quantiz(C,partition,codebook);
        
        %find and correct offset
        
        offset=0;
        for j=1:1:N
            if C(j)==0
                offset=-quant(j);
                break;
            end
        end
        quant=quant+offset;
        C=quant;
    end
    % Putting all together all the chunks
    Cchunks=[Cchunks, C]; 
    Lchunks=[Lchunks, L];
    
    Csize=[Csize length(C)];
    
    Encoder = round((i/N)*100); %indicator of progess
end

Cchunks=Cchunks(2:length(Cchunks));
Csize=[Csize(2) Csize(N+1)];
Lsize=length(L);
Lchunks=[ Lchunks(2:Lsize+1), Lchunks((N-1)*Lsize+1:length(Lchunks))];

%Indicators

PERF0mean=PERF0mean/N;
PERFL2mean=PERFL2mean/N;
n_avg=n_avg/N;
n_max;

%% Decoder

xdchunks=0;

for i=1:1:N;
    
    if i==N;
        Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(2)+(Csize(1)*(i-1))]);
        %Compander (expander)
        if compander=='on '
            if max(Cframe)==0
            else
             Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
            end
        end
        xd = waverec(Cframe,Lchunks(Lsize+2:length(Lchunks)),wavelet);
    else
        Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(1)*i]);
        %Compander (expander)
        if compander=='on '
            if max(Cframe)==0
            else
                Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
            end
        end
    xd = waverec(Cframe,Lchunks(1:Lsize),wavelet);
    end
    
    xdchunks=[xdchunks xd];
    Decoder = round((i/N)*100); %indicator of progess
    
end

xdchunks=xdchunks(2:length(xdchunks));

%%

%creating audio files with compressed schemes

audiowrite('output1.wav',xdchunks,samp_f);
[x,samp_f] = audioread('output1.wav');
len=length(x);
t=0:1/samp_f:(length(x)-1)/samp_f;
plot(t,xdchunks);
title('Commpressed signal');
xlabel('time In Seconds');
ylabel('Amplitude');


%% Calculating the comparison metrics

% reading the reconstructed data and orignal data

[y1,samp_f1]=audioread(file);
[y2,samp_f2]=audioread('output1.wav');

%%  calculating the size of the data samples

[r,c]=size(y1);
[c2x,c2y]=size(y1);

%% Error calculation

error = (sum(y1(2)-y2).^2)/(r*c);
MSE=sqrt(error);

%% PSNR calculation

MAXVAL=255;
% calculating in DB
PSNR = 20*log10(MAXVAL/MSE); 


if(MSE <= 0)
    PSNR = 99;
end

%%  compression ratio calculation

% Size of the compressed file
directory2 = dir('output1.wav');
s = directory2.bytes;
s1 = s/1024;

%  SIze of orignal file
directory = dir(file);
S = directory.bytes;
S1 = S/1024;

% Ratio
Compressionratio = S1/s1;

% displaying all the values
disp('Compression ratio:');
disp(Compressionratio);
disp('Peak Signal to noise ratio:');
disp(PSNR);
disp('Mean Square Root  Error:');
disp(MSE);


%% Important functions used
% audioread
% audiowrite
% length
% ceil
% sqrt
% wavedec
% ddencmp
% wdencmp
% fft
% compand
% waverec


