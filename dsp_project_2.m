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

%% Initializing the data for compression

% reading the data from windows XP Startup.wav for compressing the signals
[Data, samp_f] = audioread('Windows XP Startup.wav');
%chosing a block size 
wsize = 8192; % std value for testing

%initializing compressed matrice

Comp_data2 = [];
Comp_data4 = [];
Comp_data8 = [];


%changing compression percentages
half_Samples = wsize / 2;
quater_Samples = wsize / 4;
eigth_Samples = wsize / 8;

%% compression Algorithm

for i=1:wsize:length(Data)-wsize
    
    % DCT  = Discrete cosine transform
    % IDCT = Inverse Discrete cosine transform
    
    windowDcT = dct(Data(i:i+wsize-1));
    
    Comp_data2(i:i+wsize-1) = idct(windowDcT(1:half_Samples), wsize);
    Comp_data4(i:i+wsize-1) = idct(windowDcT(1:quater_Samples), wsize);
    Comp_data8(i:i+wsize-1) = idct(windowDcT(1:eigth_Samples), wsize);
end

%%  Writing the the compresssed data in another file

% using audiowrite - in 3 files we are writing compressed data for
% differernt compression ratio

%compressed file
[x,samp_f] = audioread('output1.wav');
t=0:1/samp_f:(length(x)-1)/samp_f;
plot(t,x);
title('Commpressed signal ');
xlabel('time In Seconds');
ylabel('Amplitude');
figure;


% Commpression Percentage windowsize/2
audiowrite('output3.wav',Comp_data2,samp_f)
[x,samp_f] = audioread('output3.wav');
t=0:1/samp_f:(length(x)-1)/samp_f;
subplot(3,1,1)
plot(t,x);
title('Commpression Percentage windowsize/2 ');
xlabel('time In Seconds');
ylabel('Amplitude');

% Commpression Percentage windowsize/4
audiowrite('output4.wav',Comp_data4,samp_f)
[x,samp_f] = audioread('output4.wav');
t=0:1/samp_f:(length(x)-1)/samp_f;
subplot(3,1,2);
plot(t,x);
title('Commpression Percentage windowsize/4 ');
xlabel('time In Seconds');
ylabel('Amplitude');


% Commpression Percentage windowsize/8
audiowrite('output5.wav',Comp_data8,samp_f)
[x,samp_f] = audioread('output5.wav');
t=0:1/samp_f:(length(x)-1)/samp_f;
subplot(3,1,3)
plot(t,x);
title('Commpression Percentage windowsize/8');
xlabel('time In Seconds');
ylabel('Amplitude');

%% Calculating the comparison metrics

% reading the reconstructed data and orignal data

[y1,samp_f1]=audioread(file);
[y2,samp_f2]=audioread('output3.wav');

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
directory2 = dir('output3.wav');
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

%% functions used
% Audioread
% Audiowrite
% Length
% DCT
% IDCT
% Sqrt
