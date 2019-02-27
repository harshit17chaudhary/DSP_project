# DSP_project
Audio Compression using wavelets in MATLAB
Algorithms compared:
  Haar Algorithms 
  Daubenches Algorithms 
Parameters of Comparison
  PSNR : Peak signal-to-noise ratio   
  NRMSE: Normalised root-mean-square error 
  Compression ratios
MATLAB Inbuilt functions used:
  audioread
  Audiowrite
  length
  ceil
  Sqrt
  DCT
  IDCT
  wavedec
  ddencmp
  wdencmp
  fft
  compand
  waverec

Concepts involved:
  Discrete Cosine Transform
  Inverse Discrete Cosine Transform
  Windowing

Results:
Haar Algorithm

PSNR :              55.4912
MSE:                0.4285
Compression Ratio:  1.9998

Daubenches Algorithm

SNR :               69.7492
MSE:                0.0830
Compression Ratio:  2.1594

Conclusion:

1:After analysing the observations we see that the compression ratio of the haar algorithm is less than daubenches algorithm so the audio compressed by daubenche’s algorithm will require lesser space than that by haar’s algorithm. 

2:Daubenches algorithm is best suited for lossless compression of speech signals as it has more PSNR  and substantially low NMSE. 




