clear all
nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nBitPerSym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym = 1; % number of symbols

freqOffsetkHz_v = [-200:10:200];
EbN0dB = 30; % bit to noise ratio
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

for ii = 1:length(freqOffsetkHz_v)

   % Transmitter
   ipBit = ones(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod = 2*ipBit-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbolsa

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nSym,6) ipMod(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   xt = [xt(:,[49:64]) xt];

   % Concatenating multiple symbols to form a long vector
   xt = reshape(xt.',1,nSym*80);

   % Adding frequency offset
   xt = xt.*exp(j*2*pi*freqOffsetkHz_v(ii)*(1e3/20e6)*[0:length(xt)-1]); 

   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,nSym*80) + j*randn(1,nSym*80)];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB/20)*nt;

   % Receiver
   yt = reshape(yt.',80,nSym).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix

   % converting to frequency domain
   yF = sqrt(64/80)*(sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).'; 
   yMod = yF(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]); 

   err = (yMod - ipMod) ;
   err = err(:).';
   errdB = 10*log10(err*err'/length(err));
   errdB_v(ii) = errdB;

   delta = freqOffsetkHz_v(ii)/312.5;
   theoryErr =sum(1./(j*2*pi*([-5:5]+delta)).*(exp(j*2*pi*([5:5]+delta))-1));
   theoryErr = (theoryErr-1);
   theoryErr_v(ii) = theoryErr;
   theoryErrdB(ii)   = 10*log10(theoryErr*theoryErr'/length(theoryErr));
   theoryErrdB(21) = -EbN0dB;

end


close all; figure
plot(freqOffsetkHz_v./312.5,theoryErrdB,'bs-','LineWidth',2);
hold on;
plot(freqOffsetkHz_v./312.5,errdB_v,'mx-','LineWidth',2);
axis([-0.7 0.7 -30 10])
grid on
legend('theory','simulation');
xlabel('freqency offset/subcarrier spacing')
ylabel('Error, dB')
title('Error magnitude with frequency offset')

