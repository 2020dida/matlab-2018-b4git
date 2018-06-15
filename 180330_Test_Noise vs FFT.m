N=[1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 8192*2 8192*4 8192*8 8192*16 8192*32];
for p=1:length(N)
Scale=1;
XNoise=randn(N(p),1)'*Scale;
XNoise_ZP=XNoise;
XNoise_ZP((length(XNoise_ZP)+1):4*length(XNoise_ZP))=0;
subplot(1,1,1)
plot(XNoise)

FFT_noise=fft(XNoise);
FFT_noise_ZP=fft(XNoise_ZP);

subplot(3,1,1)
plot(XNoise)

subplot(3,1,2)
plot(abs(FFT_noise))
subplot(3,1,3)
plot(abs(FFT_noise_ZP))

STDFFT_real(p)=std(real(FFT_noise));
STDFFT_abs(p)=std(abs(FFT_noise));

%std(abs(FFT_noise_ZP))
end
%
subplot(1,1,1)
plot(log10(N),log10(STDFFT_real),log10(N),log10(Scale*((N/2).^0.5)))
%©Ò¥HSpectral Noise = Noise * (N^0.5)/2
%%
subplot(2,1,1)
plot(N,STDFFT_real.^2,N,(Scale*((N/2).^0.5)).^2)

subplot(2,1,2)
plot(N,STDFFT_abs.^2,N,(Scale*((N/4.68).^0.5)).^2)