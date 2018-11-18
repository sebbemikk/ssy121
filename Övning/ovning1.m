clc
cla 
clear all
close all


Ts=0.025;


fsamp=1/(2*Ts)
tcont = [-6*Ts:0.0001:6*Ts];
tstep = -6*Ts:1/fsamp:6*Ts;
p = @(t) 6*sinc(t/(2*Ts)).^2 + sinc(t/(6*Ts)).^6 + 6*sinc(t/(10*Ts)).^10;


pstep=p(tstep);
pcont=p(tcont);

stem(tstep,p(tstep))
figure()
%plot(tcont,p(tcont))
N = length(pstep);
P=fftshift(fft(pstep))
df = fsamp/N;

xf = (-floor(N/2):1:ceil(N/2)-1)*df
%freq=linspace(-1000,1000,length(ft));
plot(xf,20*log10(P))


p_rec = zeros(1, length(tcont))

for n=-floor(N/2):floor(N/2)
    currentIndex = n + floor(N/2) +1
    
    p_rec = p_rec + p(currentIndex) * sinc(fs(ts-n*Tcont))
end

plot(tcont
