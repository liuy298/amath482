clear all; close all; clc;
load handel

v = y';
L=9; n=length(v);
v(end) = [];
t = (1:length(v))/Fs;
k = (2*pi/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k);

subplot(2,1,1)
plot(t,v)
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

vt=fft(v);
subplot(2,1,2)
plot(ks,fftshift(abs(vt))/max(abs(vt)))
xlabel('Frequency [\omega]');
ylabel('Amplitude');
title('Signal in Frequency Domain');
saveas(gcf,'signal.png')

% Gabor, Mexican hat wavelet, and Shannon wavelet
a= 1000;
tslide = 0:0.5:9;
vgt_spec=[];
for ii = 1:length(tslide)
    % Gabor filter
    g = exp(-a*(t-tslide(ii)).^2);
    % Mexican hat wavelet
    %g = (2/(sqrt(3*a)*pi^(1/4)))*(1-((t-tslide(ii)).^2/a)).*exp(-(t-tslide(ii)).^2/(2*a^2));
    % Shannon
    %g = abs(t - tslide(ii)) < a;
    vg = g.*v;
    vgt = fft(vg);
    vgt_spec = [vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,t,g)
    xlabel('time (sec)'), ylabel('v(t),g(t)')
    subplot(3,1,2),plot(t,vg)
    xlabel('time (sec)'), ylabel('v(t)g(t)')
    subplot(3,1,3),plot(t,fftshift(abs(vgt))/max(abs(vgt)))
    xlabel('frequency (\omega)'), ylabel('FFT(vg)')
    axis([0 9 0 1])
    drawnow
    pause(0.1)
end


figure(2)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time [sec]');
ylabel('frequency [\omega]');
colormap(hot)
drawnow, hold on
saveas(gcf,'spectrograms.png')