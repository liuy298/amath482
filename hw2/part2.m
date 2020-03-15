clear all; close all; clc;

%[y,Fs] = audioread('music1.wav');
[y,Fs] = audioread('music2.wav');
L = length(y)/Fs; % record time in seconds
%p8 = audioplayer(y,Fs);
%playblocking(p8);
v = y.';
n = length(v);
t = (1:length(v))/Fs;
k=(2*pi/L)*[0:n/2-1 -n/2:-1];
ks=fftshift(k);

%plot(ks,fftshift(abs(fft(v))));

tslide = 0:0.2:L;
spec = [];
Hertz = [];
a = 60;
for ii = 1:length(tslide)
    g = exp(-a*(t-tslide(ii)).^2);
    vg = g.*v;
    vgt = fft(vg);
    spec = [spec; fftshift(abs(vgt))/max(abs(vgt))];
    [M,I] = max(abs(vgt));
    Hertz = [Hertz; abs(k(I))];
end

figure(2)
pcolor(tslide,ks/(2*pi),spec.'), shading interp
xlabel('Time[sec]');
ylabel('Frequency[Hz]');
%set(gca,'Ylim',[200 350]);
%yticks(linspace(200,350,16));
set(gca,'Ylim',[750 1100]);
yticks(linspace(750,1100,36));
colormap(hot);
%saveas(gcf,'piano.png');
saveas(gcf,'record.png');