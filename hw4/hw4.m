%% Test1
clear all; close all; clc;
music1 = [];
for i = 1:3
    [y, Fs] = audioread(strcat('a00',int2str(i),'.mp3'));
    y = y(:,1)+y(:,2);
    y = resample(y,20000,Fs);
    Fs = 20000;
    part = [];
    for j = 1:5:45
        five = y(Fs*j:Fs*(j+5),1);
        song_spec = abs(spectrogram(five));
        song_spec = reshape(song_spec,1,16385*8);
        part = [part song_spec'];
    end
    music1 = [music1 part];
end

music2 = [];
for i = 1:3
    [y, Fs] = audioread(strcat('e00',int2str(i),'.mp3'));
    y = y(:,1)+y(:,2);
    y = resample(y,20000,Fs);
    Fs = 20000;
    part = [];
    for j = 1:5:120
        five = y(Fs*j:Fs*(j+5),1);
        song_spec = abs(spectrogram(five));
        song_spec = reshape(song_spec,1,16385*8);
        part = [part song_spec'];
    end
    music2 = [music1 part];
end

music3 = [];
for i = 1:3
    [y, Fs] = audioread(strcat('l00',int2str(i),'.mp3'));
    y = y(:,1)+y(:,2);
    y = resample(y,20000,Fs);
    Fs = 20000;
    part = [];
    for j = 1:5:15
        five = y(Fs*j:Fs*(j+5),1);
        song_spec = abs(spectrogram(five));
        song_spec = reshape(song_spec,1,16385*8);
        part = [part song_spec'];
    end
    music3 = [music1 part];
end

musci1_wav = mu_wavelet(music1);
musci2_wav = mu_wavelet(music2);
musci3_wav = mu_wavelet(music3);

feature = 40;
[U,S,V,result,w,sort1,sort2,sort3] = mu_trainer(musci1_wav,musci2_wav,musci3_wav,feature);
result;
[threshold1, threshold2] = mu_threshold(sort1,sort3,sort2); % adjust by hand


%%
figure()
subplot(1,3,1)
histogram(sort1,27); hold on, plot([threshold1 threshold1],[0 10],'r')
set(gca,'Xlim',[0 3],'Ylim',[0 10],'Fontsize',14)
title('epic')
xlabel('threshold 1')
subplot(1,3,2)
histogram(sort2,51); hold on, plot([threshold2 threshold2],[0 10],'r')
set(gca,'Xlim',[0 3],'Ylim',[0 10],'Fontsize',14)
title('relax')
xlabel('threshold 2')
subplot(1,3,3)
histogram(sort3,30); hold on, plot([threshold2 threshold2],[0 10],'r')
set(gca,'Xlim',[0 3],'Ylim',[0 10],'Fontsize',14)
title('classical')
xlabel('threshold 2')
saveas(gcf,'t1_threshold.png');
%% 

test101 = [];
[y, Fs] = audioread('a004.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:45
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
test101 = [test101 part];

% TestNum=size(TestSet,2);

test101_wav = mu_wavelet(test101);
test101Mat = U'*test101_wav;  % PCA projection
pva101 = w'*test101Mat;  % LDA projection

TestNum=size(test101,2);

ResVec = (pva101>threshold1);

disp('Number of mistakes')
labels = [1 1 1 1 1 1 1 1 1];
errNum = sum(abs(ResVec-labels));

disp('Rate of success');
sucRate = 1-errNum/TestNum

%% Test2
clear all; close all; clc;
music1 = [];
[y, Fs] = audioread('m001.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music1 = [music1 part];

music2 = [];
[y, Fs] = audioread('m002.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music2 = [music2 part];

music3 = [];
[y, Fs] = audioread('m003.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music3 = [music3 part];

musci1_wav = mu_wavelet(music1);
musci2_wav = mu_wavelet(music2);
musci3_wav = mu_wavelet(music3);

feature = 20;
[U,S,V,result,w,sort1,sort2,sort3] = mu_trainer(musci1_wav,musci2_wav,musci3_wav,feature);
result;
[threshold1, threshold2] = mu_threshold(sort3,sort2,sort1); % adjust by hand
%%
figure()
subplot(1,3,1)
histogram(sort3,40); hold on, plot([threshold1 threshold1],[0 10],'r')
set(gca,'Xlim',[-2 0],'Ylim',[0 10],'Fontsize',14)
title('Gopala')
xlabel('threshold 1')
subplot(1,3,2)
histogram(sort2,40); hold on, plot([threshold1 threshold1],[0 10],'r')
set(gca,'Xlim',[-2 1],'Ylim',[0 10],'Fontsize',14)
title('Kunda')
xlabel('threshold 1')
subplot(1,3,3)
histogram(sort1,40); hold on, plot([threshold2 threshold2],[0 10],'r')
set(gca,'Xlim',[3 7],'Ylim',[0 10],'Fontsize',14)
title('Nandana')
xlabel('threshold 2')
saveas(gcf,'t2_threshold.png');
%%
test2 = [];
[y, Fs] = audioread('m003.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 200:5:240
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
test2 = [test2 part];

% TestNum=size(TestSet,2);

test2_wav = mu_wavelet(test2);
test2Mat = U'*test2_wav;  % PCA projection
pva2 = w'*test2Mat;  % LDA projection


TestNum=size(test2,2);

ResVec = (pva2>threshold1);

disp('Number of mistakes')
labels = [1 1 1 1 1 1 1 1 1];
errNum = sum(abs(ResVec-labels));

disp('Rate of success');
sucRate = 1-errNum/TestNum

%% Test3
clear all; close all; clc;
music1 = [];
[y, Fs] = audioread('m001.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music1 = [music1 part];

[y, Fs] = audioread('m002.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music1 = [music1 part];


music2 = [];
[y, Fs] = audioread('c001.wav');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 400:5:500
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music2 = [music2 part];

[y, Fs] = audioread('c002.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music2 = [music2 part];

music3 = [];
[y, Fs] = audioread('z001.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:200
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music3 = [music3 part];

[y, Fs] = audioread('z002.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 1:5:50
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
music3 = [music3 part];

musci1_wav = mu_wavelet(music1);
musci2_wav = mu_wavelet(music2);
musci3_wav = mu_wavelet(music3);

feature = 20;
[U,S,V,result,w,sort1,sort2,sort3] = mu_trainer(musci1_wav,musci2_wav,musci3_wav,feature);
result;
[threshold1, threshold2] = mu_threshold(sort2,sort3,sort1); % adjust by hand
%%
figure()
subplot(1,3,1)
histogram(sort2,80); hold on, plot([threshold1 threshold1],[0 10],'r')
set(gca,'Xlim',[-9 -2],'Ylim',[0 10],'Fontsize',14)
title('meditation')
xlabel('threshold 1')
subplot(1,3,2)
histogram(sort1,61); hold on, plot([threshold2 threshold2],[0 10],'r')
set(gca,'Xlim',[-2 1],'Ylim',[0 10],'Fontsize',14)
title('pop')
xlabel('threshold 2')
subplot(1,3,3)
histogram(sort3,50); hold on, plot([threshold2 threshold2],[0 10],'r')
set(gca,'Xlim',[-2 0],'Ylim',[0 10],'Fontsize',14)
title('electric')
xlabel('threshold 2')
saveas(gcf,'t3_threshold.png');

%%
test3 = [];
[y, Fs] = audioread('m003.mp3');
y = y(:,1)+y(:,2);
y = resample(y,20000,Fs);
Fs = 20000;
part = [];
for j = 200:5:240
    five = y(Fs*j:Fs*(j+5),1);
    song_spec = abs(spectrogram(five));
    song_spec = reshape(song_spec,1,16385*8);
    part = [part song_spec'];
end
test3 = [test3 part];

% TestNum=size(TestSet,2);

test2_wav = mu_wavelet(test3);
test2Mat = U'*test2_wav;  % PCA projection
pva2 = w'*test2Mat;  % LDA projection


TestNum=size(test3,2);

ResVec = (pva2>threshold1);

disp('Number of mistakes')
labels = [1 1 1 1 1 0 1 1 0];
errNum = sum(abs(ResVec-labels));

disp('Rate of success');
sucRate = 1-errNum/TestNum

%% Build a function that takes in the music piece 
%  and returns a matrix with the cod_edge for each
function muData = mu_wavelet(mufile)

    [m,n] = size(mufile);
    nw = m/2; % wavelet resolution
    muData = zeros(nw,n);
    
    for k = 1:n
       x = mufile(:,k); 
       [cA,cD] = dwt(x,'haar');
       cod_cA1 = rescale(abs(cA));
       cod_cD1 = rescale(abs(cD));
       cod_edge = cod_cA1+cod_cD1;
       muData(:,k) = reshape(cod_edge,nw,1);
    end

end


function [U,S,V,result,w,sort1,sort2,sort3] = mu_trainer(mu1,mu2,mu3,feature)
    n1 = size(mu1,2); n2 = size(mu2,2); n3 = size(mu3,2);
    
    [U,S,V] = svd([mu1 mu2 mu3],'econ');
    
    genre = S*V'; % projection onto principal components
    U = U(:,1:feature);
    genre1 = genre(1:feature,1:n1);
    genre2 = genre(1:feature,n1+1:n1+n2);
    genre3 = genre(1:feature,n1+n2+1:n1+n2+n3);
    
    m1 = mean(genre1,2);
    m2 = mean(genre2,2);
    m3 = mean(genre3,2);
    
    Sw = 0; % within class variances
    for k=1:n1
        Sw = Sw + (genre1(:,k)-m1)*(genre1(:,k)-m1)';
    end
    for k=1:n2
        Sw = Sw + (genre2(:,k)-m2)*(genre2(:,k)-m2)';
    end
    for k=1:n3
        Sw = Sw + (genre3(:,k)-m3)*(genre3(:,k)-m3)';
    end
    mA = mean([m1 m2 m3]);
    Sb = length(genre1)*(m1-mA)*(m1-mA)' + length(genre2)*(m2-mA)*(m2-mA)'+ length(genre3)*(m3-mA)*(m3-mA)'; % between class 
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    v1 = w'*genre1; 
    v2 = w'*genre2;
    v3 = w'*genre3;
    
    
    array = sort([mean(v1) mean(v2) mean(v3)]);
    result=[0 0 0];
    if array(1) == mean(v1)
        result(1) = 1;
        if array(2) == mean(v2)
            result(2) = 2;
            result(3) = 3;
        elseif array(2) == mean(v3)
            result(2) = 3;
            result(3) = 2;
        end
    elseif array(1) == mean(v2)
        array(1) = 2;
        if array(2) == mean(v1)
            result(2) = 1;
            result(3) = 3;
        elseif array(2) == mean(v3)
            result(2) = 3;
            result(3) = 1;
        end
    else
        array(1) = 3;
        if array(2) == mean(v2)
            result(2) = 2;
            result(3) = 1;
        elseif array(2) == mean(v1)
            result(2) = 1;
            result(3) = 2;
        end
    end
    
    sort1 = sort(v1);
    sort2 = sort(v2);
    sort3 = sort(v3);
end
function [threshold1, threshould2] = mu_threshold(temp1,temp2,temp3)
    t1 = length(temp1);
    t2 = 1;
    while temp1(t1)>temp2(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    threshold1 = (temp1(t1)+temp1(t2))/2;
    
    t3 = length(temp2);
    t4 = 1;
    while temp2(t3)>temp3(t4)
        t3 = t3-1;
        t4 = t4+1;
    end
    threshould2 = (temp2(t3)+temp3(t3))/2;
end