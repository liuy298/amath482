clear all; close all; clc; 
load cam1_1;
load cam2_1;
load cam3_1;
load cam1_2;
load cam2_2;
load cam3_2;
load cam1_3;
load cam2_3;
load cam3_3;
load cam1_4;
load cam2_4;
load cam3_4;

numFrames11 = size(vidFrames1_1, 4);
numFrames21 = size(vidFrames2_1, 4);
numFrames31 = size(vidFrames3_1, 4);
numFrames12 = size(vidFrames1_2, 4);
numFrames22 = size(vidFrames2_2, 4);
numFrames32 = size(vidFrames3_2, 4);
numFrames13 = size(vidFrames1_3, 4);
numFrames23 = size(vidFrames2_3, 4);
numFrames33 = size(vidFrames3_3, 4);
numFrames14 = size(vidFrames1_4, 4);
numFrames24 = size(vidFrames2_4, 4);
numFrames34 = size(vidFrames3_4, 4);
%% Test1
y11 = [];
x11 = [];

filter = zeros(480,640);
filter(:,300:360) = 1;

for j=1:numFrames11
    gray = rgb2gray(vidFrames1_1(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x11(j) = [mean(X)];
    y11(j) = [mean(Y)];
end
[a,b] = max(y11(1:30));
y11 = y11(30:end);
x11 = x11(30:end);
%%
y21 = [];
x21 = [];

filter = zeros(480,640);
filter(:,250:360) = 1;
for j=1:numFrames21
    gray = rgb2gray(vidFrames2_1(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x21(j) = [mean(X)];
    y21(j) = [mean(Y)];
end
[a,b] = max(y21(1:40));
y21 = y21(b:end);
x21 = x21(b:end);

%%
x31 = [];
y31 = [];
filter = zeros(480,640);
filter(200:350,230:490) = 1;
for j=1:numFrames31
    gray = rgb2gray(vidFrames3_1(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x31(j) = [mean(X)];
    y31(j) = [mean(Y)];
end
[a,b] = max(y31(1:40));
y31 = y31(b:end);
x31 = x31(b:end);

%% Plot

x11 = x11(1:min([length(x11),length(x21),length(x31)]));
x21 = x21(1:length(x11));
x31 = x31(1:length(x11));
y11 = y11(1:length(x11));
y21 = y21(1:length(x11));
y31 = y31(1:length(x11));

X = [x11;y11; x21; y21; x31; y31];
[~, n] = size(X);
mn = mean(X, 2);
X = X - repmat(mn,1,n);
[u, s, v] = svd(X, 'econ');
figure()
plot(diag(s)./sum(diag(s)), 'ro')
xlabel('Principal component'); ylabel('Energy weight');
v = v*s;
saveas(gcf,'e1.png');
figure()
subplot(2,1,1)
plot(v(:,1),'Linewidth', 2); hold on;
plot(v(:,2),'Linewidth', 2);
plot(v(:,3),'Linewidth', 2);
plot(v(:,4),'Linewidth', 2);
plot(v(:,5),'Linewidth', 2); 
plot(v(:,6),'Linewidth', 2);
xlabel('Time(frame numbers)')
ylabel('Displacement (pixels)')
title('Case 1')
legend('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
subplot(2,1,2)
plot(1:length(X(1,:)), X(1,:),1:length(X(1,:)), X(2,:),'Linewidth', 2);
ylabel("Displacement (pixels)"); xlabel("Time (frames)"); 
title("Case 1: Original displacement across Z axis and XY-plane (cam 1)");
legend("Z", "XY")
saveas(gcf,'c1.png');
%% Test2
y12 = [];
x12 = [];

filter = zeros(480,640);
filter(170:430,300:450) = 1;

for j=1:numFrames12
    gray = rgb2gray(vidFrames1_2(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x12(j) = [mean(X)];
    y12(j) = [mean(Y)];
end
[a,b] = max(y12(1:30));
y12 = y12(b:end);
x12 = x12(b:end);
%%
y22 = [];
x22 = [];

filter = zeros(480,640);
filter(50:475,165:425) = 1;
for j=1:numFrames22
    gray = rgb2gray(vidFrames2_2(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x22(j) = [mean(X)];
    y22(j) = [mean(Y)];
end
[a,b] = max(y22(1:30));
y22 = y22(b:end);
x22 = x22(b:end);

%%
x32 = [];
y32 = [];
filter = zeros(480,640);
filter(200:350,230:490) = 1;
for j=1:numFrames32
    gray = rgb2gray(vidFrames3_2(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x32(j) = [mean(X)];
    y32(j) = [mean(Y)];
end
[a,b] = max(y32(1:40));
y32 = y32(b:end);
x32 = x32(b:end);
%% Plot
x12 = x12(1:min([length(x12),length(x22),length(x32)]));
x22 = x22(1:length(x12));
x32 = x32(1:length(x12));
y12 = y12(1:length(x12));
y22 = y22(1:length(x12));
y32 = y32(1:length(x12));

X = [x12; y12; x22; y22; x32; y32];
[~, n] = size(X);
mn = mean(X, 2);
X = X-repmat(mn,1,n);
[u, s, v] = svd(X, 'econ');
figure()
plot(diag(s)./sum(diag(s)), 'ro')
xlabel('Principal component'); ylabel('Energy weight');
v = v*s;
saveas(gcf,'e2.png');
figure()
subplot(2,1,1)
plot(v(:,1),'Linewidth', 2); hold on;
plot(v(:,2),'Linewidth', 2);
plot(v(:,3),'Linewidth', 2);
plot(v(:,4),'Linewidth', 2);
plot(v(:,5),'Linewidth', 2); 
plot(v(:,6),'Linewidth', 2);
xlabel('Time(frame numbers)')
ylabel('Displacement (pixels)')
title('Case 3')
legend('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
subplot(2,1,2)
plot(1:length(X(1,:)), X(1,:),1:length(X(1,:)), X(2,:),'Linewidth', 2);
ylabel("Displacement (pixels)"); xlabel("Time (frames)"); 
title("Case 2: Original displacement across Z axis and XY-plane (cam 1)");
legend("Z", "XY")
saveas(gcf,'c2.png');
%% Test3
y13 = [];
x13 = [];

filter = zeros(480,640);
filter(150:450,200:450) = 1;

for j=1:numFrames13
    gray = rgb2gray(vidFrames1_3(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x13(j) = [mean(X)];
    y13(j) = [mean(Y)];
end
[~,b] = max(y13(1:30));
y13 = y13(b:end);
x13 = x13(b:end);
%%
y23 = [];
x23 = [];

filter = zeros(480,640);
filter(100:430,165:425) = 1;
for j=1:numFrames23
    gray = rgb2gray(vidFrames2_3(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x23(j) = [mean(X)];
    y23(j) = [mean(Y)];
 end
[a,b] = max(y23(1:50));
y23 = y23(b:end);
x23 = x23(b:end);
%%
x33 = [];
y33 = [];
filter = zeros(480,640);
filter(160:365,235:595) = 1;
for j=1:numFrames33
    gray = rgb2gray(vidFrames3_3(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x33(j) = [mean(X)];
    y33(j) = [mean(Y)];
end
[a,b] = max(y33(1:30));
y33 = y33(b:end);
x33 = x33(b:end);
%% Plot
x13 = x13(1:min([length(x13),length(x23),length(x33)]));
x23 = x23(1:length(x13));
x33 = x33(1:length(x13));
y13 = y13(1:length(x13));
y23 = y23(1:length(x13));
y33 = y33(1:length(x13));

X = [x13; y13; x23; y23; x33; y33];
[m, n] = size(X);
mn = mean(X, 2);
X = X-repmat(mn,1,n);
[u, s, v] = svd(X, 'econ');
figure()
plot(diag(s)./sum(diag(s)), 'ro')
xlabel('Principal component'); ylabel('Energy weight');
v = v*s;
saveas(gcf,'e3.png');
figure()
subplot(2,1,1)
plot(v(:,1),'Linewidth', 2); hold on;
plot(v(:,2),'Linewidth', 2);
plot(v(:,3),'Linewidth', 2);
plot(v(:,4),'Linewidth', 2);
plot(v(:,5),'Linewidth', 2); 
plot(v(:,6),'Linewidth', 2);
xlabel('Time(frame numbers)')
ylabel('Displacement (pixels)')
title('Case 3')
legend('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
subplot(2,1,2)
plot(1:length(X(1,:)), X(1,:),1:length(X(1,:)), X(2,:),'Linewidth', 2);
ylabel("Displacement (pixels)"); xlabel("Time (frames)"); 
title("Case 3: Original displacement across Z axis and XY-plane (cam 1)");
legend("Z", "XY")
saveas(gcf,'c3.png');
%% Test4
y14 = [];
x14 = [];

filter = zeros(480,640);
filter(220:450,270:470) = 1;

for j=1:numFrames14
    gray = rgb2gray(vidFrames1_4(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x14(j) = [mean(X)];
    y14(j) = [mean(Y)];
end
[~,b] = max(y14(1:30));
y14 = y14(b:end);
x14 = x14(b:end);
%%
y24 = [];
x24 = [];

filter = zeros(480,640);
filter(50:400,170:425) = 1;
for j=1:numFrames24
    gray = rgb2gray(vidFrames2_4(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x24(j) = [mean(X)];
    y24(j) = [mean(Y)];
 end
[a,b] = max(y24(1:30));
y24 = y24(b:end);
x24 = x24(b:end);
%%
x34 = [];
y34 = [];
filter = zeros(480,640);
filter(100:300,270:510) = 1;
for j=1:numFrames34
    gray = rgb2gray(vidFrames3_4(:,:,:,j));
    gray_f = double(gray).*filter;
    [a,~] = max(gray_f(:));
    [X,Y] = find(gray_f > a*11/12);
    x34(j) = [mean(X)];
    y34(j) = [mean(Y)];
end
[a,b] = max(y34(1:30));
y34 = y34(b:end);
x34 = x34(b:end);
%% Plot
x14 = x14(1:min([length(x14),length(x24),length(x34)]));
x24 = x24(1:length(x14));
x34 = x34(1:length(x14));
y14 = y14(1:length(x14));
y24 = y24(1:length(x14));
y34 = y34(1:length(x14));

X = [x14; y14; x24; y24; x34; y34];
[m, n] = size(X);
mn = mean(X, 2);
X = X-repmat(mn,1,n);
[u, s, v] = svd(X, 'econ');
figure()
plot(diag(s)./sum(diag(s)), 'ro')
xlabel('Principal component'); ylabel('Energy weight');
v = v*s;
saveas(gcf,'e4.png');
figure()
subplot(2,1,1)
plot(v(:,1),'Linewidth', 2); hold on;
plot(v(:,2),'Linewidth', 2);
plot(v(:,3),'Linewidth', 2);
plot(v(:,4),'Linewidth', 2);
plot(v(:,5),'Linewidth', 2); 
plot(v(:,6),'Linewidth', 2);
xlabel('Time(frame numbers)')
ylabel('Displacement (pixels)')
title('Case 4')
legend('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')
subplot(2,1,2)
plot(1:length(X(1,:)), X(1,:),1:length(X(1,:)), X(2,:),'Linewidth', 2);
ylabel("Displacement (pixels)"); xlabel("Time (frames)"); 
title("Case 4: Original displacement across Z axis and XY-plane (cam 1)");
legend("Z", "XY")
saveas(gcf,'c4.png');