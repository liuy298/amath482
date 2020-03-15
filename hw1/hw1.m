clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes

x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

utsum = zeros(n,n,n);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    utsum = utsum + fftn(Un);
end
utave = abs(fftshift(utsum))/20;
utnorm = utave/max(utave(:));


[ix,iy,iz] = ind2sub(size(utnorm), find(utnorm==max(utnorm(:))));
kx0 = Kx(ix,iy,iz);
ky0 = Ky(ix,iy,iz);
kz0 = Kz(ix,iy,iz);

tau=0.2;
filter  = exp(-tau * ((Kx - kx0).^2 + (Ky - ky0).^2 + (Kz - kz0).^2));
path = zeros(20,3);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    final_fft = fftn(Un);
    final_unft = fftshift(final_fft);
    final_utf = filter.*final_unft;
    final_uf = abs(ifftn(final_utf));
    marble = final_uf/max(final_uf(:));
    [x,y,z] = ind2sub(size(marble),find(marble==max(marble(:))));
    path(j,1) = X(x,y,z);
    path(j,2) = Y(x,y,z);
    path(j,3) = Z(x,y,z);
end

plot3(path(:,1),path(:,2),path(:,3),'bo-');
xlabel('x')
ylabel('y')
zlabel('z')
grid on
hold on
plot3(path(20,1),path(20,2),path(20,3),'r*');
saveas(gcf,'marble.png')