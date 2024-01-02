%% 自己写的根据角谱迭代法进行相位恢复
clc,clear
close all
%% 参数初始化
lambda=632.8e-6;%波长
d=20;%衍射距离mm
N=288;%像素
PIESIZE=8e-3;%像素大小
L=N*PIESIZE;%长宽
k=2*pi/lambda;%波矢
step=100;
loss=ones(step,1);%MSE
psn=zeros(step,1);%psnr
b=1.1;%修正量
a=0.8;
gk=zeros(N,N);
minloss=1;
%% 读入
A0=im2double(imread('dog=20mm.tif'));
% A0=double(imread('4.25d=38_456.tif'));
% A0=sqrt(A0);
% A0=A0./max(max(A0));
A=ones(N,N);
phasek=2*pi.*rand(N,N);
phasek1=phasek;
Ei=A.*exp(1j.*phasek);%初始的物面
%figure;
%imshow(angle(Ei));

figure;
%% 频域初始化
[x,y,~]=size(Ei);
fX=[0:fix(x/2),ceil(x/2)-1:-1:1]./L;
fY=[0:fix(y/2),ceil(y/2)-1:-1:1]./L;
[fx,fy]=meshgrid(fX,fY);

%% 角谱传播函数
f=fx.^2+fy.^2;
H=exp(1j*k*d.*sqrt(1-(lambda*lambda).*(f)));
HB=1./H;

tic
%% 开始迭代
for n=1:step
    EOO=ifft2((fft2(Ei)).*H);
    AOO=abs(EOO).^2;
    AOO=AOO./max(max(AOO));
    EO=A0.*exp(1j.*angle(EOO));%新相位 像面
    Eii=ifft2((fft2(EO)).*HB);
    faik=angle(Eii);  %新相位 物面
    faik=faik./max(max(faik));
    beitak=(phasek-phasek1);
    %abk=(A0-AOO)./pi./2;
    if n>1
       gk1=gk;
       gk=faik-phasek;
       rk=sum((gk.*gk1),"all")/(sum((gk1.^2),"all"));%abs
       phasek1=phasek;
       phasek=faik+beitak*rk;
       phasek=phasek./max(max(phasek));
    else 
        gk=faik-phasek;
        phasek=faik;
    end
    
    Ei=exp(1j*phasek);
    loss(n)=immse(A0,AOO);
    psn(n)= 10 * log10(1/loss(n));
   imshow(faik);%物面的相位
    if loss(n)<minloss
        
         %imwrite(faik,fullfile([num2str(loss(n)) '.tif']))
    end
    minloss=min(loss);
   
end
toc
 figure;
imshow(A0);
title('原图');
 figure;
 imshow(faik);
 title('恢复');
faik=im2uint8(faik);
figure,imshow(faik);
% imshow(abs(EO));
% imshow(angle(Ei));
imwrite(faik,'dog=20mr.tif');
test = A0.*(1j.*angle(EOO));
test = ifft2(test);
imshow(test);
%% 保存数据
% save('MSE.txt','loss','-ascii');
% save('PSNR.txt','psn','-ascii');