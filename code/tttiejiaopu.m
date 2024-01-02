clear all;
clc;
close all;

%% TIE强度传输方程进行相位恢复
piesize=8e-3;%像素大小
L= 288*piesize;   %相位物体的长 
W= 288*piesize;   %相位物体的宽 512
lambda=632.8e-6 ;%光源波长
k=2*pi/lambda;
d=20;%衍射距离

I1=im2double(imread('dog=19.9mm.tif'));
I2=im2double(imread('dog=20mm.tif'));
I3=im2double(imread('dog=20.1mm.tif'));
diaotaz=0.1;%衍射图之间的距离
D=(I3-I1)./(2*diaotaz);

[x, y, ~] = size(I2);              %size返回行数和列数。r_dim为行数，c_dim为列数，~占位，表示只要行和列的值
fX = [0:fix(x/2),ceil(x/2)-1:-1:1]/L;%fix向下取整，ceil向上取整  %linspace(-x/2,x/2-1,x)/L;
fY = [0:fix(y/2),ceil(y/2)-1:-1:1]/L;%linspace(-y/2,y/2-1,y)/L;%
[fy,fx] = meshgrid(fY,fX); 

q=fx.^2+fy.^2;
pesai=ifft2(q).*fft2(k.*D);

H=exp(1j*k*d.*sqrt(1-(lambda*lambda).*(q)));
HB=1./H;

a=ifft2(gradient((gradient(pesai))./I2));
a=a/max(max(a));
 
m=ifft2(q);
phase=-m.*a;
phase=phase/max(max(phase));

EO=I2.*exp(1j.*phase);
Ei=ifft2((fft2(EO)).*HB);
faik=angle(Ei);
faik=abs(faik);
faik(find(faik==0))=0.01;%去除矩阵中等于0的点

faik=faik/max(max(faik));
% faik=sqrt(faik);
faik=mat2gray(faik);
figure;
imshow(faik);
imwrite(faik,'btie2.tif');



%% 加速角谱迭代法 开始迭代 tie恢复的相位作为角谱的输入相位迭代

step=2000;
loss=ones(step,1);%MSE
psn=zeros(step,1);%psnr
N=288;%
gk=zeros(N,N);
minloss=1;
%% 读入
A0=I2;
%A0=sqrt(A0);
%A0=A0./max(max(A0));
A=ones(N,N);
phasek=2*pi.*rand(N,N);
phasek1=phasek;
Ei=A.*exp(1j.*faik);%初始的物面
figure;
tic
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
   imshow(faik);
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
 imshow(AOO);
 title('恢复');
faik=im2uint8(faik);
figure,imshow(faik);
% imshow(abs(EO));
% imshow(angle(Ei));
imwrite(faik,'tie+jp_bd=2.tif')