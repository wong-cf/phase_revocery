close all;clear all;clc; %
iterative=3000;            %设迭代次数为300次吧
imagename='dog=20mm.tif';    %你想要提取相位的图像名称
phaseimage='b.tif';  %要保存的相位图像名称
figure(1),imshow(imagename);
%空域输入图像的幅度（是已知的，也就是清晰的图像，它的灰度就是幅值）和相位图像（待恢复）
known_abs_spatial=imread(imagename);            %作为输入图像的幅度，是已知的
%known_abs_spatial =rgb2gray(known_abs_spatial);%注意要用单通道图像做实验，如果你读取的是彩色图像，那就吧这行取消注释变成灰度图像吧
known_abs_spatial=im2double(known_abs_spatial); %将图像灰度映射到0～1
unknown_phase=known_abs_spatial;                %Peppers图像作为输入图像的相位，也即为待恢复的数据，
                                                %要求它和known_abs_spatial大小一致，所以这里直接赋值就好了
unknown_phase=im2double(unknown_phase);         %将图像灰度映射到0～1
unknown_phase2=unknown_phase*2*pi;              %相位范围映射到0-2*pi
unknown_phase2(unknown_phase2>pi)=unknown_phase2(unknown_phase2>pi)-2*pi;%进一步映射至[-pi,+pi]
[width,length]=size(known_abs_spatial);         %获取图像的大小
input=known_abs_spatial.*exp(1i*unknown_phase2); %最终输入图像:幅度*e^(i*相位角度)，它是复数图像
known_abs_fourier=abs(fft2(input));             %先将input图像进行傅立叶变换，然后取模，就是傅氏变换后的幅度
%以下开始迭代求相位
phase_estimate=pi*rand(width,length);           %这是生成了一副大小为(width*length)的图像
                                                %它的像素值是在[0,pi]范围内随机生成的。
figure(2),imshow(phase_estimate)
%以下开始迭代
for p=1:iterative
    signal_estimate_spatial=known_abs_spatial.*exp(1i*phase_estimate);   %Step 1  构造estimated signal：还是幅度*e^(i*相位角度)变成复数形式
    temp1=fft2(signal_estimate_spatial);                                %傅立叶变换到频域
    temp_ang=angle(temp1);                                              %求相位弧度，它的范围是[-pi,pi]
    signal_estimate_fourier=known_abs_fourier.*exp(i*temp_ang);         %Step 2  替换傅氏变换后的幅度，产生estimate Fourier transform
    temp2=ifft2(signal_estimate_fourier);                               %Step 3  对Step 2产生的estimate Fourier transform进行傅立叶反变换，又变换到空域了
    phase_estimate=angle(temp2);                                        %Step 4:estimated phase
%     IS=abs(abs(temp2)-abs(temp1)).^2;
%     MSE=sum(IS(:))/256^2%计算均方误差
   
   
end
%以上循环就是通过随便预设一个相位图像，在循环中不断调整逼近真实的相位，直到满足条件（也就是我们求的相位和真实相位非常接近的时候）
%不过这里我们只需要设定一个比较大的循环就可以了，基本上都可以满足条件了，这个激光原理就讲过了。
phase_estimate(phase_estimate<0)=phase_estimate(phase_estimate<0)+2*pi; %把estimate_phase从[-pi,+pi]，映射到[0,2pi]
retrieved=phase_estimate/(2*pi);%再映射到[0,1]

%     IS=abs(abs(temp2)-abs(temp1)).^2;
%     MSE1=sum(IS(:))/256^2%计算均方误差
%     figure(2);plot(log10(MSE),'LineWidth',1.5);
%     xlabel('Iterative number');%迭代的次数
%     ylabel('Logarithm of Mean Square Error');%均方误差的对数
   
figure (3)
imshow(retrieved,[]);title('相位图像')%显示我们提取到的相位图像
a=min(min(retrieved));
 retrieved=retrieved/max(max(retrieved));
% imwrite(retrieved,phaseimage)

Uz=known_abs_spatial.*exp(1j.*retrieved);
piesize=8e-3;%像素大小
L= 288*piesize;   %相位物体的长 
W= 288*piesize;   %相位物体的宽 512
lambda=632.8e-6 ;%光源波长
k=2*pi/lambda;
d=20;%衍射距离

[x, y, ~] = size(known_abs_spatial);              %size返回行数和列数。r_dim为行数，c_dim为列数，~占位，表示只要行和列的值
fX = [0:fix(x/2),ceil(x/2)-1:-1:1]/L;%fix向下取整，ceil向上取整  %linspace(-x/2,x/2-1,x)/L;
fY = [0:fix(y/2),ceil(y/2)-1:-1:1]/L;%linspace(-y/2,y/2-1,y)/L;%
[fy,fx] = meshgrid(fY,fX); 
q=fx.^2+fy.^2;
H=exp(1j*k*d.*sqrt(1-(lambda*lambda).*(q)));
HB=1./H;

Eii=ifft2((fft2(Uz)).*HB);
phase1=angle(Eii);
phase1(find(phase1==0))=0.01;
imshow(phase1,[]);

%% 加速角谱迭代法 开始迭代 tie恢复的相位作为角谱的输入相位迭代

step=500;
loss=ones(step,1);%MSE
psn=zeros(step,1);%psnr
N=288;%
gk=zeros(N,N);
minloss=1;
%% 读入
A0=known_abs_spatial;
%A0=sqrt(A0);
%A0=A0./max(max(A0));
A=ones(N,N);
phasek=retrieved;
phasek1=phasek;
Ei=A.*exp(1j.*phase1);%初始的物面
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
 imshow(faik);
 title('恢复');
faik=im2uint8(faik);
figure,imshow(faik);
% imshow(abs(EO));
% imshow(angle(Ei));
imwrite(faik,'GS+jp_dogd=20.tif')
