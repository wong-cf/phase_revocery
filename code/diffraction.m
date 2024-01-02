%% 角谱传播理论进行衍射图的仿真
clc,clear
close all
%% 初始化参数
d= 20.1%传播距离 mm  300000
piesize=8e-3;%像素大小
L= 288*piesize;   %相 位物体的长 
W= 288*piesize;   %相位物体的宽 512
lambda=632.8e-6 ;%光源波长
%% 构建相位物体  数据预处理
object_amplitude = ones(288,288);  %振幅              %im2double(imread());
object_phase     =  im2double(imread('doggray.tif'));  figure,imshow(object_phase);title('原图');%直接读入im2double(imread('fruit.tif'));
%i=mat2gray(object_phase);figure,imshow(object_phase);
phase=imresize(object_phase,[288,288])./max(max(object_phase));figure;imshow(phase);title('归一化');
%phase=pi*imresize(object_phase,[288,288])./max(max(object_phase));figure;imshow(phase);title('归一化投影到pi');
%phase=2*pi*imresize(object_phase,[288,288])./max(max(object_phase));figure;imshow(phase2);title('归一化投影到2pi');
%phase=*pi*object_phase;figure;imshow(phase3);title('直接乘以2pi');
U0 = exp(1j.*phase);figure,imshow(U0);title('显示实部');%%
a=angle(U0);figure,imshow(a);title('显示相位');
%(1j*pi*2*object_phase);%exp(1j*pi*(object_phase/max(object_phase))); %归一化的复振幅2*pi*
%% 初始化频率   空间频率的计算公式
[x, y, ~] = size(U0);              %size返回行数和列数。r_dim为行数，c_dim为列数，~占位，表示只要行和列的值
fX = [0:fix(x/2),ceil(x/2)-1:-1:1]/L;%fix向下取整，ceil向上取整  %linspace(-x/2,x/2-1,x)/L;
fY = [0:fix(y/2),ceil(y/2)-1:-1:1]/L;%linspace(-y/2,y/2-1,y)/L;%
[fy,fx] = meshgrid(fY,fX); 

%% 对图像做傅里叶变换
U0_fft = fft2(U0);%(U0);

%% 角谱传播
f = sqrt(fx.^2+fy.^2);
%circ_f = circ(f.*lambda);    % 高于1/λ的频率无法传播
Uz_fft = U0_fft.*exp(1j*2*pi*d.*sqrt(1./lambda.^2-f.^2));%.*circ_f;

%% 复原图像
Uz = ifft2(Uz_fft);

%%  显示图像
I=abs(Uz);%幅度
%I=I.*I;%强度
I=I./max(max(I));%归一化
%%imshow(I,[])重要 自动调整显示的灰度，后面用imcrop保存
%I=mat2gray(I);%%%%数据类型转换，非常重要，把double矩阵转化为灰度图，解决了上面imshow无法保存的问题以及【】的替代
 maxium1=max(max(I)); %取二维图像的最大值
 minium1=min(min(I)); %取二维图像的最小值
%  I=I-minium1;
%  I=I./max(max(I));
I=im2uint8(I);
figure;
imshow(I);%colorbar;
%% 保存图片
% I=imcrop(I,[0 0 280 280]);
%     I=imcrop(I);       %直接在图片上进行裁剪！
% imshow(I);
imwrite(I,'dog=20.1mm.tif');



