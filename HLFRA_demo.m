%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Region-edge-based active contours driven by hybrid and local 
%   fuzzy region-based energy for image segmentation"(HLFRA)
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 23th, Oct, 2018
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
addpath 'images'
ImgID = 6;
Img = imread([num2str(ImgID),'.bmp']);

[M,N,L] = size(Img);
u = zeros(M,N);

%setting the initial level set function 'u':
u(:,:) = 0.3;
u(40:70,40:60) = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Additional noise:
%  Gaussian noise: Img_gray = imnoise(Img,'gaussian',0,0.02);
%  speckle noise:  Img_gray = imnoise(Img,'speckle',0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ImgID
    case 1
        Img = imnoise(Img,'speckle',0.02);
        Img_gray = rgb2gray(Img);        
        iterNum = 40;
        rad = 3;
        lambda1 = 1;
        lambda2 = 1;
        alpha1 = .001;
        alpha2 = .001;
        belta1 = 1;
        belta2 = 1;
    case 2
        Img_gray = rgb2gray(Img);
        iterNum = 100;
        rad = 3;
        lambda1 = 1.8;
        lambda2 = 1;
        alpha1 = 1;
        alpha2 = 1;
        belta1 = 1;
        belta2 = 1;
    case 3
        Img_gray = Img;        
        iterNum = 40;
        rad = 3;
        lambda1 = .1;
        lambda2 = 1;
        alpha1 = .5;
        alpha2 = .5;
        belta1 = 1;
        belta2 = 1;
    case 4
        Img_gray = Img;
        Img = imnoise(Img,'speckle',0.02);
        iterNum = 40;
        rad = 3;
        lambda1 = 1.8;
        lambda2 = 1;
        alpha1 = 0.3;
        alpha2 = 0.3;
        belta1 = 1;
        belta2 = 1;
    case 5
        Img_gray = Img;
        iterNum = 40;
        rad = 3;
        lambda1 = 1.8;
        lambda2 = 1;
        alpha1 = .01;
        alpha2 = .01;
        belta1 = 1;
        belta2 = 1;
    case 6
        Img_gray = Img;
        iterNum = 40;
        rad = 3;
        lambda1 = 1;
        lambda2 = 1;
        alpha1 = 1;
        alpha2 = 1;
        belta1 = 1;
        belta2 = 1;
    otherwise       
        Img_gray = Img;        
        iterNum = 100;
        rad = 3;
        lambda1 = .1;
        lambda2 = 1;
        alpha1 = .5;
        alpha2 = .5;
        belta1 = 1;
        belta2 = 1;
end

sigma = 3;
Ksigma = fspecial('gaussian',sigma,1.5); % Caussian kernel   

diswght = disweight(rad);

[Ix,Iy] = gradient(double(Img_gray));
f = Ix.^2+Iy.^2;
g = 1./(1+f);  % edge indicator function.


saliency = imfilter(Img_gray,diswght,'replicate');

figure;subplot(2,2,1);imshow(Img);hold on;%axis off,axis equal

title('Initial contour');
[c,h] = contour(u-0.5,[0 0],'r','LineWidth',1);

subplot(2,2,2);
imshow(saliency,[]);hold on;
subplot(2,2,3);

pause(0.1);
tic;
% start level set evolution
for n=1:iterNum
    u = HLFRA_v1(double(Img_gray),u,Ksigma,lambda1,lambda2,alpha1,alpha2,belta1,belta2,diswght);
    if mod(n,5)==0
        pause(0.1);
        imshow(Img, []);hold on;axis off,axis equal
        [c,h] = contour(u-0.5,[0 0],'r','LineWidth',0.5);
        iterNum=[num2str(n), 'iterations'];
        title(iterNum);
        hold off;
    end
end

seg = ((u-0.5)>0);
subplot(2,2,4),imshow(seg);
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);
time = toc;
figure;
mesh(u-0.1);
title('Final level set function');


