%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Region-edge-based active contours driven by hybrid and local 
%   fuzzy region-based energy for image segmentation"(HLFRA)
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 23th, Oct, 2018
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = HLFRA_v1(Img,u0,Ksigma,lambda1,lambda2,alpha1,alpha2,belta1,belta2,diswght)

  epsilon = 1;
  timestep = 0.02;
  mu = 1.5; % level set regularization term, please refer to "Chunming Li 
  nu = 5;%length term

  u1 = u0.^2;
  u2 = (1-u0).^2;
  
  Iu1 = Img.*u1;
  Iu2 = Img.*u2;
  
  c1 = sum(sum(Iu1))/sum(sum(u1));
  c2 = sum(sum(Iu2))/sum(sum(u2));  

  Ku1 = imfilter(u1,diswght,'replicate'); 
  Ku2 = imfilter(u2,diswght,'replicate'); 
         
  KI1 = imfilter(Iu1,diswght,'replicate');
  KI2 = imfilter(Iu2,diswght,'replicate');
         
  fb = sum(sum(KI1))/sum(sum(Ku1));
  fs = sum(sum(KI2))/sum(sum(Ku2));  
  
   
 
  DH1 = (Img-fb).^2;
  DH2 = (Img-fs).^2;
  
  F1_old = (Img-c1/2-fb/2).^2.*u1; %The old energy defined in Eq. (26)
  F2_old = (Img-c2/2-fs/2).^2.*u2; %The old energy defined in Eq. (27)
  
  DHW1 = imfilter(DH1,diswght,'replicate');  
  DHW2 = imfilter(DH2,diswght,'replicate');
  
  F3_old = DHW1.*u1; %The old energy defined in Eq. (28)
  F4_old = DHW2.*u2; %The old energy defined in Eq. (29)
  
  % Define in Eq.(21)
  un= 1./(1+(lambda1*((Img-c1/2-fb/2).^2)+alpha1*DHW1)./(lambda2*((Img-c2/2-fs/2).^2)+alpha2*DHW2));
  
  un1 = un.^2;
  un2 = (1-un).^2;
  
  LIun1 = Img.*un1;
  LIun2 = Img.*un2;
  
  KIn1 = imfilter(LIun1,diswght,'replicate'); 
  KIn2 = imfilter(LIun2,diswght,'replicate'); 
  
  Kun1 = imfilter(un1,diswght,'replicate');
  Kun2 = imfilter(un2,diswght,'replicate');
  
  nc1 = sum(sum(LIun1))./sum(sum(un1));
  nc2 = sum(sum(LIun2))./sum(sum(un2));

  ns1 = sum(sum(KIn1))./sum(sum(Kun1));
  ns2 = sum(sum(KIn2))./sum(sum(Kun2));
  
  F1_new = (Img-nc1/2-ns1/2).^2.*un1; %The new energy defined in Eq. (30)
  F2_new = (Img-nc2/2-ns2/2).^2.*un2; %The new energy defined in Eq. (31)
  
  NDH1 = (Img-ns1).^2;
  NDH2 = (Img-ns2).^2;
  
  F3_new = imfilter(NDH1,diswght,'replicate').*un1; %The new energy defined in Eq. (32)
  F4_new = imfilter(NDH2,diswght,'replicate').*un2; %The new energy defined in Eq. (33)
  
  deltaF = lambda1*(F1_new-F1_old)+lambda2*(F2_new-F2_old)+alpha1*(F3_new - F3_old)+alpha2*(F4_new - F4_old); %The change of energy functional
 
  idx = find(deltaF<0);
  u0(idx)=un(idx); 
 
  u = u0;   
  u = imfilter(u,Ksigma,'replicate'); 
  Delta = Dirac(u,epsilon);
  K=curvature_central(u);
  P=mu*(4*del2(u) - K);
  L=nu.*Delta.*K;%length term   
  u = u+timestep*(belta1*L+belta2*P);
end

function K = curvature_central(u)
   [bdx,bdy]=gradient(u);
   mag_bg=sqrt(bdx.^2+bdy.^2)+1e-10;
   nx=bdx./mag_bg;
   ny=bdy./mag_bg;
   [nxx,nxy]=gradient(nx);
   [nyx,nyy]=gradient(ny);
   K=nxx+nyy;
end

function f = Dirac(x, epsilon)
   f=(epsilon/pi)./(epsilon^2.+x.^2);
end




