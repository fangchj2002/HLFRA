%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Region-edge-based active contours driven by hybrid and local 
%   fuzzy region-based energy for image segmentation"(HLFRA)
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 23th, Oct, 2018
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wght = disweight(rad)
    
    dist = ones(2*rad+1,2*rad+1);
    
    for m=1:2*rad+1
        for n=1:2*rad+1
            dist(m,n) = 1/(1+sqrt((m-rad)^2+(n-rad)^2));
        end
    end
    wght = dist/sum(sum(dist));      
end