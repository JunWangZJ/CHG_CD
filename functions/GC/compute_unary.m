function D_min = compute_unary(rgb_pt, gmm_params)
%COMPUTE_UNARY Part of GrabCut. Compute unary (data) term
%
% Inputs:
%   - rgb_pt: a RGB point   单个点像素        
%   - gmm_params: GMM parameters
%
% Output:
%   - D_min: minimum D
%
% Author:
%   Xiuming Zhang
%   xiuming6zhang[on]gmail.com
%   Dept. of ECE, National University of Singapore
%   April 2015
%

% The pixel's D to all Gaussians in this GMM
pix_D = zeros(size(gmm_params, 1), 1);

% For every Gaussian
for idx = 1:size(gmm_params, 1)    % 多类别  多均值多sigma
    
    mu = gmm_params{idx, 2};
    sigma = gmm_params{idx, 3};
    
    diff_i=rgb_pt-mu;
    
    val = diff_i*inv(sigma).*diff_i+log(det(sigma));  % 当前类别下的能量值
    
    pix_D(idx) = val;    %都存起来
    
end


D_min = min(pix_D);   % 选最小类别值
