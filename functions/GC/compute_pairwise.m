function pairwise = compute_pairwise(im_sub, gamma)
%COMOUTE_PAIRWISE Part of GrabCut. Compute the pairwise terms.
%
% Inputs:
%   - im_sub: 2D subimage, on which Graph Cut is performed
%   - gamma: gamma parameter
%
% Output:
%   - pairwise: a dense no_edgesx6 matrix of doubles. Each row is of the
%format [i, j, e00, e01, e10, e11] where i and j are neighbours and the four
%coefficients define the interaction potential
%
% Author:
%   Xiuming Zhang
%   xiuming6zhang[on]gmail.com
%   Dept. of ECE, National University of Singapore
%   April 2015
%

% Get image dimensions
[im_h, im_w, ~] = size(im_sub);

%------- Compute \beta

beta = compute_beta(im_sub);

%------- Set pairwise

pairwise = zeros((im_h-1)*(im_w-1)*2+(im_h-1)+(im_w-1), 6);   %初始化矩阵

% Loop through all the pixels (nodes) and set pairwise
idx = 1;
for y = 1:im_h
    for x = 1:im_w
        % Current node
        node = (x-1)*im_h+y;
        color = double(im_sub(y, x));
        
        % Right neighbor
        if x < im_w % Has a right neighbor
            node_r = (x+1-1)*im_h+y;
            color_r = double(im_sub(y, x+1));
            pairwise(idx, 1) = node;
            pairwise(idx, 2) = node_r;
            pairwise(idx, 3) = compute_V(color, 0, color_r, 0, gamma, beta);
            pairwise(idx, 4) = compute_V(color, 0, color_r, 1, gamma, beta);
            pairwise(idx, 5) = compute_V(color, 1, color_r, 0, gamma, beta);
            pairwise(idx, 6) = compute_V(color, 1, color_r, 1, gamma, beta);
            idx = idx+1;
        end
        
        % Down neighbor
        if y < im_h % Has a down neighbor
            node_d = (x-1)*im_h+y+1;
            color_d = double(im_sub(y+1, x));
            pairwise(idx, 1) = node;
            pairwise(idx, 2) = node_d;
            pairwise(idx, 3) = compute_V(color, 0, color_d, 0, gamma, beta);
            pairwise(idx, 4) = compute_V(color, 0, color_d, 1, gamma, beta);
            pairwise(idx, 5) = compute_V(color, 1, color_d, 0, gamma, beta);
            pairwise(idx, 6) = compute_V(color, 1, color_d, 1, gamma, beta);
            idx = idx+1;
        end
    end
end

end


function beta = compute_beta(im_sub)

% Get image dimensions
[im_h, im_w, ~] = size(im_sub);

beta_sum = 0;
cnt = 0;

for y = 1:im_h
    for x = 1:im_w
        % Current node
        color = double(im_sub(y, x));
        
        % Right neighbor
        if x < im_w % Has a right neighbor
            color_r = double(im_sub(y, x+1));
            beta_sum = beta_sum+norm(color-color_r)^2; 
            % norm -- If X is a vector, this is equal to the Euclidean distance
            cnt = cnt+1;
        end
        % Down neighbor
        if y < im_h % Has a down neighbor
            color_d = double(im_sub(y+1, x));
            beta_sum = beta_sum+norm(color-color_d)^2;
            cnt = cnt+1;
        end
    end
end

beta = 1/(2*(beta_sum/cnt));

end


function V = compute_V(color1, label1, color2, label2, gamma, beta)

V = gamma*double(label1~=label2)*exp(-beta*(norm(color1-color2)^2));   
% color1 当前节点像素值， color2邻域节点像素值   2-norm二范数  因为是图像网格所以不考虑距离问题

end
