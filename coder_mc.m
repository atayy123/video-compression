% Motion compensated video coding
% parameters: video frames and stepsize
% returns: distortion (mse) and rate (entropy per frame)
function [d, entro] = coder_mc(Y, stepsize, en)
num_frames = length(Y);
lambda = 0.0003*(stepsize^2);

% write functions for 8x8 2d dct
A = dctmtx(8);
dct = @(block_struct) A * block_struct.data * A';
% prepare sorting the coefficients for entropy calculation
ent = @(block_struct) reshape(transpose(block_struct.data), 1, []);

% first, apply inter-mode to all frames to get the rate of inter mode
% for the first frame, apply intra mode (transform and quantize)
sel = Y{1};
% divide 1. frame to blocks and get parameters for blocks
[rows, cols] = size(sel);
frame_size = rows*cols;
y_l = rows/16;
x_l = cols/16;
numblocks = x_l * y_l;
% apply intra mode
transform = blockproc(sel, [8 8], dct);
quantized = stepsize * round(transform/stepsize);

% create memory for motion vectors and residual images
motionVec = zeros(2, numblocks, num_frames-1);
% 99 blocks from 49 frames, 256 coefficients
% entropy for residual blocks
residual_VLC = zeros(99*49, 256);

% storage for previous image
storage = quantized;

for f = 2:num_frames
    count = 0;
    % DCT transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    % residual image
    res_full = zeros(rows,cols);
    % reconstructed image
    recons = zeros(rows,cols);
    % blockwise operation
    for i = 1:x_l
        for j = 1:y_l
            count = count +1;
            % get one block of the current frame
            x_ind = 16*(i-1)+1;
            y_ind = 16*(j-1)+1;
            block = transform(y_ind:y_ind+15, x_ind:x_ind+15);
            % initialize search
            lowest = inf;
            sel_dx = 0; 
            sel_dy = 0;
            for dx = -10:10
                for dy = -10:10
                    % check if the block is out of bounds
                    if (x_ind + dx > 0) && (y_ind + dy > 0) && (rows > y_ind + dy + 15) && (cols > x_ind + dx + 15)
                        % calculate error
                        err = immse(block, storage(y_ind+dy:y_ind+dy+15,x_ind+dx:x_ind+dx+15));
                        % if error is lowest, save the motion vector
                        if err < lowest
                            lowest = err;
                            sel_dx = dx;
                            sel_dy = dy;
                        end
                    end
                end
            end
            % calculate the residual image
            residual = block-storage(y_ind+sel_dy:y_ind+sel_dy+15,x_ind+sel_dx:x_ind+sel_dx+15);
            % quantize the residual
            residual_q = stepsize * round(residual/stepsize);
            % save the reconstructed block
            res_full(y_ind:y_ind+15, x_ind:x_ind+15) = residual_q;
            recons(y_ind:y_ind+15, x_ind:x_ind+15) = residual_q + storage(y_ind+sel_dy:y_ind+sel_dy+15,x_ind+sel_dx:x_ind+sel_dx+15);
            % save the motion vectors
            motionVec(:, count, f-1) = [sel_dx, sel_dy];
        end
    end
    % prepare entropy
    entro_temp = blockproc(res_full, [16 16], ent);
    entro_temp = reshape(entro_temp', 256, frame_size/256)';
    residual_VLC(((f-1)*99)+1:99*f,:) = entro_temp;
    % store the current frame
    storage = recons;
end

% calculate the entropy
mot_entro = entropy(mat2gray(motionVec));
res_entro = 0;
% sum over the entropies of all coefficients of all frames to get the
% total entropy per block
for j = 1:256
    % select jth coefficient from all blocks
    sel = residual_VLC(:,j);
    res_entro = res_entro + entropy(mat2gray(sel))/256;
end

% rates for different modes, since there are 3 modes, 2 bits needed to
% represent the mode
R_intra = en(log2(stepsize)-2)*256 + 2;
R_mc = 256*res_entro + mot_entro + 2; % entropy of motion vectors + entropy of residual frame + 2
R_copy = 2;

% after calculating the rates, start video coding with motion compensation
% for the first frame, apply intra mode (transform and quantize)
sel = Y{1};
%apply intra mode
transform = blockproc(sel, [8 8], dct);
quantized = stepsize * round(transform/stepsize);

% calculate distortion for the first frame
d = immse(transform, quantized)/num_frames;
% calculate the rate for first frame
ent = numblocks * R_intra;

% storage for previous image
storage = quantized;

% code all frames
for f = 2:num_frames
    % dct transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    % reconstructed image
    recons = zeros(rows,cols);
    % decide on the mode for each block based on the lagrangian for each
    % mode
    for i = 1:x_l
        for j = 1:y_l
            % get one block of the current frame
            x_ind = 16*(i-1)+1;
            y_ind = 16*(j-1)+1;
            block = transform(y_ind:y_ind+15, x_ind:x_ind+15);
            % initialize search
            lowest = inf;
            sel_dx = 0; 
            sel_dy = 0;
            for dx = -10:10
                for dy = -10:10
                    % check if the block is out of bounds
                    if (x_ind + dx > 0) && (y_ind + dy > 0) && (rows > y_ind + dy + 15) && (cols > x_ind + dx + 15)
                        % calculate error
                        err = immse(block, storage(y_ind+dy:y_ind+dy+15,x_ind+dx:x_ind+dx+15));
                        % if error is lowest, save the motion vector
                        if err < lowest
                            lowest = err;
                            sel_dx = dx;
                            sel_dy = dy;
                        end
                    end
                end
            end
            % calculate the residual image
            residual = block-storage(y_ind+sel_dy:y_ind+sel_dy+15,x_ind+sel_dx:x_ind+sel_dx+15);
            % quantize the residual
            residual_q = stepsize * round(residual/stepsize);
            % calculate the reconstructed block
            recons_mc = residual_q + storage(y_ind+sel_dy:y_ind+sel_dy+15,x_ind+sel_dx:x_ind+sel_dx+15);
            
            % make decision on mode
            % calculate Lagrangian for copy mode
            d_copy = immse(block, storage(y_ind:y_ind+15,x_ind:x_ind+15));
            L_copy = d_copy + lambda * R_copy;

            % calculate Lagrangian for mc mode
            d_mc = immse(block, recons_mc);
            L_mc = d_mc + lambda * R_mc;

            % calculate Lagrangian for intra mode
            quantized = stepsize * round(block/stepsize);
            d_intra = immse(quantized, block);
            L_intra = d_intra + lambda * R_intra;
            
            % select frame based on minimal value of Lagrangian
            L = [L_copy, L_mc, L_intra];
            [~, ind] = min(L);
  
            if ind == 2
                % mc mode
                d = d + d_mc/(num_frames*numblocks);
                ent = ent + R_mc;
                recons(y_ind:y_ind+15, x_ind:x_ind+15) = recons_mc;
            elseif ind == 3
                % intra mode
                d = d + d_intra/(num_frames*numblocks); % bc loops over all frames and all blocks
                ent = ent + R_intra;
                % change the stored block
                recons(y_ind:y_ind+15, x_ind:x_ind+15) = quantized;
            else
                % copy mode
                d = d + d_copy/(num_frames*numblocks);
                ent = ent + R_copy;
                recons(y_ind:y_ind+15, x_ind:x_ind+15) = storage(y_ind:y_ind+15,x_ind:x_ind+15);
            end
        end
    end
    storage = recons;
end
% summed entropy over all frames has to be averaged for one frame
entro = ent/(num_frames);
end

% notes
% 1. calculate cost function for all three modes:
% J = D + lambda*R
% 
% see part 2 and 3 for calculation
% motion-compensation: 
% - rate (entropy calculation): 
% prediction (motion vectors for each block) + residual (coded by intra-frame coder, reestimate statistics) + number of bit for intermode (1 or 2?)
% - distortion: 
% mse (original, reconstructed image)
% 
% - residual: (original - reconstructed image) 
% - reconstructed image: 	sum of motion compensated block (for each block take the block from the previous image shifted by vector)
% 			and residual block (original - motion compensation)
% - reestimate statistics: 
% 	perform motion compensation for all frames
% 	compute coefficients of all residuals
% 	calculate average entropy for one frame (take entropies of all coefficients individually and add the up) 
% 	--> used for calculation of all frames
% 
% - motion estimator: select best motion vector for each block (minimizes the mse of the prediction)
% 	--> make a block procedure for current image 
% 	--> loop over all possible displacement windows and take the same block shifted by all displacement vectors (from the previous image)
% 	--> calculate MSEs 
% 	--> pick vector for the block with the lowest MSE
% 
% - entropy for motion vectors: 
% 
% 2. mode selection through if-else
