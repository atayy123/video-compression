% Motion compensated video coding
% parameters: video frames and stepsize
% returns: distortion (mse) and rate (entropy per frame)
function [d, entro] = coder_mc(Y, stepsize, en)
num_frames = length(Y);
lambda = 0.0005*(stepsize^2);

% write functions for 8x8 2d dct
A = dctmtx(8);
dct = @(block_struct) A * block_struct.data * A';

% for the first frame, apply intra mode (transform and quantize)
sel = Y{1};
transform = blockproc(sel, [8 8], dct);
quantized = stepsize * round(transform/stepsize);

% divide 1. frame to blocks and get parameters for blocks
[rows, cols] = size(sel);
y_l = rows/16;
x_l = cols/16;
numblocks = x_l * y_l;

% rates for different modes (select to)
R_intra = en(log2(stepsize)-2)*256 + 1;
R_mc = 1; % entropy of motion vectors + entropy of residual frame + 1

% needed parameters for Lagrangian cost function
% calculate distortion for the first frame
d = immse(transform, quantized)/num_frames;
% calculate the rate for first frame
ent = numblocks * R_intra;

% storage for previous image (TODO: maybe adapt and take whole picture)
storage = cell(y_l, x_l); % just frame Y{f-1}

% code all frames with either mc or intra-coder
for f = 2:num_frames
    % dct transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    % decide on the mode for each block based on the lagrangian for each
    % mode
    for i = 1:x_l
        for j = 1:y_l
            % preparation for calculating J for mc

            % calculation of mc vectors (maybe write extra function(previous_image, current block, measurement window))
            % calculation of the residual
            % reconstructed image = residual + blocks from previous image
            % shifted by motion vector

            % calculation of the rate (entropy of motion vectors + entropy of residual frame + 1)
            % calculation of the distortion (immse(original,reconstructed))
            
            % get one block of the current frame
            block = transform(16*(j-1)+1:16*j, 16*(i-1)+1:16*i);
            
            % make decision on mode
            % calculate Lagrangian for mc mode
            d_mc = immse(block, storage{j,i});
            L_mc = d_mc + lambda * R_mc;
            % calculate Lagrangian for intra mode
            quantized = stepsize * round(block/stepsize);
            d_intra = immse(quantized, block);
            L_intra = d_intra + lambda * R_intra;
            
            % select frame based on minimal value of Lagrangian
  
            % mc mode
            % update distortion (distortion of first frame through quantization + distortion through copying)
            if L_mc <= L_intra 
                d = d + d_mc/(num_frames*numblocks);
                ent = ent + R_mc;
            % intra mode
            else
                d = d + d_intra/(num_frames*numblocks); % bc loops over all frames and all blocks
                ent = ent + R_intra;
                % change the stored block
                storage{j,i} = quantized;
            end
        end
    end
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
% 
% questions:
% does the intermode replace the copy mode? --> displacement vector (0,0)
