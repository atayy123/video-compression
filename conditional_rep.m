% conditional replenishment video coder function 
% parameters: video frames, stepsize, entropy of intra frame coding 
% returns: distortion (mse) and rate (entropy per frame)
function [d, entro] = conditional_rep(Y, stepsize, en)
num_frames = length(Y);

% write functions for 8x8 2d dct
A = dctmtx(8);
dct = @(block_struct) A * block_struct.data * A';

% for the first frame, apply intra mode (ransform and quantize)
sel = Y{1};
transform = blockproc(sel, [8 8], dct);
quantized = stepsize * round(transform/stepsize);

% calculate distortion for the first frame
d = immse(transform, quantized)/num_frames;
% divide frame to blocks and get parameters for blocks
[rows, cols] = size(sel);
y_l = rows/16;
x_l = cols/16;
numblocks = x_l * y_l;

% rates for different modes (select to)
R_intra = en(log2(stepsize)-2)*256 + 1;
R_copy = 1;
% calculate the rate for first frame
ent = numblocks * R_intra;

storage = cell(y_l, x_l);
% block by block operation to store the blocks of the first frame
for i = 1:x_l
    for j = 1:y_l
        % flatten the block for easier calculation
        block = quantized(16*(j-1)+1:16*j, 16*(i-1)+1:16*i);
        % store the blocks
        storage{j,i} = block;
    end
end

lambda = 0.0005*(stepsize^2);
% conditional replenishment for other frames
for f = 2:num_frames
    % dct transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    % decide on the mode for each block based on the lagrangian for each
    % mode
    for i = 1:x_l
        for j = 1:y_l
            % make decision on mode

            % get one block of the current frame
            block = transform(16*(j-1)+1:16*j, 16*(i-1)+1:16*i);
            % calculate Lagrangian for copy mode
            d_copy = immse(block, storage{j,i});
            L_copy = d_copy + lambda * R_copy;
            % calculate Lagrangian for intra mode
            quantized = stepsize * round(block/stepsize);
            d_intra = immse(quantized, block);
            L_intra = d_intra + lambda * R_intra;
            
            % select frame based on minimal value of Lagrangian
  
            % copy mode
            % update distortion (distortion of first frame through quantization + distortion through copying)
            if L_copy <= L_intra 
                d = d + d_copy/(num_frames*numblocks);
                ent = ent + R_copy;
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


