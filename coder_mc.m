% Motion compensated video coding
% parameters: video frames, stepsize and entropy of intra-frame coding
% returns: distortion (mse) and rate (entropy per frame)

% TODO: add copy mode
function [d, entro] = coder_mc(Y, stepsize, en_intra)
num_frames = length(Y);
lambda = 0.0005*(stepsize^2);
mv = {}; % motion vectors

% write functions for 8x8 2d dct
A = dctmtx(8);
dct = @(block_struct) A * block_struct.data * A';
inverse_dct = @(block_struct) A' * block_struct.data * A;

% for the first frame, apply intra mode (transform and quantize)
sel = Y{1};
transform = blockproc(sel, [8 8], dct);
quantized = stepsize * round(transform/stepsize);

% divide 1. frame to blocks and get parameters for blocks (for copy mode)
[rows, cols] = size(sel);
y_l = rows/16;
x_l = cols/16;
numblocks = x_l * y_l;

% rates for different modes (per block)
R_intra = en_intra*256 + 1;

% needed parameters for Lagrangian cost function
% calculate distortion for the first frame
d = immse(transform, quantized)/num_frames;
% calculate the rate for first frame
ent = numblocks * R_intra;

% perform inter mode on all frames to get statistics
for f = 2:num_frames
    % dct transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    quantized = stepsize * round(transform/stepsize);
    
    % get statistics of mc (residual and motion vectors)
    for i = 1:x_l
        for j = 1:y_l
            % calculation of mc vectors (done for pixels, doesn't make sense for coefficients)
            mv{end+1} = motion_estimator(Y{f},Y{f-1},[i,j]);
            % calculation of the residual: original block - mc block
            x_mc_start = 16*(i-1)+1 + mv{end}(1);
            x_mc_end = 16*i + mv{end}(1);
            y_mc_start = 16*(j-1)+1 + mv{end}(2);
            y_mc_end = 16*j + mv{end}(2);
            residual{i,j} = Y{f}(16*(j-1)+1:16*j, 16*(i-1)+1:16*i) - Y{f-1}(y_mc_start:y_mc_end, x_mc_start:x_mc_end);
            % reconstructed image = residual (transform, quantized, inverse dct) + block from previous image
            % image shifted by motion vector
            % dct and quantization of residual to obtain reconstructed res
            % (res quantized)
            res_transformed = blockproc(residual{i,j},[8,8],dct);
            res_quantized = stepsize * round(res_transformed/stepsize);
            res_reconstructed = blockproc(res_quantized,[8,8],inverse_dct);
            mc_block = Y{f-1}(y_mc_start:y_mc_end, x_mc_start:x_mc_end);
            reconstructed{i,j} = res_reconstructed + mc_block;
        end
    end
end
% R_mc
ent_mv = entropy_vectors(mv); % calculate entropy of motion_vectors
ent_residual = 0; % interframe coding
R_mc = ent_mv + ent_residual;

% code all frames with either mc or intra-coder
for f = 2:num_frames
    % dct transform of the current frame
    sel = Y{f};
    transform = blockproc(sel, [8 8], dct);
    % decide on the mode for each block based on the lagrangian
    % iterate over all blocks
    for i = 1:x_l
        for j = 1:y_l
            % TODO: calculation of the rate (entropy of motion vectors + entropy of residual frame + 1)
            R_mc = 0;
            % get one block of the current frame
            original_block = transform(16*(j-1)+1:16*j, 16*(i-1)+1:16*i);
            % calculation of the distortion (immse(original,reconstructed))
            d_mc = immse(original_block, reconstructed{i,j});
            quantized = stepsize * round(original_block/stepsize);
            d_intra = immse(quantized, original_block);
            % calculate Lagrangians
            J_mc = d_mc + lambda * R_mc;
            J_intra = d_intra + lambda * R_intra;
            
            % select frame based on minimal value of Lagrangian
  
            % mc mode
            % update distortion (distortion of first frame through quantization + distortion through mc)
            if J_mc <= J_intra 
                d = d + d_mc/(num_frames*numblocks);
                ent = ent + R_mc;
            % intra mode
            else
                d = d + d_intra/(num_frames*numblocks); % bc loops over all frames and all blocks
                ent = ent + R_intra;
            end
        end
    end
end
% summed entropy over all frames has to be averaged for one frame
entro = ent/(num_frames);
end
