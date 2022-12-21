% intra frame video coder function for whole frame sequence
% parameters: video frames and stepsize
% returns: distortion and rate
function [d, entro] = intraframe_coder(Y, stepsize)
    % write functions for 8x8 2d dct
    A = dctmtx(8);
    dct = @(block_struct) A * block_struct.data * A';
    % prepare sorting the coefficients for entropy calculation
    ent = @(block_struct) reshape(transpose(block_struct.data), 1, []);
    
    num_frames = length(Y);
    
    % transform all frames of the video
    transformedY = cell(num_frames,1);
    % calulate frame size
    [n,m] = size(Y{1});
    frame_size = n*m;
    
    % 99 blocks from 50 frames (4950), 256 coefficients
    full_VLC = zeros(4950, 256);
    
    % calculate average distortion
    d = 0;
    
    % transform frames and get the VLC for all frames
    for i = 1:num_frames
        % transform frame
        transformed = blockproc(Y{i}, [8 8], dct);
        % quantize coefficients
        result = stepsize * round(transformed/stepsize);
        % calculate distortion using Parseval's theorem and add it up
        d = d + immse(result, transformed)/num_frames; % ???? devide later?
        transformedY{i} = result;
        % prepare entropy
        entro_temp = blockproc(result, [16 16], ent);
        entro_temp = reshape(entro_temp', 256, frame_size/256)';
        % 99 Blocks per image for all frames entro_temp values for all
        % coefficients within one block
        full_VLC(((i-1)*99)+1:99*i,:) = entro_temp;
    end
    
    % calculate the entropy
    entro = 0;
    % sum over the entropies of all coefficients of all frames to get the
    % total entropy per block
    for j = 1:256
        % select jth coefficient from all blocks
        sel = full_VLC(:,j);
        entro = entro + entropy(mat2gray(sel))/256;
    end
end
