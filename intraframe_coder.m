% intra frame video coder function that inputs the video frames and
% stepsize, and outputs distortion and rate
function [d, entro] = intraframe_coder(Y, stepsize)
    % write functions for 8x8 2d dct
    A = dctmtx(8);
    dct = @(block_struct) A * block_struct.data * A';
    % prepare sorting the coefficients for entropy calculation
    ent = @(block_struct) reshape(transpose(block_struct.data), 1, []);
    
    numframe = length(Y);
    
    % transform all frames of the video
    transformedY = cell(numframe,1);
    
    [n,m] = size(Y{1});
    sz = n*m;
    
    full_VLC = zeros(4950, 256);
    
    % calculate average distortion
    d = 0;
    
    % transform frames and get the VLC
    for i = 1:numframe
        transformed = blockproc(Y{i}, [8 8], dct);
        result = stepsize * round(transformed/stepsize);
        % calculate distortion using Parseval's theorem
        d = d + immse(result, transformed)/numframe;
        transformedY{1} = result;
        entro_temp = blockproc(result, [16 16], ent);
        entro_temp = reshape(entro_temp', 256, sz/256)';
        full_VLC(((i-1)*99)+1:99*i,:) = entro_temp;
    end
    
    % calculate the entropy
    entro = 0;
    % sum over the entropies of all coefficients of all frames to get the total entropy
    for j = 1:256
        sel = full_VLC(:,j);
        entro = entro + entropy(mat2gray(sel))/256;
    end
end
