% calculate entropy of cell array with matrices
function ent = entropy_vectors(vectors)
    % vectors = {[1 1],[1 1],[0 0],[0 0]};
    [dimension, num_v] = size(vectors);
    int = zeros(1,num_v);
    for i=1:num_v
        % map each vector to an integer number
        int(i) = (vectors{i}(1)+10) * 100 + (vectors{i}(2)+10);
    end
    ent = entropy(mat2gray(int));
end