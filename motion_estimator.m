% estimate best vector based on SSD (sum of squared distances)
% parameters: current block, previous frame
% returns: optimal motion compensation vector for each block (2D)
function motion_vector = motion_estimator(current_frame, previous_frame, current_block_position)
m = current_block_position(1); % x-coordinate
n = current_block_position(2); % y-coordinate
best_ssd = -1;
[rows, cols] = size(current_frame);

% iterate over ssquare measurement window from -10 to 10
for i = -10:10 % x-coordinate
    for j = -10:10 % y-coordinate
        x_start = 16*(m-1)+1+i;
        x_end = 16*m+i;
        y_start = 16*(n-1)+1+j;
        y_end = 16*n+j;
        % only take into account if motion vector points inside the image
        % (using short circuit boolean expressions)
        if x_start >= 1 && x_end <= cols && y_start >= 1 && y_end <= rows
            block_p = previous_frame(y_start:y_end, x_start:x_end);
            block_c = current_frame(16*(n-1)+1:16*n, 16*(m-1)+1:16*m);
            % compare frames
            difference = double(block_c) - double(block_p);
            ssd = sum(difference(:).^2);
            % find mv with minimal ssd
            if ssd <= best_ssd || best_ssd == -1
                best_ssd = ssd;
                best_mv = [i j];
            end
        end
    end
end
motion_vector = best_mv;
end