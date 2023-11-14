% Load original image
% Comment out the image to not be used
% image = 'DSC_1568a.tif';
image = 'DSC_1994a.tif';
disp(['Display image: ' num2str(image)]);
original = imread(image);

% Get image dimensions
block_size = 8;
[rows, cols, ~] = size(original);

% Edit and display the loss parameter p to adjust compression degree
p =32;
disp(['Loss Degree (p): ' num2str(p)]);

% Initialize result image
result = zeros(rows, cols, 3, 'uint8');

% Loop through the image in blocks of the specified size
for row = 1:block_size:rows
    for col = 1:block_size:cols
        % Extract a single block to process
        block_row = row:min(row + block_size - 1, rows);
        block_col = col:min(col + block_size - 1, cols);
        X1 = original(block_row, block_col, :);

        % Get DCT matrix
        n = block_size;
        DCTMatrix = DCT(n);

        % Get quantization matrix by using loss parameter (p)
        Q = quant(p);

        % Compress and then decompress the block
        Yq = compress(X1, DCTMatrix, Q);
        Xf = decompress(Yq, DCTMatrix, Q);

        % Put the block back into the result image
        result(block_row, block_col, :) = Xf;
    end
end

% Show the result color image
imshow(result);
title('Processed Color Image');

% Calculate Discrete Cosine Transform coefficient matrix
% Takes matrix length n as parameter
% Used from text p.528
function C = DCT(n)
    C = zeros(n, n);
    for i = 1:n
        for j = 1:n
            C(i, j) = cos((i - 1) * (2 * j - 1) * pi / (2 * n));
        end
    end
    C = sqrt(2 / n) * C;
    C(1, :) = C(1, :) / sqrt(2);
end

% Calculate Linear Quantization matrix
% Takes in loss parameter p to determine compression degree
% Used from text p.532
function Q = quant(p)
    Q = p * 8 ./ hilb(8);
end

% Compress the image data for RGB blocks of data
% Takes data block (X), DCT matrix (C), and quantization matrix (Q)
% Used from text p.533, edited to accommodate RGB values
function Yq = compress(X, C, Q)
    Yq = zeros(size(X));
    for color = 1:3  % Iterate over the colors R, G, B
        Xd = double(X(:,:,color));
        Xc = Xd - 128;
        Y = C * Xc * C';
        Yq(:,:,color) = round(Y ./ Q);
    end
end

% Decompress the image data for RGB blocks to get processed image data
% Takes compressed data block (Yq), DCT matrix (C), and quant. matrix (Q)
% Used from text p.533, edited to accommodate RGB values
function Xf = decompress(Yq, C, Q)
    Xf = zeros(size(Yq));
    for color = 1:3  % Iterate over the colors R, G, B
        Ydq = Yq(:,:,color) .* Q;
        Xdq = C' * Ydq * C;
        Xe = Xdq + 128;
        Xf(:,:,color) = uint8(Xe);
    end
end