% Showing the gray image
image = 'DSC_1568a.tif';
x = imread(image);
x = double(x);
r = x(:,:,1);
g = x(:,:,2);
b = x(:,:,3);
xgray = 0.2116*r + 0.7152*g + 0.0722*b;
xgray = uint8(xgray);

% Define block size
block_size = 8;

% Get image dimensions
[rows, cols] = size(xgray);

% Initialize the result image
result_image = zeros(rows, cols, 'uint8');

% Loop through the image in 8x8 blocks
for row = 1:block_size:rows
    for col = 1:block_size:cols
        % Extract an 8x8 block
        block_row = row:min(row + block_size - 1, rows);
        block_col = col:min(col + block_size - 1, cols);
        X1 = xgray(block_row, block_col);

        % Get DCT matrix
        n = block_size;
        DCTMatrix = DCT(n);

        % Get quantization matrix
        % Edit loss parameter to change degree of compression
        Q = quant(1);

        % Compress and then decompress the block
        Yq = compress(X1, DCTMatrix, Q);
        Xf = decompress(Yq, DCTMatrix, Q);

        % Put the block back into the result image
        result_image(block_row, block_col) = Xf;
    end
end

% Show the result image
imagesc(result_image, [0,255]);
colormap('gray');
title('Processed Image');

% Discrete Cosine Transform that takes in matrix length n
function C = DCT(n)
    % Initialize the DCT coefficients matrix
    C=zeros(n, n);
    for i=1:n
        for j=1:n
            C(i,j)=cos((i-1)*(2*j-1)*pi/(2*n));
        end
    end
    C=sqrt(2/n)*C;
    C(1,:)=C(1,:)/sqrt(2);
end

% Linear quantization matrix
% Define p as the loss parameter, that determines visual accuracy
function Q = quant(p)
    Q=p*8./hilb(8);
end

% X is the original matrix, C is the DCT matrix,
% Q is the linear quantization
function Yq = compress(X, C, Q)
    Xd=double(X);
    Xc=Xd-128;
    Y=C*Xc*C';
    Yq=round(Y./Q);
end

% Yq is the compressed matrix, C is the DCT matrix
% Q is the linear quantization
function Xf = decompress(Yq, C, Q)
    Ydq=Yq.*Q;
    Xdq=C'*Ydq*C;
    Xe=Xdq+128;
    Xf=uint8(Xe);
end