function [B,H,z] = nrldpc_Hmatrix(BG)
    
    % load the base graph given in input (BG)
    load(sprintf('%s.txt',BG),BG);
    if strcmp(BG, 'NR_2_6_52')
        z = 52;
        BG = NR_2_6_52;
    else 
        z = 352;
        BG = NR_1_5_352;
    end
    B = BG;
    [n, m] = size(B);
    % Initialize H matrix
    H = zeros(n*z, m*z);
    Iz = eye(z);
    I0 = zeros(z);

    % Lifting
    for i=1:n
        rw = (i-1)*z+(1:z);
        for j=1:m
            col = (j-1)*z+(1:z);
            if B(i, j)==-1
                H(rw, col) = I0;
            else
                H(rw, col) = circshift(Iz, -B(i, j));
                % -B(i, j) because we want to shift Iz up (or right) and circshift default shifts the matrix downward
            end
        end
    end
end