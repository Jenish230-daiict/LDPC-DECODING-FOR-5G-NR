function cword = nr5g_encoder(B, z, msg)
    [n, m] = size(B);
    cword = zeros(1, m*z);
    
    % Place the message bits in the beginning of the codeword
    cword(1:(m-n)*z) = msg;
    
    tm = zeros(1, z);
    for i=1:4
        for j=1:m-n
            % Multiply message blocks with base matrix entries
            tm = mod(tm+mul_sh(msg((j-1)*z+1:j*z), B(i, j)), 2);
        end

    end
    if B(2, m-n+1)==-1
        p1_sh = B(3, m-n+1);
    else
        p1_sh = B(2, m-n+1);
    end

    % Compute parity p1
    cword((m-n)*z+1:(m-n+1)*z) = mul_sh(tm, z-p1_sh); % got p1
    % z-p1_sh, because Ik inverse will be I(z-k) (Ik is k right shifted identity)
    
    % Compute parity p2, p3, p4
    for i=1:3
        tm = zeros(1, z);
        for j=1:m-n+i
            % Multiply message blocks with base matrix entries
            tm = mod(tm+mul_sh(cword((j-1)*z+1:j*z), B(i, j)), 2);
        end
        cword((m-n+i)*z+1:(m-n+i+1)*z) = tm;
    end
    
    % Compute remaining parity from p5 to pn
    for i=5:n
        tm = zeros(1, z);
        for j=1:m-n+4
            tm = mod(tm+mul_sh(cword((j-1)*z+1:j*z), B(i, j)), 2);
        end
        cword((m-n+i-1)*z+1:(m-n+i)*z) = tm;
    end

    % totalparity = n*z;
    % info = m-n-2;
    % nbRM = ceil(info/coderate)+2;
    % nBlocklen = nbRM*z;
    % 
    % needed_p = totalparity - (m*z - nBlocklen);
    % total_bits = m*z-n*z+needed_p;
    % 
    % c_word = zeros(1, total_bits);
    % c_word(1:total_bits) = cword(1:total_bits);

    %this puncturing will be done seperatly in simulation after getting full codeword
end