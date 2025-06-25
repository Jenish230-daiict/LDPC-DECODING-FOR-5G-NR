baseGraph = 'NR_2_6_52';
% Coderates
coderate = [1/4 1/3 1/2 3/5];
Eb_no_db = 0:0.5:10;
colors = lines(length(Eb_no_db));

% Lifting
[B, Hfull, z] = nrldpc_Hmatrix(baseGraph);

% for storing the outputs
decoding_error = zeros(length(coderate), length(Eb_no_db));
bit_error = zeros(length(coderate), length(Eb_no_db));

% Number of simulations
nsim = 1000;

% Maximum iterations
max_it = 20;
iterations = 1:1:max_it;

for rr=1:length(coderate)
    [n, m] = size(B);
    cr = coderate(rr);

    % Adjusting H matrix for specific coderate
    totalparity = n*z;
    info = m-n-2;
    needed_blocks = ceil(info/cr)+2;
    nBlocklen = needed_blocks*z;
    needed_p = totalparity - (m*z - nBlocklen);
    total_bits = n*z-m*z+nBlocklen;
    H = Hfull(:, 1:nBlocklen);
    H = H(1:total_bits, :);


    [row, col] = size(H);
    infob = col-row;

    % Mapping for which check nodes connectd to a VNi
    vn_to_cn = cell(col, 1);
    % Mapping for which variable nodes connectd to a CNi
    cn_to_vn = cell(row, 1);
    % VN->CN and CN->VN msg storing matrix L
    L = zeros(row, col);

    % Mapping
    for i=1:col
        for j=1:row
            if H(j, i)==1
                vn_to_cn{i, 1} = [vn_to_cn{i, 1} j];
            end
        end
    end
    for i=1:row
        for j=1:col
            if H(i, j)==1
                cn_to_vn{i, 1} = [cn_to_vn{i, 1} j];
            end
        end
    end


    % To store output for iteration success (prob. of getting success on iteration i)
    itr_success = zeros(length(Eb_no_db), max_it);
    
    for eb=1:length(Eb_no_db)
        SNR = Eb_no_db(eb);
        SNRL = 10^(SNR/10);
        sigma = sqrt(1/(2*SNRL*cr));
        success = 0;
        error = 0;

        vn_sum = zeros(1, col);
        for sim=1:nsim
            % Generating random msg
            org = randi([0 1], 1, (m-n)*z);
            % Encoding of msg
            encoded_msg = nr5g_encoder(B, z, org);

            % Puncturing
            encoded_msg = encoded_msg(1:nBlocklen);

            % BPSK modulation
            modulated = 1-2.*encoded_msg;
            % Adding noise
            noise = sigma*randn(1, nBlocklen);
            recevied_sig = modulated+noise;

            % Demodulation
            recevied = (recevied_sig<0);
            % prev to check that decode codeword is similar as previous iteration?
            prev = recevied;
            % estimated codeword
            c_aprox = zeros(1, col);

            for it=1:max_it
                if it==1
                    % For first iteration every VN will send its value to connected CNs
                    for i=1:col
                        for j=vn_to_cn{i, 1}
                            % msg to ith VN to jth CN
                            L(j, i) = recevied_sig(1, i);
                        end 
                    end
                else
                    for i=1:col
                        for j=vn_to_cn{i, 1}
                            % subtracting the value sent by the jth cn to avoid positive feedback loop
                            total = vn_sum(1, i)-L(j, i);
                            % msg to ith VN to jth CN
                            L(j, i) = total;
                        end 
                    end
                end

                % Min-sum algo
                for i=1:row
                    mini1 = 1e10;
                    mini2 = 1e10;
                    ind = -1;
                    total_sign = 1;

                    for j=cn_to_vn{i, 1}
                        val = abs(L(i, j));
                        if val<=mini1
                            mini2 = mini1;
                            mini1 = val;
                            ind = j;
                        elseif val<=mini2
                            mini2 = val;
                        end

                        if(L(i, j)~=0)
                            if(L(i, j)<0)
                                total_sign = total_sign*-1;
                            end
                        end
                    end

                    for j=cn_to_vn{i, 1}
                        if j~=ind
                            L(i, j) = total_sign*sign(L(i, j))*mini1;
                        else 
                            L(i, j) = total_sign*sign(L(i, j))*mini2;
                        end
                    end

                end

                % Add values sent by all CNs
                for i=1:col
                    sum_1 = recevied_sig(1, i);
                    tm = L(:, i);
                    sum_1 = sum_1+sum(tm);
                    vn_sum(1, i) = sum_1;
                end

                % Estimate codeword
                c_aprox = (vn_sum<0);

                %check if its the original msg
                check = 1;
                for i=1:infob
                    if(c_aprox(i) ~= org(i)) 
                        check = 0;
                        break
                    end
                end
                if check==1
                    success = success+1;
                    % if yes then we got the success
                    % And if we get success in this iteration we will also get success in remaining iterations
                    for j=it:max_it
                        itr_success(eb, j) = itr_success(eb, j)+1;
                    end
                    break;
                end

                % check if decoded codeword is same as previous iteration
                check2 = 1;
                for i=1:col
                    if c_aprox(1, i) ~= prev(1, i)
                        check2 = 0;
                        break;
                    end
                end
                % if yes then break
                if check2==1
                    break;
                end

                % set prev to decoded codeword
                prev = c_aprox;
            end
            % Count for error bits
            for i=1:col
                if c_aprox(1, i)~=encoded_msg(1, i)
                    error = error+1;
                end
            end
        end
        % calculation for decoding error for coderate coderate(rr) and SNR Eb_no_db(eb)
        decoding_error(rr, eb) = (nsim-success)/nsim;
        % calculation for bit error for coderate coderate(rr) and SNR Eb_no_db(eb)
        bit_error(rr, eb) = error/(nsim*col);

        % we want to plot in logarthmic scale so if bit_error is 0 then log(0) will -inf(which will be ignored in graph) so make it so small value
        if bit_error(rr, eb)==0
            bit_error(rr, eb) = 1e-305;
        end 
    end
    % disp(decoding_error);

    

    % Iteration success probability (Performace graph)
    figure;
    for i=1:length(Eb_no_db)
        plot(iterations,itr_success(i, :)./nsim,'Color',colors(i,:));
        xlabel("Iteration Number");
        ylabel("Success Probability");
        title(['Iteration Success Probability for soft decision decoding, Coderate = ', num2str(cr)]);
        legend('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','8.5','9.0','9.5','10.0'); 
        
        grid on;
        hold on;
    end

end

% Performance Graphs
% Decoding error probability
for i=1:length(coderate)
    figure;
    plot(Eb_no_db, decoding_error(i, :), 'LineWidth', 2);
    xlabel("Eb/No (dB)");
    ylabel("Decoding error probability");
    title(['Soft Decision Decoding Error Probability, Coderate = ', num2str(coderate(i))]);
    grid on;
end


% Bit error probability with normalization comparison
for j=1:length(coderate)
    figure;
    %shannon 

    r = coderate(j);            
    N = 512;             
    EbNo = 10.^(Eb_no_db/10);
    PN_e = zeros(size(EbNo));
    log2e = log2(exp(1));

    for i = 1:length(EbNo)
        P = r * EbNo(i);
        C = log2(1 + P);
        V = (log2e)^2 * (P * (P + 2)) / (2 * (P + 1)^2);
        NA_term = sqrt(N / V) * (C - r + log2(N)/(2*N));
        PN_e(i) = qfunc(NA_term);
    end
    shannonLimit_dB = 10 * log10((2^r - 1)/r);

    semilogy(Eb_no_db, PN_e, 'r-', 'LineWidth', 2); 
    hold on;
    xline(shannonLimit_dB, '--b');
    hold on;
    semilogy(Eb_no_db, bit_error(j, :), 'b-', 'LineWidth', 2);
    legend('Normal Approximation', 'Shannon Limit', 'Simulation');
    grid on;
    hold on;    
    xlabel("Eb/No (dB)");
    ylabel("Bit error probability");
    title(['Soft Decision Bit Error Probability, Coderate = ', num2str(coderate(j))]);
    grid on;

    ylim([1e-30, 100000]);  % 100000, so that legend don't cover the actual graph (So that graph is visible)
    xlim([shannonLimit_dB-0.2, 10]); 
end


% Decoding error probability comparison for all coderates
figure;
for i=1:length(coderate)
    plot(Eb_no_db, decoding_error(i, :), 'LineWidth', 2);
    xlabel("Eb/No (dB)");
    ylabel("Decoding error probability");
    title('Soft Decision Decoding Error Probability Comparison');
    legend('Coderate = 1/4', 'Coderate = 1/3', 'Coderate = 1/2', 'Coderate = 3/5');
    grid on;
    hold on;
end

% Bit error probability comparison for all coderates
figure;
for i=1:length(coderate)
    plot(Eb_no_db, bit_error(i, :), 'LineWidth', 2);
    xlabel("Eb/No (dB)");
    ylabel("Bit error probability");
    title('Soft Decision Bit Error Probability Comparison');
    legend('Coderate = 1/4', 'Coderate = 1/3', 'Coderate = 1/2', 'Coderate = 3/5');
    grid on;
    hold on;
end
