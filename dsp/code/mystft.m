function [STFT, f, t] = mystft(x, win, hop, nfft, fs)
    x = x(:);
    xlen = length(x);
    wlen = length(win);

    % 窗口数目 L
    L = 1+fix((xlen-wlen)/hop);
    STFT = zeros(nfft, L);
    
    % STFT
    for l = 0:L-1
        % 加窗
        xw = x(1+l*hop : wlen+l*hop).*win;
        
        % FFT计算
        X = fft(xw, nfft);
        X = fftshift(X);

        STFT(:, 1+l) = X(1:nfft);
    end
    
    % 取每个窗口中点的时间点
    t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
    %f = (0:nfft-1)*fs/nfft;
    % 频率 (fftshift之后的)
    f = (-nfft/2:nfft/2-1) * (fs/nfft);
    
end
