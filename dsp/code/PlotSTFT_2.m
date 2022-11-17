function PlotSTFT_2(T, F, S, win)
    wlen = length(win);
    C = sum(win)/wlen;
    S = abs(S)/wlen/C;
    
    S = 20*log10(S + 1e-6);

    %figure(1)
    surf(T, F, S)
    shading interp;
    axis tight;
    view(0, 90);
    
    xlabel('Time, s');
    ylabel('F');
    
    
    hcol = colorbar;
    
    ylabel(hcol, 'dB');
end
