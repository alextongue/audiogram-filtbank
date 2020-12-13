function fhandleout = ploterror(ref_freqs, ref_mag, query_ir, fs, fhandlein)
    if isempty(fhandlein)
        fhandleout = figure;
    else
        fhandleout = figure(fhandlein);
    end
    
    logplot = true; % (! HARD CODED FOR NOW !)
    fft_len = 2^ceil(log2(numel(query_ir)));
    
    % interpolate
    query_tf = fft(query_ir,fft_len);
    query_mag = 20*log10(abs(query_tf));
    fftfreqs = (0:(fft_len-1)).*fs/fft_len;
    
    ref_mag_interp = interp1(ref_freqs,ref_mag,fftfreqs,[],'extrap');
    
    error_tf = query_mag - ref_mag_interp;
    
    subplot(2,1,1);
    plot(fftfreqs,ref_mag_interp,'k',...
        'linewidth',1.5,'displayname','Audiogram Ref.');
    hold on;
    plot(fftfreqs,query_mag,'k:',...
        'linewidth',1.5,'displayname','after optimization');
    hold off;
    if logplot
        set(gca,'XScale','log');
    end
    xlim([20,fs/2]);
    grid on;
    legend('location','best');
    
    subplot(2,1,2);
    plot(fftfreqs,zeros(size(fftfreqs)),'k',...
        'linewidth',1.5,'displayname','Audiogram Ref');
    hold on;
    plot(fftfreqs,error_tf,'k:',...
        'linewidth',1.5,'displayname','after optimization')
    hold off;
    if logplot
        set(gca,'XScale','log');
    end
    xlim([20,fs/2])
    ylim([-10,10]);
    grid on;
    legend('location','best');
    title('dB Error');
    
    
end


