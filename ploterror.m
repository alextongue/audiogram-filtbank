function fhandleout = ploterror(ref_freqs, ref_gains, query_ir, fs, fhandlein, dispname, lineformat, plotref)
    if isempty(fhandlein)
        fhandleout = figure;
    else
        fhandleout = figure(fhandlein);
    end
    
    if isempty(lineformat)
        lineformat = 'k:';
    end
    
    
    logplot = true; % (! HARD CODED FOR NOW !)
    fft_len = 2^ceil(log2(numel(query_ir)));
    
    % interpolate
    query_tf = fft(query_ir,fft_len);
    query_mag = 20*log10(abs(query_tf));
    fftfreqs = (0:(fft_len-1)).*fs/fft_len;
    ref_mag_interp = interp1(ref_freqs,ref_gains,fftfreqs,[],'extrap');
    
    error_tf = query_mag - ref_mag_interp;
    
    subplot(2,1,1);
    hold on;
    if plotref
        plot(fftfreqs,ref_mag_interp,'k','linewidth',1.5,'displayname','Audiogram (Reference)');
    end
    plot(fftfreqs,query_mag,lineformat,'linewidth',1.5,'displayname',dispname);
    hold off;
    if logplot
        set(gca,'XScale','log');
    end
    xlim([20,fs/2]);
    xlabel('Frequency');
    ylabel('[dB]');
    grid on;
    legend('location','best');
    title('Audiogram vs. Approximations');
    
    subplot(2,1,2);
    hold on;
    if plotref
        plot(fftfreqs,zeros(size(fftfreqs)),'k','linewidth',1.5,'displayname','Audiogram (Reference)');
    end
    plot(fftfreqs,error_tf,lineformat,'linewidth',1.5,'displayname',dispname)
    hold off;
    if logplot
        set(gca,'XScale','log');
    end
    xlim([20,fs/2])
    ylim([-15,15]);
    xlabel('Frequency');
    ylabel('[dB]');
    grid on;
    legend('location','best');
    title('dB Error');
    
    
end


