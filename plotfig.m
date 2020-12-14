function fhandleout = plotfig(plotdata, fs, domain, fhandlein, plotname, linetype)

    if isempty(fhandlein)
        fhandleout = figure;
    else
        fhandleout = figure(fhandlein);
    end
    
    
    % plot
    if strcmpi(domain, 'imp')
        plotx = (1:numel(plotdata))-1;
    end
    
    if strcmpi(domain, 'mag') ...
        || strcmpi(domain, 'maglog')
        plotx = (0:(numel(plotdata)-1)) .* (fs/numel(plotdata));
        plotdata = 20*log10(abs(fft( plotdata, numel(plotdata) )));
        %plotdata = AKfractOctSmooth(plotdata,'schaerer',fs,12,'amp','oct',[],false);
    end
    
    if strcmpi(domain, 'grp')
        plotx = (1:numel(plotdata))-1;
        %plotx = (0:(numel(plotdata)-1)) .* (fs/numel(plotdata));
        %plotdata = 20*log10(abs(fft( plotdata, numel(plotdata) )));
        plotdata = grpdelay(plotdata,1,numel(plotx));
    end

    if strcmpi(domain, 'phase')
        plotx = (0:(numel(plotdata)-1)) .* (fs/numel(plotdata));
        plotdata = unwrap(angle(fft( plotdata, numel(plotdata) )));
        %plotdata = grpdelay(plotdata,1,numel(plotx));
    end

    
    hold on;
    if ~isempty(linetype)
        plot(plotx, plotdata, linetype, 'linewidth', 1.5, 'displayname', plotname);
    else
        plot(plotx, plotdata, 'linewidth', 1.5, 'displayname', plotname);
    end
    hold off;
    
    % format
    if strcmpi(domain, 'mag')
        xlim([0,0.5]*fs);
        ylim([-35,5]);
        xlabel('Frequency');
        ylabel('Amplitude [dB]');
        grid on;
    end
    
    if strcmpi(domain, 'maglog')
        xlim([0.005,0.5]*fs);
        ylim([-20,80]);
        set(gca, 'XScale', 'log');
        xlabel('Frequency');
        ylabel('Amplitude [dB]');
        grid on;
    end
    legend('location','northwest');
    
end