function C = gram2gain(B,compgains,compfreqs,fs,opt_method,thresh)

    % Input params / prep.
    assert(numel(compgains)==numel(compfreqs), 'number of gains ~= number of frequencies');
    compgains = compgains(:);
    compfreqs = compfreqs(:);

    % Find FFT size.
    maxlen = 0;
    for ii = 1:numel(B)
        if numel(B(ii).ir) > maxlen
            maxlen = numel(B(ii).ir);
        end
    end
    fftlen = 2^ceil(log2(maxlen));
    fftfreqs = ((0:(fftlen-1))*fs/fftlen);
    fftfreqs = fftfreqs(1:fftlen/2); % single sided spectrum
    for ii = 1:numel(B)
        B(ii).tf = fft(B(ii).ir, fftlen);
        B(ii).tfhalf = B(ii).tf(1:(fftlen/2));
    end
    
    
    % Find optimal gains.
    Bgains = ones(numel(B),1);
    switch opt_method
        case 1 % method 1: relative influence
            % Find filterbanks that contribute to each audiogram frequency.
            matchidx = zeros(size(compfreqs));
            for ii = 1:numel(compfreqs)
                diff = fftfreqs - compfreqs(ii);
                [~,idx] = min(diff.^2);
                matchidx(ii) = idx; % store dft idx nearest to audiogram freq
            end
            contrib = zeros(numel(B),numel(compfreqs));
            for ii = 1:numel(B) % ***
                for jj = 1:numel(compfreqs)
                    contrib(ii,jj) = (abs(B(ii).tf(matchidx(jj)))); % linear magn.
                end
            end
            if ~isempty(thresh) % contribution discard threshold
                contrib(20*log10(contrib) < thresh) = 0;
            end
            Cgains = contrib;
            for ii = 1:numel(B)
                Cgains(ii,:) = contrib(ii,:).*compgains'; % contrib (lin/%) * compgains (dB)
            end
            contribsums = sum(contrib,2);
            Cgains = Cgains./contribsums; % normalize
            Cgains(isnan(Cgains)) = 0;
            Cgains = sum(Cgains,2);
            Cgains = 10.^(Cgains./20);
            
        case 2 % method 2: brute force
            Cgains = zeros(numel(B),1);
            agram_interp = interp1(compfreqs,compgains,fftfreqs,[],'extrap');
            boost_upper_init = 10^(100/20); % linear
            boost_lower_init = 10^(-100/20); % linear
            converg_thresh = 1e-8;
            
            for ii = 1:numel(B)
                pb_idxs = 20*log10(abs(B(ii).tfhalf)) > thresh; % passband threshold
                B(ii).cropfreqs = fftfreqs(pb_idxs);
                B(ii).crop = B(ii).tf(pb_idxs);
                agram_crop = agram_interp(pb_idxs);
                error_ref = 20*log10(abs( B(ii).crop )) - agram_crop; % initial error
                errsum_ref = sum(error_ref);
                boost_upper = boost_upper_init;
                boost_lower = boost_lower_init;
                iter = 1;
                while boost_upper-boost_lower > converg_thresh
                    boost_curr = (boost_upper + boost_lower) / 2;
                    Btemp = B(ii).crop * boost_curr;
%                    fprintf('iter %d: trying %.2fdB boost...',iter,20*log10(boost_curr));
                    error_curr = 20*log10(abs( Btemp )) - agram_crop;
                    errsum_curr = sum(error_curr);
%                    fprintf('error=%.4fdB\n',errsum_curr);
                    if errsum_curr > 0 %overshoot condition
                        boost_upper = boost_curr;
                    else %undershoot condition
                        boost_lower = boost_curr;
                    end
                    errsum_ref = errsum_curr;
                    iter = iter+1;
                end
                Cgains(ii) = boost_curr;
                fprintf('B%d : iter=%d, boost=%.2fdB, error=%.4f\n',...
                    ii, iter-1, 20*log10(boost_curr), errsum_ref)
            end
            
        otherwise
            error('Method ID not found');
    end
    
    % TODO: iterative; for each B, iteratively boost and query the
    % frequencies affected until minimum statistics is met.    
    % ITERATE BY: usable -6dB bandwidth? -40dB/thresh bandwidth? hard code?
    
    % apply gains
    C = B;
    for ii = 1:numel(B)
        C(ii).ir = B(ii).ir .* Cgains(ii);
        C(ii).tf = B(ii).tf .* Cgains(ii);
    end
    
end