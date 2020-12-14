function C = gram2gain(B,compgains,compfreqs,fs,opt_method,thresh)
% (! iterative method assumes unimodal or monotonic (LP/HP/BP shape) of filter bank !)

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
    if opt_method==1 % method 1: relative influence
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

    elseif opt_method==2 || opt_method==3 % 2=brute force; 3=brute w/slope
        Cgains = zeros(numel(B),1);
        agram_interp = interp1(compfreqs,compgains,fftfreqs,[],'extrap');
        boost_upper_init = 10^(100/20); % linear
        boost_lower_init = 10^(-10/20); % linear
        converg_thresh = 1e-8;

        for ii = 1:numel(B)
            pb_idxs = 20*log10(abs(B(ii).tfhalf)) > thresh; % passband threshold
            B(ii).cropfreqs = fftfreqs(pb_idxs);
            B(ii).crop = B(ii).tf(pb_idxs);
            agram_crop = agram_interp(pb_idxs);
            
            if opt_method==3 % SLOPE DETECTION (assumes unimodal)
                dy = agram_crop(end) - agram_crop(1);
                dydx = dy/numel(agram_crop);
                halfidx = round(numel(agram_crop)/2);
                if dy > 0 % positive slope == discard RHS
                    croprange = 1:halfidx;
                else % negative/constant slope == discard LHS
                    croprange = (numel(agram_crop)-halfidx+1):numel(agram_crop);
                end
                if (ii>1 && ii<8) && abs(dydx>0.5)
                    fprintf('REACHED SLOPE %.4f\n',dydx);
                    B(ii).cropfreqs = B(ii).cropfreqs(croprange);
                    B(ii).crop = B(ii).crop(croprange);
                    agram_crop = agram_crop(croprange);
                end
                
            end
            error_ref = 20*log10(abs( B(ii).crop )) - agram_crop; % initial error
            errsum_ref = sum(error_ref)./numel(error_ref);
            boost_upper = boost_upper_init;
            boost_lower = boost_lower_init;
            boost_log = zeros(100,1);
            errsum_log = zeros(100,1);
            iter = 0;
            while boost_upper-boost_lower > converg_thresh
                iter = iter+1;
                boost_log(iter) = (boost_upper + boost_lower) / 2;
                Btemp = B(ii).crop * boost_log(iter);
%                    fprintf('iter %d: trying %.2fdB boost...',iter,20*log10(boost_log(iter)));
                error_curr = 20*log10(abs( Btemp )) - agram_crop;
                errsum_curr = sum(error_curr)./numel(error_curr);
               %fprintf('error=%.4fdB\n',errsum_curr);
                if errsum_curr > 0 %overshoot condition
                    boost_upper = boost_log(iter);
                else %undershoot condition
                    boost_lower = boost_log(iter);
                end
                errsum_log(iter) = errsum_ref;
                errsum_ref = errsum_curr;
            end
            Cgains(ii) = boost_log(iter);
            fprintf('C%d : iter=%d, boost=%.2fdB, error=%.4f\n',...
                ii, iter, 20*log10(boost_log(iter)), errsum_ref)
        end
        errsum_log = [errsum_log(1:(iter));errsum_curr];

        figure;
        subplot(1,2,1)
        plot(1:(iter), 20*log10(abs(boost_log(1:(iter)))),'k', 'linewidth',1.5);
        grid on;
        xlabel('Iteration');
        ylabel('[dB]');
        title('Gain v. Iteration');
        subplot(1,2,2)
        plot(0:(iter), errsum_log,'k', 'linewidth',1.5);
        grid on;
        xlabel('Iteration');
        title('MMSE v. Iteration');
        

    end
        
    % apply gains
    C = B;
    for ii = 1:numel(B)
        C(ii).ir = B(ii).ir .* Cgains(ii);
        C(ii).tf = B(ii).tf .* Cgains(ii);
    end
    
end