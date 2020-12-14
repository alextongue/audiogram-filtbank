function C = gram2gain_old(B,compgains,compfreqs,fs)

    assert(numel(compgains)==numel(compfreqs), 'number of gains ~= number of frequencies');
    infl_thresh = -60; % influence threshold

    maxlen = 0;
    for ii = 1:numel(B)
        if numel(B(ii).ir) > maxlen
            maxlen = numel(B(ii).ir);
        end
    end
    
    fftlen = 2^ceil(log2(maxlen));
    fftfreqs = ((0:(fftlen-1))*fs/fftlen);
    
    fftfreqs = fftfreqs(1:fftlen/2);
    
    matchidx = zeros(size(compfreqs));
    for ii = 1:numel(compfreqs)
        diff = fftfreqs - compfreqs(ii);
        [~,idx] = min(diff.^2);
        matchidx(ii) = idx; % store dft idx nearest to audiogram freq
    end
    
    contrib = zeros(numel(B),numel(compfreqs));
    for ii = 1:numel(B)
        B(ii).tf = fft(B(ii).ir, fftlen);
        for jj = 1:numel(compfreqs)
            infl_temp = (abs(B(ii).tf(matchidx(jj))));
            contrib(ii,jj) = infl_temp;
        end
    end
    if ~isempty(infl_thresh)
        contrib(20*log10(contrib) < infl_thresh) = 0;
    end
    contribsums = sum(contrib,2);

    Cgains = contrib;
    for ii = 1:numel(B)
        Cgains(ii,:) = contrib(ii,:).*compgains;
    end
    
    Cgains = sum(Cgains,2) ./ contribsums;
    Cgains(isnan(Cgains)) = 0;
    
    C = B;
    for ii = 1:numel(B)
        C(ii).ir = B(ii).ir .* 10^(Cgains(ii)/20);
        C(ii).tf = B(ii).tf .* 10^(Cgains(ii)/20);
    end
    
end