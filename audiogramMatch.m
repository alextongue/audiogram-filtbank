function gain = audiogramMatch(thresh)

    % initialize temporary gain vector
    temp_gain = ones(1, 6);
    % establish frequencies
    freqs = [250 500 1000 2000 4000 8000];
    % find thresholds > 15dB

    % Convert dB
    thresh_conv = 10.^(thresh/20);
    % Find midpoint between octave thresholds
    samples = [375, 750, 1500, 3000, 6000];
    vq = interp1(freqs, thresh_conv, samples);
%     temp_gain = [temp_gain(1) vq]; % for 15dB thresh
    temp_gain = [thresh_conv(1) vq];
    % Find gain at 6kHz and 7kHz via interpolation
    vq = interp1(freqs, temp_gain, [6000 7000]);
    new_gain = [temp_gain(1:5) vq temp_gain(6)];
    gain = 20*log10(new_gain);

    % 15dB thresh discard
%     idx = find(thresh <= 15); % threshold at 250Hz will be gain if below 15dB
%     gain(idx) = 1;
%     if ~ismember(1, idx)
%         gain(1) = thresh(1);
%     end
    
end