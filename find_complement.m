function hout = find_complement(h)
    centeridx = ceil((numel(h)-1)/2+1);
    hdly = zeros(size(h));
    hdly(centeridx) = 1;
    hout = hdly - h;
end