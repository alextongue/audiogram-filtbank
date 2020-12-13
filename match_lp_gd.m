function h_out = match_lp_gd(h_short, h_long)
    

    N_short = numel(h_short)-1;
    N_long = numel(h_long)-1;
    [~,maxidx_short] = max(h_short);
    [~,maxidx_long] = max(h_long);
    maxidx_diff = maxidx_long-maxidx_short; 
    
    assert(N_long > N_short,'make sure first filter is shorter');
    
    h_out = zeros(size(h_long));
    inrange = (maxidx_diff+1):(maxidx_diff+numel(h_short));
    h_out(inrange) = h_short;
    
end