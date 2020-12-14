function [B,P] = make_fb_nonunif(fs, H_order, H_bw, F_order, F_bw)
    % prototype filters are lin. phase (grpdelay == N/2)
    Hproto = dsp.FIRHalfbandDecimator(...
        'Specification','Filter order and transition width', ...
        'FilterOrder',H_order, ...
        'TransitionWidth',H_bw, ...
        'SampleRate',fs);

    Fproto = dsp.FIRHalfbandDecimator(...
        'Specification','Filter order and transition width', ...
        'FilterOrder',F_order, ...
        'TransitionWidth',F_bw, ...
        'SampleRate',fs);

    % upsample and find complement
    H = Hproto.coeffs.Numerator;
    Hc = find_complement(H);
    H2 = upsample(H,2);
    H4 = upsample(H,4);
    H8 = upsample(H,8);

    F = Fproto.coeffs.Numerator;
    Fc = find_complement(F);
    F2 = upsample(F,2);
    F4 = upsample(F,4);
    F8 = upsample(F,8);

    % derive P(1).ir-P(8).ir
    P(1).ir = conv(conv(conv(H8, F4), F2), F);
    P(8).ir = conv(conv(conv(H8, F4), F2), Fc);
    P(2).ir = conv(conv(H4,F2),F);
    P(7).ir = conv(conv(H4,F2),Fc);
    P(3).ir = conv(H2, F);
    P(6).ir = conv(H2, Fc);
    P(4).ir = H;
    P(5).ir = Hc;

    % add delays
    B(1).ir = P(1).ir;
    B(8).ir = P(8).ir;
    P(2).ir = match_lp_gd(P(2).ir,P(1).ir);
    P(7).ir = match_lp_gd(P(7).ir,P(8).ir);
    P(3).ir = match_lp_gd(P(3).ir,P(2).ir);
    P(6).ir = match_lp_gd(P(6).ir,P(7).ir);
    P(4).ir = match_lp_gd(P(4).ir,P(3).ir);
    P(5).ir = match_lp_gd(P(5).ir,P(6).ir);

    % subtract
    B(2).ir = P(2).ir - P(1).ir;
    B(7).ir = P(7).ir - P(8).ir;
    B(3).ir = P(3).ir - P(2).ir;
    B(6).ir = P(6).ir - P(7).ir;
    B(4).ir = P(4).ir - P(3).ir;
    B(5).ir = P(5).ir - P(6).ir;
end