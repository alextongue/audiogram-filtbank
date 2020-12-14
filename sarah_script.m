clearvars; clc;

%% params

audiogram_freqs = [250,500,1000,2000,4000,8000];
audiogram_data = [-10, -15, -22, -40, -50, -65];

audiogram_comp = audiogramMatch(0-audiogram_data);

fs = 16000;

H_order = 18;
F_order = 38;

H_bw = 0.25*fs; % bandwidth (Hz)
F_bw = 0.125*fs; % bandwidth (Hz)

% M = 2; % dyadic
% H_f0 = 0.45*fs; % passband edge


%% input
% filename
fname = fullfile(pwd, 'FB_male_female_single-talk_seq.wav');
frange = [1,48000];

% sinusoid params 
Nx = 0.25; % input length (sec)
xfreqs = [250,500,1000,2000,4000];

if isempty(fname)
    t = 0:(Nx/fs):(Nx-Nx/fs);
    n = (0:Nx*fs)./fs;

    x = zeros(size(n));
    for ii = 1:numel(xfreqs)
        x = x + sin(2*pi*xfreqs(ii).*n);
    end
    x = x./(2^(numel(xfreqs)+6));

    [~] = plotfig(x, fs, 'imp', [], 'input');
    [~] = plotfig(x, fs, 'mag', [], 'inputfft');
else
    [x,xfs] = audioread(fname,frange);
    if xfs>fs
        if ~mod(xfs,fs)
            x = downsample(x, xfs/fs);
            xfs = fs;
        else
            warning('cannot downsample input sig to fs=%d',fs);
        end
    end
end

%% make filts

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

%% add gains
C = gram2gain_old(B,audiogram_comp,[250 500 1000 2000 4000 6000 7000 8000],fs);

% C = B;
% for ii=1:numel(C)
%     C(ii).ir = B(ii).ir .* audiogram_comp(ii);
% end

Csum = (C(1).ir+C(2).ir+C(3).ir+C(4).ir+C(5).ir+C(6).ir+C(7).ir+C(8).ir);
Bsum = (B(1).ir+B(2).ir+B(3).ir+B(4).ir+B(5).ir+B(6).ir+B(7).ir+B(8).ir);

%% process
Y(1).ir = conv(x,C(1).ir);
Y(2).ir = conv(x,C(2).ir);
Y(3).ir = conv(x,C(3).ir);
Y(4).ir = conv(x,C(4).ir);
Y(5).ir = conv(x,C(5).ir);
Y(6).ir = conv(x,C(6).ir);
Y(7).ir = conv(x,C(7).ir);
Y(8).ir = conv(x,C(8).ir);

sumtest = (Y(1).ir+Y(2).ir+Y(3).ir+Y(4).ir+Y(5).ir+Y(6).ir+Y(7).ir+Y(8).ir);

%%
% ploterror(audiogram_freqs,audiogram_comp,Csum,fs,[]);
ploterror([250 500 1000 2000 4000 6000 7000 8000],audiogram_comp,Csum,fs,[]);

%%
fh.fig5 = figure('name', 'Filter banks');
[~] = plotfig(C(1).ir, fs, 'maglog', fh.fig5, 'C(1).ir');
[~] = plotfig(C(2).ir, fs, 'maglog', fh.fig5, 'C(2).ir');
[~] = plotfig(C(3).ir, fs, 'maglog', fh.fig5, 'C(3).ir');
[~] = plotfig(C(4).ir, fs, 'maglog', fh.fig5, 'C(4).ir');
[~] = plotfig(C(5).ir, fs, 'maglog', fh.fig5, 'C(5).ir');
[~] = plotfig(C(6).ir, fs, 'maglog', fh.fig5, 'C(6).ir');
[~] = plotfig(C(7).ir, fs, 'maglog', fh.fig5, 'C(7).ir');
[~] = plotfig(C(8).ir, fs, 'maglog', fh.fig5, 'C(8).ir');
[~] = plotfig(Csum, fs, 'maglog', fh.fig5, 'fbsum');
hold on;
plot([250 500 1000 2000 4000 6000 7000 8000],audiogram_comp,'kx:', 'linewidth', 1.5);
hold off;
ylim([0,100]);


%%

fh.fig4 = plotfig(x, fs, 'maglog', [], 'X');
[~] = plotfig(sumtest, fs, 'maglog', fh.fig4, 'Y');
ylim([-20,100]);