clearvars; clc;

%% params

agram_freqs = [250,500,1000,2000,4000,8000];
agram_data = [-10, -15, -22, -40, -50, -65]; % monotonic slope
%agram_data = [-15, -26, -35, -55, -95, -65]; % complex slope
%agram_data = [0, -10, -20, -30, -40, -50]; % 10dB/oct falling threshold
%agram_data = [-50, -40, -30, -20, -10, 0]; % 10dB/oct rising threshold
%agram_data = [-20,-20,-20,-20,-20,-20]; % constant
%agram_data = [-20,0,-20,0,-20,0]; % zigzag

agram_comp = 0 - agram_data; % compensation

fs = 16000;

H_order = 18;
F_order = 38;

H_bw = 0.25*fs; % bandwidth (Hz)
F_bw = 0.125*fs; % bandwidth (Hz)

% M = 2; % dyadic
% H_f0 = 0.45*fs; % passband edge


%% File input.

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

%% Derive filter banks.
[B,P] = make_nonunif_fb(fs,H_order,H_bw,F_order,F_bw);
Bsum = (B(1).ir+B(2).ir+B(3).ir+B(4).ir+B(5).ir+B(6).ir+B(7).ir+B(8).ir);

%% add gains

% hand-tuned gains
%{
C(1).ir = 10^(audiogram_comp(1)/20)*B(1).ir;
C(2).ir = 10^((audiogram_comp(3)-4)/20)*B(2).ir;
C(3).ir = 10^((audiogram_comp(3)+3)/20)*B(3).ir;
C(4).ir = 10^(45/20)*B(4).ir;
C(5).ir = 10^(53/20)*B(5).ir;
C(6).ir = B(6).ir;
C(7).ir = B(7).ir;
C(8).ir = 10^(65/20)*B(8).ir;
%}

% interp (sarah)
agram_freqs2 = [250 500 1000 2000 4000 6000 7000 8000];
agram_comp2 = audiogramMatch(agram_comp);
C1 = gram2gain(B,agram_comp2,agram_freqs2,fs,1,-60); % optimization method 1
C1sum = (C1(1).ir+C1(2).ir+C1(3).ir+C1(4).ir+C1(5).ir+C1(6).ir+C1(7).ir+C1(8).ir);

% iterative
C2 = gram2gain(B,agram_comp,agram_freqs,fs,2,-3); % optimization method 2
C2sum = (C2(1).ir+C2(2).ir+C2(3).ir+C2(4).ir+C2(5).ir+C2(6).ir+C2(7).ir+C2(8).ir);

% iterative w/slope detection
C3 = gram2gain(B,agram_comp,agram_freqs,fs,3,-6); % optimization method 2
C3sum = (C3(1).ir+C3(2).ir+C3(3).ir+C3(4).ir+C3(5).ir+C3(6).ir+C3(7).ir+C3(8).ir);



%% process
Y = B;
for ii = 1:numel(C1)
    Y(ii).ir = conv(x,C1(ii).ir);
end
Ysum = (Y(1).ir+Y(2).ir+Y(3).ir+Y(4).ir+Y(5).ir+Y(6).ir+Y(7).ir+Y(8).ir);

%%
% Prototype filters
%{
fh.fig1 = figure();
[~] = plotfig(H, fs, 'mag', fh.fig1, 'H(z)');
[~] = plotfig(F, fs, 'mag', fh.fig1, 'F(z)');

fh.fig2 = figure('name', 'Prototype and derived filters H');
[~] = plotfig(H, fs, 'mag', fh.fig2, 'H(z)');
[~] = plotfig(H2, fs, 'mag', fh.fig2, 'H(z^2)');
[~] = plotfig(H4, fs, 'mag', fh.fig2, 'H(z^4)');
[~] = plotfig(H8, fs, 'mag', fh.fig2, 'H(z^8)');

fh.fig3 = figure('name', 'Proto and derived filters F');
[~] = plotfig(F, fs, 'mag', fh.fig3, 'F(z)');
[~] = plotfig(F2, fs, 'mag', fh.fig3, 'F(z^2)');
[~] = plotfig(F4, fs, 'mag', fh.fig3, 'F(z^4)');
[~] = plotfig(F8, fs, 'mag', fh.fig3, 'F(z^8)');

fh.fig4 = figure('name', 'Filter banks pre-subtraction');
[~] = plotfig(P(1).ir, fs, 'mag', fh.fig4, 'P(1).ir');
[~] = plotfig(P(2).ir, fs, 'mag', fh.fig4, 'P(2).ir');
[~] = plotfig(P(3).ir, fs, 'mag', fh.fig4, 'P(3).ir');
[~] = plotfig(P(4).ir, fs, 'mag', fh.fig4, 'P(4).ir');
[~] = plotfig(P(5).ir, fs, 'mag', fh.fig4, 'P(5).ir');
[~] = plotfig(P(6).ir, fs, 'mag', fh.fig4, 'P(6).ir');
[~] = plotfig(P(7).ir, fs, 'mag', fh.fig4, 'P(7).ir');
[~] = plotfig(P(8).ir, fs, 'mag', fh.fig4, 'P(8).ir');
%}

% Flat filterbanks (B)

fh.fig1 = figure('name', 'Filter banks');
[~] = plotfig(B(1).ir, fs, 'mag', fh.fig1, 'B1','k');
[~] = plotfig(B(2).ir, fs, 'mag', fh.fig1, 'B2','k');
[~] = plotfig(B(3).ir, fs, 'mag', fh.fig1, 'B3','k');
[~] = plotfig(B(4).ir, fs, 'mag', fh.fig1, 'B4','k');
[~] = plotfig(B(5).ir, fs, 'mag', fh.fig1, 'B5','k');
[~] = plotfig(B(6).ir, fs, 'mag', fh.fig1, 'B6','k');
[~] = plotfig(B(7).ir, fs, 'mag', fh.fig1, 'B7','k');
[~] = plotfig(B(8).ir, fs, 'mag', fh.fig1, 'B8','k');
hold on;
plot(agram_freqs,zeros(size(agram_freqs)),'r', 'linewidth', 1.5, 'displayname','Sum');
plot(agram_freqs,zeros(size(agram_freqs)),'kx', 'linewidth', 1.5, 'displayname','Sum');
hold off;
legend hide;
title('Filter Banks, flat response');
%}

%%
fh.fig2 = figure('name', 'Matching Errors');
[~] = ploterror(agram_freqs, agram_comp, C1sum, fs, fh.fig2, 'Interpolated', 'r:', true);
[~] = ploterror(agram_freqs, agram_comp, C2sum, fs, fh.fig2, 'Iterative SE', 'b:', false);
[~] = ploterror(agram_freqs, agram_comp, C3sum, fs, fh.fig2, 'Iterative SE w/slope', 'g:', false);

%%
fh.fig5 = figure('name', 'Filter banks');
[~] = plotfig(C3(1).ir, fs, 'maglog', fh.fig5, 'C1', 'k');
[~] = plotfig(C3(2).ir, fs, 'maglog', fh.fig5, 'C2', 'k');
[~] = plotfig(C3(3).ir, fs, 'maglog', fh.fig5, 'C3', 'k');
[~] = plotfig(C3(4).ir, fs, 'maglog', fh.fig5, 'C4', 'k');
[~] = plotfig(C3(5).ir, fs, 'maglog', fh.fig5, 'C5', 'k');
[~] = plotfig(C3(6).ir, fs, 'maglog', fh.fig5, 'C6', 'k');
[~] = plotfig(C3(7).ir, fs, 'maglog', fh.fig5, 'C7', 'k');
[~] = plotfig(C3(8).ir, fs, 'maglog', fh.fig5, 'C8', 'k');
[~] = plotfig(C3sum, fs, 'maglog', fh.fig5, 'Sum', 'g');
hold on;
% plot(audiogram_freqs,audiogram_comp,'kx:', 'linewidth', 1.5); % alex
plot(agram_freqs,agram_comp,'kx:', 'linewidth', 1.5); % sarah
hold off;
ylim([10,100]);
legend hide;
title('Audiogram Match (Iterative Optimization)');

%%
%{
fh.fig6 = figure('name', 'Filter banks IR');
[~] = plotfig(B(1).ir, fs, 'phase', fh.fig6, 'B(1).ir',[]);
[~] = plotfig(B(2).ir, fs, 'phase', fh.fig6, 'B(2).ir',[]);
[~] = plotfig(B(3).ir, fs, 'phase', fh.fig6, 'B(3).ir',[]);
[~] = plotfig(B(4).ir, fs, 'phase', fh.fig6, 'B(4).ir',[]);
[~] = plotfig(B(5).ir, fs, 'phase', fh.fig6, 'B(5).ir',[]);
[~] = plotfig(B(6).ir, fs, 'phase', fh.fig6, 'B(6).ir',[]);
[~] = plotfig(B(7).ir, fs, 'phase', fh.fig6, 'B(7).ir',[]);
[~] = plotfig(B(8).ir, fs, 'phase', fh.fig6, 'B(8).ir',[]);
%}
%%
fh.fig4 = figure('name','Signal In/Out');
[~] = plotfig(x, fs, 'maglog', fh.fig4, 'X', []);
[~] = plotfig(Ysum, fs, 'maglog', fh.fig4, 'Y', []);
title('Signal In/Out');
ylim([-20,100]);