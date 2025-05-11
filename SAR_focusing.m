

% --------------- Set parameters and open the file ------------------------

% Define system (and data) parameters
sys = struct( ...
  'fc',       5.3e9,     ...  % Carrier frequency   [Hz]
  'fs',     100.0e6,     ...  % Sampling frequency  [Hz]
  'del',      7.2277e-5, ...  % Range offset        [s] 
  'pule',    17.84e-6,   ...  % Pulse length        [s]
  'pubw',   100e6,       ...  % Pulse bandwidth     [s]
  'azdx',     1.5,       ...  % Azimuth spacing     [m]
  'alt',   7500,         ...  % Flight altitude     [m]
  'nra',   2048,         ...  % Data samples (range)
  'naz',   1024          ...  % Data lines   (azimuth)
);

% Read raw SAR data - I and Q 8-bit samples
fid = fopen('fl067_m0953_foulum_chh.win');
iqdata = fread(fid, [2*sys.nra,sys.naz], 'int8');
fclose(fid);

% Convert I/Q data to complex numbers
raw = complex(iqdata(1:2:end,:), iqdata(2:2:end,:));

% Clean-up
clear fid iqdata 

%%
%----------------------- SAR exercise 1 -----------------------------------
%--------------------- Range Compression ----------------------------------

% Making a chirp for the matched filter.
B = sys.pubw;
T = sys.pule;
f_s = sys.fs;

N = round(T * f_s/2);
M = 2 * N + 1; % Length of the chirp.
chirp_rate = B / T;
delta_t = 1 / f_s;
t = (-N : N)/f_s;
chirp_phase = pi * chirp_rate * t.^2;

% This is our chirp.
u_range = exp(1i * chirp_phase);

% Padding the chirp with zeros.
w_range = [u_range zeros(1, sys.nra - M)];

% It needs to be transposed to be a column vector,
% avoids having to do a for loop to multiply each column.
w_range = w_range.';

% Correlate columns with filter, use fft and ifft.
H_ra = fft(w_range);
raw_range_compressed = ifft(fft(raw).* conj(H_ra));

% Delete rows corresponding to throwaway region.
% Discard the last M minus 1 samples.
% The extra minus 1 is to discard row 265.
discard_ra = M - 1;
raw_range_compressed(size(raw, 1) - (discard_ra - 1):end , :) = [];


%%
%-------------------------- SAR exercise 2 --------------------------------
%------------------------ Azimuth Compression -----------------------------
lambda = 3e8 / sys.fc;
range_offset = sys.del*3e8 / 2;
azimuth_spacing = sys.azdx;
range_pixel_spacing = 3e8 / (2 * sys.fs);

% The r_0 is changing for each row of the matrix so we have to do
% a loop where each row is correlated with the filter.

% Declare a matrix to put the results in.
raw_azimuth_compressed = zeros(size(raw_range_compressed));

for i = 1:size(raw_range_compressed, 1)
    r_0 = range_offset + (i - 1) * range_pixel_spacing;
    L = (lambda * r_0) / (2 * azimuth_spacing);
    
    % We need an even number of points on each side of zero. So an odd
    % number in total.
    x = -L/2 : azimuth_spacing : L/2;
    chirp_rate = -2 / (lambda * r_0);
    chirp_phase = pi * chirp_rate * x.^2;
    
    % This is our chirp.
    u_azimuth = exp(1i * chirp_phase);

    % Padding the chirp with zeros. The chirp should wrap around.
    % So there are: number of rows - P zeros in the middle.
    N = floor(L / azimuth_spacing / 2);
    P = 2 * N + 1;
    w_azimuth = [u_azimuth(N + 1:P) zeros(1, sys.naz - P) u_azimuth(1:N)];

    % Correlate rows with filter.
    H_az = fft(w_azimuth);
    raw_azimuth_compressed(i,:) = ifft(fft(raw_range_compressed(i,:)).* conj(H_az));
end

% Discard the first M and the last M samples.
% The extra minus 1 is to discard column 885.
discard_az = N;
raw_azimuth_compressed(:, 1:N) = [];
raw_azimuth_compressed(:, size(raw_azimuth_compressed, 2) - (discard_az - 1):end) = [];


%%
%--------------------------- SAR exercise 3 -------------------------------
%------------------------- FFT interpolation ------------------------------

% Isolate patch that includes the central peak and a few of the side-lobes
% of the corner reflector.

peak_raw = max(max(raw_azimuth_compressed));
[peak_range_raw, peak_azimuth_raw]=find(raw_azimuth_compressed==peak_raw);

% Upsampling factor.
up_sampling_factor = 64;

% How much do we want on each side of the peak.
range_delta = peak_range_raw - up_sampling_factor:peak_range_raw + up_sampling_factor;
azimuth_delta = peak_azimuth_raw - up_sampling_factor:peak_azimuth_raw + up_sampling_factor;

patch = raw_azimuth_compressed(range_delta,azimuth_delta);

% Do 2D interpolation by first doing a 1D interpolation in the range
% dimension and then in the azimuth. Use interpft.
% The upsampling should be by a known factor. To upsample by a factor 8
% find the size of the patch in the dimension you are upsampling and
% multiply it by 8.
N_range = size(patch, 1);
N_azimuth = size(patch, 2);
interpolation_1D = interpft(patch, N_range * up_sampling_factor, 1);
interpolation_2D = interpft(interpolation_1D, N_azimuth * up_sampling_factor, 2);

% Take range and azimuth profiles of the power across the peak.
% Find the index of the max value.
peak = max(max(interpolation_2D));
[peak_range, peak_azimuth]=find(interpolation_2D==peak);


%%
%--------------------------------------------------------------------------
%---------------------------- Plotting ------------------------------------
%--------------------------------------------------------------------------

% This is just a plot to show the range compressed image before cutting the
% throw-away region.
figure();
imagesc(abs(ifft(fft(raw).* conj(H_ra))));

% title({'Range compressed', 'including throw-away'},FontSize=18);
xlabel('Azimuth (samples)', FontSize=18);
ylabel('Range (samples)', FontSize=18);
fontsize(12, "points");
set(gcf, 'Position', [0, 200, 512, 1024]);


%------------------- Plotting the Range Compression -----------------------

% Set up plot.

% Plot amplitude.
figure();
imagesc(abs(real(raw)));

% title('Raw data amplitude',FontSize=18);
xlabel('Azimuth (samples)', FontSize=18);
ylabel('Range (samples)', FontSize=18);
fontsize(16, "points");
set(gcf, 'Position', [0, 200, 600, 800]);


% Plot phase.
% t2 = nexttile(2, [2 1]);
% imagesc(angle(raw));

% title('Raw data phase',FontSize=18);
% xlabel('Azimuth (samples)', FontSize=18);
% ylabel('Range (samples)', FontSize=18);
% t2.FontSize = 20;

% Plot range compressed image.
figure();

imagesc(abs(raw_range_compressed));

% title('Range compressed',FontSize=18);
xlabel('Azimuth (samples)', FontSize=18);
ylabel('Range (samples)', FontSize=18);
fontsize(16, "points");
set(gcf, 'Position', [200, 200, 300, 400])


%------------ Plotting the Range and azimuth Compressed image -------------

% Plot azimuth compressed.
figure();

imagesc(abs(raw_azimuth_compressed).^0.3);

% title('Range and azimuth compressed',FontSize=18);
xlabel('Azimuth (samples)', FontSize=18);
ylabel('Range (samples)', FontSize=18);
fontsize(16, "points");
set(gcf, 'Position', [600, 400, 884, 264])


%------------------ Plotting interpolation --------------------------------

% Plot 2D interpolated image.
figure();

imagesc(abs(interpolation_2D));

title('2D interpolated',FontSize=18);
xlabel('Azimuth (samples)', FontSize=18);
ylabel('Range (samples)', FontSize=18);
fontsize(16,"points");


% Plot range and azimuth power profiles.
figure();
tiledlayout(1,2);
set(gcf, 'Position', get(0, 'Screensize'))

% Range power profile.
t1 = nexttile;

% The column we are interested in plotting.
range_power = interpolation_2D(:,peak_azimuth);

% The x-axis.
range_pixels = ((-peak_range + 1):(length(range_power) - peak_range))./up_sampling_factor;

% The y-values, normalized to get db.
y_azimuth =  10*log10(abs(range_power./peak).^2);

plot(range_pixels, y_azimuth);

title('Range power profile',FontSize=24);
xlabel('Range distance from peak (original pixels)', FontSize=20);
ylabel('Power (dB)', FontSize=20);
t1.FontSize = 20;
grid on;

% Plot azimuth power profile.
t2 = nexttile;

% The row we are interested in plotting.
azimuth_power = interpolation_2D(peak_range,:);

% The x-axis.
azimuth_pixels = ((-peak_azimuth + 1):(length(range_power) - peak_azimuth))./up_sampling_factor;

% The y-values, normalized to get db.
y_azimuth =  10*log10(abs(azimuth_power./peak).^2);

plot(azimuth_pixels, y_azimuth);

title('Azimuth power profile',FontSize=24);
xlabel('Azimuth distance from peak (original pixels)', FontSize=20);
ylabel('Power (dB)', FontSize=20);
t2.FontSize = 20;
grid on

axis(t1, [-8 8 -60 10]);
axis(t2, [-8 8 -60 10]);