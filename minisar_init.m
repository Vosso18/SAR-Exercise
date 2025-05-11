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
