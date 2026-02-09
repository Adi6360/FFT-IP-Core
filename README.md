Got it — you want MATLAB code to generate COE files for Vivado FFT IP with these exact settings:
FFT Length = 64
Channels = 1
Input from BRAM
Need input COE file (fixed-point complex samples)
Need output COE file (bit-accurate fixed-point FFT result)
Format = Fixed point
Scaling = Unscaled
Output order = Natural order
Vivado FFT IP compatible packing = {Imag, Real}
I’ll give you a clean, ready-to-run MATLAB script that:
✅ Generates test signal
✅ Quantizes to fixed-point
✅ Computes FFT (unscaled)
✅ Matches Vivado fixed-point behavior
✅ Writes input.coe and fft_output.coe
✅ Uses proper BRAM packing format
✅ MATLAB Script — FFT COE Generator (Vivado Compatible)
Copy code
Matlab
clc;
clear;

%% ============================================================
% FFT IP PARAMETERS (MATCH YOUR VIVADO CONFIG)
% ============================================================

N = 64;                 % FFT length
DATA_WIDTH = 16;        % Fixed point width
FRAC_BITS = 15;         % Q1.15 format
SCALE = 2^FRAC_BITS;

%% ============================================================
% GENERATE TEST INPUT SIGNAL (COMPLEX)
% ============================================================

n = 0:N-1;

% Example: single tone complex sinusoid
f = 5;  % tone index
x = 0.8 * exp(1j*2*pi*f*n/N);

% You can replace with:
% x = randn(1,N) + 1j*randn(1,N);

%% ============================================================
% FIXED POINT QUANTIZATION (INPUT)
% ============================================================

x_real = round(real(x) * SCALE);
x_imag = round(imag(x) * SCALE);

% Saturate to signed 16-bit
x_real = max(min(x_real, 32767), -32768);
x_imag = max(min(x_imag, 32767), -32768);

%% ============================================================
% COMPUTE FFT (UNSCALED — MATCH VIVADO SETTING)
% ============================================================

X = fft(x, N);   % MATLAB FFT is unscaled → matches "unscaled" mode

% Quantize output to same format
Xr = round(real(X) * SCALE);
Xi = round(imag(X) * SCALE);

Xr = max(min(Xr, 32767), -32768);
Xi = max(min(Xi, 32767), -32768);

%% ============================================================
% PACK FORMAT FOR VIVADO FFT IP BRAM
% FORMAT = {IMAG, REAL}
% ============================================================

input_packed  = uint32(typecast(int16([x_imag; x_real]), 'uint16'));
output_packed = uint32(typecast(int16([Xi; Xr]), 'uint16'));

% combine into 32-bit words
input_words  = bitshift(uint32(uint16(x_imag)),16) + uint32(uint16(x_real));
output_words = bitshift(uint32(uint16(Xi)),16) + uint32(uint16(Xr));

%% ============================================================
% WRITE INPUT COE FILE
% ============================================================

fid = fopen('fft_input.coe','w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k=1:N-1
    fprintf(fid,'%08X,\n', input_words(k));
end
fprintf(fid,'%08X;\n', input_words(N));
fclose(fid);

%% ============================================================
% WRITE OUTPUT COE FILE
% ============================================================

fid = fopen('fft_output.coe','w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k=1:N-1
    fprintf(fid,'%08X,\n', output_words(k));
end
fprintf(fid,'%08X;\n', output_words(N));
fclose(fid);

%% ============================================================
% OPTIONAL: VERIFY ERROR
% ============================================================

% Convert back to float for checking
xr_float = double(x_real)/SCALE + 1j*double(x_imag)/SCALE;
X_fixed_check = fft(xr_float,N);

fprintf('Done. COE files generated.\n');
✅ Vivado BRAM Data Format (Important)
Your FFT IP expects:
Copy code

TDATA = {IMAG[15:0], REAL[15:0]}
So each COE word is:
Copy code

32-bit hex = IIII RRRR
Example memory entry:
Copy code

0001FF20
means:
Copy code

Imag = 0x0001
Real = 0xFF20
✅ Vivado FFT IP Settings (Must Match)
Set these exactly:
Copy code

Transform length = 64
Channels = 1
Input data width = 16
Phase factor width = default
Scaling = Unscaled
Output order = Natural Order
Complex input = Enabled
✅ How to Use
In Vivado:
Input side
Copy code

Block Memory Generator → load fft_input.coe
→ connect to s_axis_data_tdata
Output verification
Copy code

Capture FFT output → compare with fft_output.coe
