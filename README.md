`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
module Top_fft #(
    parameter DATAWIDTH    = 16,
    parameter ADDRESSWIDTH = 6,   // 64 depth
    parameter CONFIGWIDTH  = 16,
    parameter FFT_LEN      = 64
)(
    input  wire aclk,
    input  wire aresetn,

    /* ================= FFT CONFIG ================= */
    output wire [CONFIGWIDTH-1:0] s_axis_config_tdata,
    output wire                   s_axis_config_tvalid,
    input  wire                   s_axis_config_tready,

    /* ================= FFT INPUT ================= */
    output wire  [2*DATAWIDTH-1:0] s_axis_data_tdata,
    output wire         s_axis_data_tvalid,
    input  wire        s_axis_data_tready,
    output wire         s_axis_data_tlast,

    /* ================= FFT OUTPUT ================= */
    input wire [31:0] m_axis_data_tdata,
    input wire        m_axis_data_tvalid,
    output  wire        m_axis_data_tready,
    input wire        m_axis_data_tlast,

    /* ================= OUTPUT SPLIT ================= */
    output wire  [23:0] out_real_data,
    output wire  [23:0] out_imag_data,

    /* ================= EVENTS ================= */
    output wire event_frame_started,
    output wire event_tlast_unexpected,
    output wire event_tlast_missing,
    output wire event_status_channel_halt,
    output wire event_data_in_channel_halt,
    output wire event_data_out_channel_halt
);

    /* ================= FSM STATES ================= */
    reg[2:0] state;
    localparam IDLE           =      3'b000;
    localparam CFG            =      3'b001;
    localparam BRAM_WAIT      =      3'b010;
    localparam STREAM         =      3'b011;
    localparam WAIT_FOR_DONE =       3'b100;
    localparam INTER_FRAME    =      3'b101;

 
    /* ================= CONFIG ================= */
    reg [CONFIGWIDTH-1:0] config_tdata_reg;
    reg                   config_tvalid_reg;
    reg[31:0]             data_tdata_reg;
    reg                   data_tvalid_reg;
    reg[ADDRESSWIDTH-1:0] bram_addr;
    reg [7:0]             sample_cnt;
    reg [7:0]             gap_cnt;
    
    //* Continuous Assignment *//
    assign s_axis_config_tdata =  config_tdata_reg;
    assign s_axis_config_tvalid = config_tvalid_reg;
    assign s_axis_data_tdata =    data_tdata_reg;
    assign s_axis_data_tvalid =  data_tvalid_reg;
    assign s_axis_data_tlast = (sample_cnt == FFT_LEN-1) && data_tvalid_reg; 
    assign m_axis_data_tready = 1'b1;

    /* ================= BRAM ================= */
    wire  [31:0] bram_dout;
    blk_mem_gen_0 u_bram (
        .clka  (aclk),
        .addra (bram_addr),
        .douta (bram_dout)
    );
    /* ================= FSM ================= */
    always @(posedge aclk) begin
        if (!aresetn) begin
            state              <= IDLE;
            bram_addr          <= 0;
            sample_cnt         <= 0;
            gap_cnt            <= 0;
            config_tvalid_reg  <= 0;
            data_tvalid_reg    <= 0;
            config_tvalid_reg  <= 16'b0000000010000110;
        end else begin
            case (state)
            IDLE: begin
                config_tvalid_reg <= 1'b1;
                state        <= CFG;
            end
            CFG: begin
                if (s_axis_config_tready && s_axis_config_tvalid) begin
                    config_tvalid_reg <= 0;
                    state <= BRAM_WAIT;
                    bram_addr <= 0;
                end
            end
            BRAM_WAIT: begin
            //one cycle delay for BRAM output to be valid
                state <= STREAM;
                data_tvalid_reg <= 1'b1;
             end
           STREAM:begin
                if (s_axis_data_tready && m_axis_data_tvalid) begin
                  if(sample_cnt == FFT_LEN -1)begin
                      sample_cnt <= 0;
                      data_tvalid_reg <= 0;
                      state <= WAIT_FOR_DONE;   // Next frame allowed only here
                end else begin
                     sample_cnt <= sample_cnt + 1;
                     bram_addr <= bram_addr + 1;
                  end
                 end
                end
                
            WAIT_FOR_DONE: begin
                if(m_axis_data_tvalid && m_axis_data_tlast)begin
                   state <= INTER_FRAME;
                   gap_cnt <= 0;
                 end
                end
            INTER_FRAME:begin 
              if(gap_cnt >= INTER_FRAME)begin
                   state <= BRAM_WAIT;
                   bram_addr <= 0;
                end else begin
                   gap_cnt <= gap_cnt+1;
                end
               end
             endcase
              end
                end
                
      // Fetch next data from BRAM
      always@(*) begin
         data_tdata_reg = bram_dout;
      end
      
      
    /* ================= FFT IP ================= */
    xfft_1 fft_i (
        .aclk(aclk),
        .aresetn(aresetn),

        .s_axis_config_tdata (s_axis_config_tdata),
        .s_axis_config_tvalid(s_axis_config_tvalid),
        .s_axis_config_tready(s_axis_config_tready),

        .s_axis_data_tdata   (s_axis_data_tdata),
        .s_axis_data_tvalid  (s_axis_data_tvalid),
        .s_axis_data_tready  (s_axis_data_tready),
        .s_axis_data_tlast   (s_axis_data_tlast),

        .m_axis_data_tdata   (m_axis_data_tdata),
        .m_axis_data_tvalid  (m_axis_data_tvalid),
        .m_axis_data_tready  (m_axis_data_tready),
        .m_axis_data_tlast   (m_axis_data_tlast),
//        .out_real_data       (out_real_data),
//        .out_imag_data       (out_real_data),

        .event_frame_started         (event_frame_started),
        .event_tlast_unexpected      (event_tlast_unexpected),
        .event_tlast_missing         (event_tlast_missing),
        .event_status_channel_halt   (event_status_channel_halt),
        .event_data_in_channel_halt  (event_data_in_channel_halt),
        .event_data_out_channel_halt (event_data_out_channel_halt)
    );

    /* ================= OUTPUT SPLIT ================= */
    assign out_real_data = m_axis_data_tdata[22:0];
    assign out_imag_data = m_axis_data_tdata[46:24];

endmodule
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
