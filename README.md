# FFT-IP-Core
1. Updated Top_fft (BRAM-based, no sine generator)
What changed
in_real_sine / in_imag_sine / in_addr are removed from functionality
BRAM supplies {imag, real} samples
Internal FSM streams BRAM → FFT
FFT config logic kept in your SEND / DONE style
Top_fft.v
`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
module Top_fft #(
    parameter DATAWIDTH     = 16,
    parameter ADDRESSWIDTH  = 16,
    parameter CONFIGWIDTH   = 48,
    parameter FFT_LEN       = 64
)(
    input  wire aclk,
    input  wire aresetn,

    /* ================= FFT CONFIG ================= */
    output wire [CONFIGWIDTH-1:0] s_axis_config_tdata,
    output wire                   s_axis_config_tvalid,
    input  wire                   s_axis_config_tready,

    /* ================= FFT INPUT ================= */
    output reg  [2*DATAWIDTH-1:0] s_axis_data_tdata,
    output reg                    s_axis_data_tvalid,
    input  wire                   s_axis_data_tready,
    output reg                    s_axis_data_tlast,

    /* ================= FFT OUTPUT ================= */
    output wire [2*DATAWIDTH-1:0] m_axis_data_tdata,
    output wire                   m_axis_data_tvalid,
    input  wire                   m_axis_data_tready,
    output wire                   m_axis_data_tlast,

    /* ================= OUTPUT SPLIT ================= */
    output wire [DATAWIDTH-1:0] out_real_sine,
    output wire [DATAWIDTH-1:0] out_imag_sine,

    /* ================= EVENTS ================= */
    output wire event_frame_started,
    output wire event_tlast_unexpected,
    output wire event_tlast_missing,
    output wire event_status_channel_halt,
    output wire event_data_in_channel_halt,
    output wire event_data_out_channel_halt
);

    /* ================= FSM ================= */
    reg [1:0] state;
    localparam IDLE   = 2'b00;
    localparam CFG    = 2'b01;
    localparam STREAM = 2'b10;
    localparam DONE   = 2'b11;

    /* ================= CONFIG ================= */
    reg [CONFIGWIDTH-1:0] config_tdata;
    reg                   config_tvalid;

    assign s_axis_config_tdata  = config_tdata;
    assign s_axis_config_tvalid = config_tvalid;

    /* ================= BRAM ================= */
    reg  [ADDRESSWIDTH-1:0] bram_addr;
    wire [31:0]             bram_dout;

    blk_mem_gen_0 u_bram (
        .clka  (aclk),
        .addra (bram_addr),
        .douta (bram_dout)
    );

    /* ================= COUNTER ================= */
    reg [6:0] sample_cnt;

    /* ================= FSM LOGIC ================= */
    always @(posedge aclk) begin
        if (!aresetn) begin
            state               <= IDLE;
            bram_addr           <= 0;
            sample_cnt          <= 0;
            config_tvalid       <= 0;
            s_axis_data_tvalid  <= 0;
            s_axis_data_tlast   <= 0;

            // FFT CONFIG (example: FWD FFT, NFFT=6 → 64)
            config_tdata <= 48'h0000_0000_0000_0106;
        end else begin
            case (state)

            IDLE: begin
                config_tvalid <= 1'b1;
                state <= CFG;
            end

            CFG: begin
                if (s_axis_config_tready) begin
                    config_tvalid <= 0;
                    state <= STREAM;
                end
            end

            STREAM: begin
                if (s_axis_data_tready) begin
                    s_axis_data_tvalid <= 1'b1;
                    s_axis_data_tdata  <= bram_dout; // {imag, real}
                    s_axis_data_tlast  <= (sample_cnt == FFT_LEN-1);

                    bram_addr  <= bram_addr + 1;
                    sample_cnt <= sample_cnt + 1;

                    if (sample_cnt == FFT_LEN-1) begin
                        sample_cnt <= 0;
                        state <= DONE;
                    end
                end
            end

            DONE: begin
                s_axis_data_tvalid <= 0;
                s_axis_data_tlast  <= 0;
                state <= IDLE;
            end

            endcase
        end
    end

    /* ================= FFT IP ================= */
    xfft_0 otfs_fft (
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

        .event_frame_started        (event_frame_started),
        .event_tlast_unexpected     (event_tlast_unexpected),
        .event_tlast_missing        (event_tlast_missing),
        .event_status_channel_halt  (event_status_channel_halt),
        .event_data_in_channel_halt (event_data_in_channel_halt),
        .event_data_out_channel_halt(event_data_out_channel_halt)
    );

    /* ================= OUTPUT SPLIT ================= */
    assign out_real_sine = m_axis_data_tdata[15:0];
    assign out_imag_sine = m_axis_data_tdata[31:16];

endmodule


2. Updated Testbench (coe_waveforms)
BRAM-driven, no sine generation
tb_Top_fft.v

`timescale 1ns / 1ps
module coe_waveforms #(
    parameter DATAWIDTH    = 16,
    parameter ADDRESSWIDTH = 16,
    parameter CONFIGWIDTH  = 48
)();

    /* ================= CLOCK & RESET ================= */
    reg aclk = 0;
    reg aresetn = 0;
    always #5 aclk = ~aclk;

    initial begin
        #50 aresetn = 1;
    end

    /* ================= FFT WIRES ================= */
    wire [CONFIGWIDTH-1:0] s_axis_config_tdata;
    wire                   s_axis_config_tvalid;
    wire                   s_axis_config_tready;

    wire [2*DATAWIDTH-1:0] s_axis_data_tdata;
    wire                   s_axis_data_tvalid;
    wire                   s_axis_data_tready;
    wire                   s_axis_data_tlast;

    wire [2*DATAWIDTH-1:0] m_axis_data_tdata;
    wire                   m_axis_data_tvalid;
    reg                    m_axis_data_tready = 1;
    wire                   m_axis_data_tlast;

    wire [DATAWIDTH-1:0] out_real_sine;
    wire [DATAWIDTH-1:0] out_imag_sine;

    /* ================= DUT ================= */
    Top_fft uut (
        .aclk(aclk),
        .aresetn(aresetn),

        .s_axis_config_tdata (s_axis_config_tdata),
        .s_axis_config_tvalid(s_axis_config_tvalid),
        .s_axis_config_tready(s_axis_config_tready),

        .s_axis_data_tdata (s_axis_data_tdata),
        .s_axis_data_tvalid(s_axis_data_tvalid),
        .s_axis_data_tready(s_axis_data_tready),
        .s_axis_data_tlast (s_axis_data_tlast),

        .m_axis_data_tdata (m_axis_data_tdata),
        .m_axis_data_tvalid(m_axis_data_tvalid),
        .m_axis_data_tready(m_axis_data_tready),
        .m_axis_data_tlast (m_axis_data_tlast),

        .out_real_sine(out_real_sine),
        .out_imag_sine(out_imag_sine),

        .event_frame_started(),
        .event_tlast_unexpected(),
        .event_tlast_missing(),
        .event_status_channel_halt(),
        .event_data_in_channel_halt(),
        .event_data_out_channel_halt()
    );

    /* ================= MONITOR ================= */
    always @(posedge aclk) begin
        if (m_axis_data_tvalid) begin
            $display("FFT OUT: REAL=%d IMAG=%d LAST=%b",
                $signed(out_real_sine),
                $signed(out_imag_sine),
                m_axis_data_tlast);
        end

        if (m_axis_data_tvalid && m_axis_data_tlast) begin
            $display("✔ Full FFT frame received");
            #100 $finish;
        end
    end

endmodule

FINAL MATLAB COE + GOLDEN INPUT GENERATION (MATCHES RTL)

% ============================================================
% PARAMETERS (MATCH RTL + FFT IP)
% ============================================================

N    = 64;     % FFT length (64-point FFT)
GRID = 64;     % Number of FFT frames (OTFS grid size)
IN_BW  = 16;   % Input bit-width (signed, BRAM format)
OUT_BW = 24;   % FFT IP output bit-width (per real/imag)


% ============================================================
% GENERATE 64x64 OTFS INPUT GRID (COMPLEX)
% Each column = one FFT frame
% ============================================================

otfs_grid = zeros(N, GRID);   % 64 rows × 64 columns

for c = 1:GRID
    for r = 1:N
        % Deterministic complex exponential
        % Ideal for FFT verification and OTFS testing
        otfs_grid(r,c) = exp(1j * 2*pi * (r-1) * (c-1) / N);
    end
end


% ============================================================
% COLUMN-WISE SERIALIZATION (OTFS CORRECT ORDER)
% 0–63, 64–127, ..., 4032–4095
% ============================================================

fft_in_serial = reshape(otfs_grid, [], 1);   % 4096×1 complex vector


% ============================================================
% FIXED-POINT SCALING (INPUT SIDE)
% Matches BRAM → FFT IP input
% ============================================================

in_scale = 2^(IN_BW-2);    % Leave 2 bits of headroom

re_in = round(real(fft_in_serial) * in_scale);
im_in = round(imag(fft_in_serial) * in_scale);

% Saturate to signed 16-bit
re_in = max(min(re_in,  2^(IN_BW-1)-1), -2^(IN_BW-1));
im_in = max(min(im_in,  2^(IN_BW-1)-1), -2^(IN_BW-1));


% ============================================================
% PACK REAL & IMAG INTO 32-BIT WORD (BRAM FORMAT)
% [31:16] → Real, [15:0] → Imag
% ============================================================

bram_words = zeros(length(re_in), 1, 'uint32');

for k = 1:length(re_in)
    bram_words(k) = bitshift(uint32(typecast(int16(re_in(k)), 'uint16')), 16) ...
                  + uint32(typecast(int16(im_in(k)), 'uint16'));
end


% ============================================================
% WRITE COE FILE (FOR blk_mem_gen_0)
% ============================================================

fid = fopen('fft_input.coe', 'w');

fprintf(fid, 'memory_initialization_radix=16;\n');
fprintf(fid, 'memory_initialization_vector=\n');

fprintf(fid, '%08X,\n', bram_words(1:end-1));
fprintf(fid, '%08X;\n',  bram_words(end));

fclose(fid);


% ============================================================
% SAVE INPUT DATA FOR MATLAB ↔ RTL GOLDEN COMPARISON
% ============================================================

save('fft_input.mat', 'fft_in_serial');

Assumptions (must match your FFT IP)
Setting
Value
FFT Length
64
Direction
Forward FFT
Scaling
Unscaled
Input width
16-bit
Output width
24-bit
Output format
Fixed-point
(This matches what you’ve been using so far.)
✅ MATLAB FFT OUTPUT COE GENERATOR (FINAL)
🔹 fft_output_golden.m

% ============================================================
% PARAMETERS (MATCH FFT IP + RTL)
% ============================================================

N       = 64;      % FFT length
GRID    = 64;      % Number of FFT frames
IN_BW   = 16;      % Input bit-width
OUT_BW  = 24;      % FFT output bit-width (per real/imag)


% ============================================================
% LOAD INPUT USED FOR RTL (SAME AS BRAM CONTENT)
% ============================================================

load('fft_input.mat');    % loads: fft_in_serial (4096x1 complex)


% ============================================================
% RESHAPE INTO FRAMES (COLUMN-WISE)
% ============================================================

fft_in = reshape(fft_in_serial, N, GRID);


% ============================================================
% PERFORM FFT (MATCH XILINX FFT IP)
% - Column-wise FFT
% - No normalization
% ============================================================

fft_out = zeros(N, GRID);

for c = 1:GRID
    fft_out(:,c) = fft(fft_in(:,c));
end


% ============================================================
% SERIALIZE OUTPUT (FRAME ORDER)
% ============================================================

fft_out_serial = reshape(fft_out, [], 1);   % 4096x1 complex


% ============================================================
% FIXED-POINT SCALING (FFT OUTPUT)
% ------------------------------------------------------------
% Xilinx FFT IP:
% Unscaled → grows by log2(N) = 6 bits
% Input scale = 2^(IN_BW-2)
% ============================================================

growth_bits = log2(N);     % = 6
out_scale = 2^(IN_BW-2) * 2^growth_bits;

re_out = round(real(fft_out_serial) * out_scale);
im_out = round(imag(fft_out_serial) * out_scale);


% ============================================================
% SATURATION TO 24-BIT SIGNED
% ============================================================

MAX_VAL =  2^(OUT_BW-1)-1;
MIN_VAL = -2^(OUT_BW-1);

re_out = max(min(re_out, MAX_VAL), MIN_VAL);
im_out = max(min(im_out, MAX_VAL), MIN_VAL);


% ============================================================
% PACK INTO 48-BIT WORD
% [47:24] → REAL
% [23:0]  → IMAG
% ============================================================

fft_out_words = zeros(length(re_out), 1, 'uint64');

for k = 1:length(re_out)
    real_u = uint32(typecast(int32(re_out(k)), 'uint32'));
    imag_u = uint32(typecast(int32(im_out(k)), 'uint32'));

    fft_out_words(k) = bitshift(uint64(real_u), 24) ...
                     + uint64(bitand(imag_u, hex2dec('FFFFFF')));
end


% ============================================================
% WRITE GOLDEN FFT OUTPUT COE
% ============================================================

fid = fopen('fft_output_golden.coe', 'w');

fprintf(fid, 'memory_initialization_radix=16;\n');
fprintf(fid, 'memory_initialization_vector=\n');

fprintf(fid, '%012X,\n', fft_out_words(1:end-1));
fprintf(fid, '%012X;\n',  fft_out_words(end));

fclose(fid);

disp('✔ fft_output_golden.coe generated successfully');


% ============================================================
% SAVE FOR MATLAB ↔ RTL NUMERICAL COMPARISON
% ============================================================

save('fft_output_golden.mat', 're_out', 'im_out');

How to use this in practice
🔹 RTL side
Capture FFT output (m_axis_data_tdata)
Store or dump it to a text/COE file
🔹 MATLAB side
Load fft_output_golden.mat
Compare:
error_real = rtl_real - re_out;
error_imag = rtl_imag - im_out;
max(abs(error_real))
max(abs(error_imag))
Expected result:
Error ≤ ±1 LSB
Pro tip (very important)
If you later enable scaling schedule in FFT IP:
out_scale must change
Otherwise MATLAB and RTL will not match
