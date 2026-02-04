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
    parameter CONFIGWIDTH   = 16,
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
    parameter CONFIGWIDTH  = 16
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

`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
module Top_fft #(
    parameter DATAWIDTH    = 16,
    parameter ADDRESSWIDTH = 6,    // 64 deep BRAM
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
    output reg  [31:0] s_axis_data_tdata,
    output reg         s_axis_data_tvalid,
    input  wire        s_axis_data_tready,
    output reg         s_axis_data_tlast,

    /* ================= FFT OUTPUT ================= */
    output wire [47:0] m_axis_data_tdata,
    output wire        m_axis_data_tvalid,
    input  wire        m_axis_data_tready,
    output wire        m_axis_data_tlast,

    /* ================= OUTPUT SPLIT ================= */
    output wire signed [22:0] out_real_sine,
    output wire signed [22:0] out_imag_sine,

    /* ================= EVENTS ================= */
    output wire event_frame_started,
    output wire event_tlast_unexpected,
    output wire event_tlast_missing,
    output wire event_status_channel_halt,
    output wire event_data_in_channel_halt,
    output wire event_data_out_channel_halt
);

    /* ================= FSM ================= */
    localparam IDLE        = 3'd0;
    localparam SEND_CFG    = 3'd1;
    localparam CFG_WAIT    = 3'd2;
    localparam SEND_DATA   = 3'd3;
    localparam WAIT_OUTPUT = 3'd4;

    reg [2:0] state;

    /* ================= CONFIG ================= */
    reg config_valid;
    assign s_axis_config_tdata  = 16'h0106;   // FFT, NFFT=6 (64)
    assign s_axis_config_tvalid = config_valid;

    /* ================= BRAM ================= */
    reg  [ADDRESSWIDTH-1:0] bram_addr;
    wire [31:0] bram_dout;
    reg  [31:0] bram_dout_r;   // pipeline register

    blk_mem_gen_0 u_bram (
        .clka  (aclk),
        .addra (bram_addr),
        .douta (bram_dout)
    );

    always @(posedge aclk)
        bram_dout_r <= bram_dout;

    /* ================= COUNTER ================= */
    reg [6:0] sample_cnt;

    /* ================= FSM ================= */
    always @(posedge aclk) begin
        if (!aresetn) begin
            state              <= IDLE;
            config_valid       <= 1'b0;
            s_axis_data_tvalid <= 1'b0;
            s_axis_data_tlast  <= 1'b0;
            bram_addr          <= 0;
            sample_cnt         <= 0;
        end else begin
            case (state)

                IDLE: begin
                    config_valid <= 1'b1;
                    state        <= SEND_CFG;
                end

                SEND_CFG: begin
                    if (s_axis_config_tready) begin
                        config_valid <= 1'b0;
                        state <= CFG_WAIT;
                    end
                end

                CFG_WAIT: begin
                    bram_addr  <= 0;
                    sample_cnt <= 0;
                    state <= SEND_DATA;
                end

                SEND_DATA: begin
                    if (s_axis_data_tready) begin
                        s_axis_data_tvalid <= 1'b1;
                        s_axis_data_tdata  <= bram_dout_r;
                        s_axis_data_tlast  <= (sample_cnt == FFT_LEN-1);

                        bram_addr  <= bram_addr + 1;
                        sample_cnt <= sample_cnt + 1;

                        if (sample_cnt == FFT_LEN-1) begin
                            s_axis_data_tvalid <= 1'b0;
                            s_axis_data_tlast  <= 1'b0;
                            state <= WAIT_OUTPUT;
                        end
                    end
                end

                WAIT_OUTPUT: begin
                    if (m_axis_data_tvalid && m_axis_data_tlast)
                        state <= IDLE;
                end

                default: state <= IDLE;
            endcase
        end
    end

    /* ================= FFT IP ================= */
    xfft_0 fft_i (
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

        .event_frame_started         (event_frame_started),
        .event_tlast_unexpected      (event_tlast_unexpected),
        .event_tlast_missing         (event_tlast_missing),
        .event_status_channel_halt   (event_status_channel_halt),
        .event_data_in_channel_halt  (event_data_in_channel_halt),
        .event_data_out_channel_halt (event_data_out_channel_halt)
    );

    /* ================= OUTPUT SPLIT ================= */
    assign out_real_sine = m_axis_data_tdata[22:0];
    assign out_imag_sine = m_axis_data_tdata[46:24];

endmodule


fft_output_data

% ============================================================
% FFT OUTPUT COE GENERATION (64-point, 48-bit)
% ============================================================

N = 64;

IN_BW  = 16;   % input width
OUT_BW = 23;   % FFT IP output width per real/imag

% ------------------------------------------------------------
% Load input used for RTL (same as input COE)
% ------------------------------------------------------------
load('fft_input.mat');   % fft_in_serial (complex, length 64)

% Undo input scaling
in_scale = 2^(IN_BW-2);
x = fft_in_serial / in_scale;

% ------------------------------------------------------------
% MATLAB FFT (golden reference)
% ------------------------------------------------------------
y = fft(x, N);

% ------------------------------------------------------------
% FFT IP GROWTH (UNSCALED MODE)
% Vivado FFT grows by log2(N) bits
% ------------------------------------------------------------
growth = log2(N);        % = 6 for N=64
y = y * (2^growth);

% ------------------------------------------------------------
% FIXED-POINT CONVERSION (TRUNCATION)
% ------------------------------------------------------------
y_re = floor(real(y));
y_im = floor(imag(y));

% Saturate to signed 23-bit
MAX =  2^(OUT_BW-1)-1;
MIN = -2^(OUT_BW-1);

y_re = max(min(y_re, MAX), MIN);
y_im = max(min(y_im, MAX), MIN);

% ------------------------------------------------------------
% PACK INTO 48-BIT WORD
% [46:24] = IMAG, [22:0] = REAL
% ------------------------------------------------------------
fft_out_words = zeros(N,1,'uint64');

for k = 1:N
    fft_out_words(k) = ...
        bitshift(uint64(typecast(int32(y_im(k)), 'uint32')), 24) + ...
        uint64(typecast(int32(y_re(k)), 'uint32'));
end

% ------------------------------------------------------------
% WRITE COE FILE
% ------------------------------------------------------------
fid = fopen('fft_output_64.coe','w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k = 1:N-1
    fprintf(fid,'%012X,\n', fft_out_words(k));
end
fprintf(fid,'%012X;\n', fft_out_words(N));
fclose(fid);

disp('✅ fft_output_64.coe generated successfully');


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
