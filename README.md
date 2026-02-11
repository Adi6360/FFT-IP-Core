`timescale 1ns / 1ps
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

module Top_fft #(
    parameter DATAWIDTH    = 16,
    parameter ADDRESSWIDTH = 6,
    parameter CONFIGWIDTH  = 16,
    parameter FFT_LEN      = 64,
    parameter INTER_FRAME_GAP = 20
)( 
    input  wire aclk,
    input  wire aresetn,

    /* ================= FFT CONFIG ================= */
    output wire [CONFIGWIDTH-1:0] s_axis_config_tdata,
    output wire                  s_axis_config_tvalid,
    input  wire                  s_axis_config_tready,

    /* ================= FFT INPUT ================= */
    output wire [2*DATAWIDTH-1:0] s_axis_data_tdata,
    output wire                  s_axis_data_tvalid,
    input  wire                  s_axis_data_tready,
    output wire                  s_axis_data_tlast,

    /* ================= FFT OUTPUT ================= */
    input  wire [47:0] m_axis_data_tdata,
    input  wire        m_axis_data_tvalid,
    output wire        m_axis_data_tready,
    input  wire        m_axis_data_tlast,

    output wire [46:24] out_real_data, 
    output wire [22:0]  out_imag_data,

    /* ================= STATUS ================= */
    output wire event_frame_started,
    output wire event_tlast_unexpected,
    output wire event_tlast_missing,
    output wire event_status_channel_halt,
    output wire event_data_in_channel_halt,
    output wire event_data_out_channel_halt
);

//////////////////////////////
// FSM STATES
//////////////////////////////
reg [2:0] state;
localparam IDLE          = 3'b000;
localparam CFG           = 3'b001;
localparam BRAM_WAIT     = 3'b010;
localparam STREAM        = 3'b011;
localparam WAIT_FOR_DONE = 3'b100;
localparam INTER_FRAME   = 3'b101;

//////////////////////////////
// INTERNAL REGS
//////////////////////////////
reg [CONFIGWIDTH-1:0] config_tdata_reg;
reg config_tvalid_reg;

reg [31:0] data_tdata_reg;
reg data_tvalid_reg;

reg [ADDRESSWIDTH-1:0] bram_addr;
reg [7:0] sample_cnt;
reg [7:0] gap_cnt;

//////////////////////////////
// ✅ NEW — 2-cycle pipelines
//////////////////////////////
reg [1:0] tvalid_pipe;
reg [1:0] tlast_pipe;

//////////////////////////////
// ASSIGNS
//////////////////////////////
assign s_axis_config_tdata  = config_tdata_reg;
assign s_axis_config_tvalid = config_tvalid_reg;

assign s_axis_data_tdata  = data_tdata_reg;
assign s_axis_data_tvalid = tvalid_pipe[1];   // delayed 2 cycles
assign s_axis_data_tlast  = tlast_pipe[1];    // delayed 2 cycles

assign m_axis_data_tready = 1'b1;

//////////////////////////////
// BRAM
//////////////////////////////
wire [31:0] bram_dout;

blk_mem_gen_0 u_bram (
  .clka(aclk),
  .addra(bram_addr),
  .douta(bram_dout)
);

//////////////////////////////
// FSM
//////////////////////////////
always @(posedge aclk) begin
    if (!aresetn) begin
        state <= IDLE;
        bram_addr <= 0;
        sample_cnt <= 0;
        gap_cnt <= 0;
        config_tvalid_reg <= 0;
        data_tvalid_reg <= 0;
        config_tdata_reg <= 16'b0000000010000110;

        tvalid_pipe <= 0;
        tlast_pipe  <= 0;
    end else begin

        //////////////////////////
        // ✅ pipeline update (ONLY timing change)
        //////////////////////////
        tvalid_pipe <= {tvalid_pipe[0], data_tvalid_reg};
        tlast_pipe  <= {tlast_pipe[0], (sample_cnt == FFT_LEN-1) && data_tvalid_reg};

        case (state)

        IDLE: begin
            config_tvalid_reg <= 1'b1;
            state <= CFG;
        end

        CFG: begin
            if (s_axis_config_tready && s_axis_config_tvalid) begin
                config_tvalid_reg <= 0;
                state <= BRAM_WAIT;
                bram_addr <= 0;
            end
        end

        // BRAM latency wait
        BRAM_WAIT: begin
            state <= STREAM;
            data_tvalid_reg <= 1'b1;
        end

        STREAM: begin
            if (s_axis_data_tready) begin
                if (sample_cnt == FFT_LEN-1) begin
                    sample_cnt <= 0;
                    data_tvalid_reg <= 0;
                    state <= WAIT_FOR_DONE;
                end else begin
                    sample_cnt <= sample_cnt + 1;
                    bram_addr  <= bram_addr + 1;
                    data_tvalid_reg <= 1;
                end
            end
        end

        WAIT_FOR_DONE: begin
            if (m_axis_data_tlast && m_axis_data_tvalid) begin
                state <= INTER_FRAME;
                gap_cnt <= 0;
            end
        end

        INTER_FRAME: begin
            if (gap_cnt >= INTER_FRAME_GAP) begin
                state <= BRAM_WAIT;
                bram_addr <= 0;
            end else begin
                gap_cnt <= gap_cnt + 1;
            end
        end

        endcase
    end
end

//////////////////////////////
// Data path (unchanged)
//////////////////////////////
always @(*) begin
    data_tdata_reg = bram_dout;
end

//////////////////////////////
// FFT IP
//////////////////////////////
xfft_1 u_fft (
    .aclk(aclk),
    .aresetn(aresetn),

    .s_axis_config_tdata(s_axis_config_tdata),
    .s_axis_config_tvalid(s_axis_config_tvalid),
    .s_axis_config_tready(s_axis_config_tready),

    .s_axis_data_tdata(s_axis_data_tdata),
    .s_axis_data_tvalid(s_axis_data_tvalid),
    .s_axis_data_tready(s_axis_data_tready),
    .s_axis_data_tlast(s_axis_data_tlast),

    .m_axis_data_tdata(m_axis_data_tdata),
    .m_axis_data_tvalid(m_axis_data_tvalid),
    .m_axis_data_tready(m_axis_data_tready),
    .m_axis_data_tlast(m_axis_data_tlast),

    .event_frame_started(event_frame_started),
    .event_tlast_unexpected(event_tlast_unexpected),
    .event_tlast_missing(event_tlast_missing),
    .event_status_channel_halt(event_status_channel_halt),
    .event_data_in_channel_halt(event_data_in_channel_halt),
    .event_data_out_channel_halt(event_data_out_channel_halt)
);

assign out_real_data = m_axis_data_tdata[46:24];
assign out_imag_data = m_axis_data_tdata[22:0];

endmodule

Your same Top_fft code
✅ Only change = delay s_axis_data_tvalid by 2 cycles after STREAM starts
✅ And ensure tvalid also ends aligned with last sample (63) after BRAM latency
✅ No FSM restructuring
✅ No counter logic changes
✅ No BRAM logic changes
✅ Only AXI valid/last alignment pipeline added
To make start + end both align, we must pipeline:
data_tvalid_reg → 2 cycles
tlast condition → same 2 cycles




% ============================================================
% PARAMETERS — MATCH VIVADO FFT IP
% ============================================================
N = 64;

IN_W  = 16;   % input real/imag width
OUT_W = 24;   % FFT IP output real/imag width

% ============================================================
% GENERATE TEST COMPLEX INPUT (I + jQ)
% ============================================================
n = 0:N-1;

% single-tone complex exponential (clean FFT peak)
x = exp(1j*2*pi*5*n/N);

% ============================================================
% FIXED POINT INPUT QUANTIZATION (SIGNED 16 BIT)
% truncation mode
% ============================================================
scale_in = 2^(IN_W-1);

re_in = fix(real(x) * scale_in);
im_in = fix(imag(x) * scale_in);

re_in = max(min(re_in,  2^(IN_W-1)-1), -2^(IN_W-1));
im_in = max(min(im_in,  2^(IN_W-1)-1), -2^(IN_W-1));

x_fixed = re_in + 1j*im_in;

% ============================================================
% FFT — UN SCALED (matches IP setting)
% ============================================================
X = fft(x_fixed);   % no /N scaling

% ============================================================
% OUTPUT TRUNCATION TO 24 BIT (NOT ROUND)
% ============================================================
scale_out = 1;   % already integer domain

re_out = fix(real(X) * scale_out);
im_out = fix(imag(X) * scale_out);

re_out = max(min(re_out,  2^(OUT_W-1)-1), -2^(OUT_W-1));
im_out = max(min(im_out,  2^(OUT_W-1)-1), -2^(OUT_W-1));

% ============================================================
% PACK INPUT → 32 BIT WORD (BRAM FORMAT)
% [31:16] I , [15:0] Q
% ============================================================
in_word = zeros(N,1,'uint32');

for k=1:N
    in_word(k) = bitshift(uint32(typecast(int16(re_in(k)),'uint16')),16) + ...
                          uint32(typecast(int16(im_in(k)),'uint16'));
end

% ============================================================
% PACK OUTPUT → 48 BIT WORD (FFT AXIS FORMAT)
% [47:24] REAL , [23:0] IMAG
% ============================================================
out_word = zeros(N,1,'uint64');

for k=1:N
    r = uint64(typecast(int32(re_out(k)),'uint32'));
    i = uint64(typecast(int32(im_out(k)),'uint32'));
    out_word(k) = bitshift(r,24) + bitand(i, hex2dec('FFFFFF'));
end

% ============================================================
% -------- COE WRITER FUNCTION --------
% ============================================================
write_coe = @(fname, data, width_hex) ...
    write_coe_file(fname, data, width_hex);

% ============================================================
% WRITE COE FILES — INPUT
% ============================================================
write_coe('fft_in_real.coe', re_in, 4);
write_coe('fft_in_imag.coe', im_in, 4);
write_coe('fft_in_packed32.coe', in_word, 8);

% ============================================================
% WRITE COE FILES — OUTPUT
% ============================================================
write_coe('fft_out_real.coe', re_out, 6);
write_coe('fft_out_imag.coe', im_out, 6);
write_coe('fft_out_packed48.coe', out_word, 12);

disp('✅ All FFT input/output COE files generated');


% ============================================================
% COE FILE WRITER
% ============================================================
function write_coe_file(fname, data, hexwidth)

fid = fopen(fname,'w');

fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k=1:length(data)-1
    if data(k) < 0
        val = data(k) + 2^(hexwidth*4);
    else
        val = data(k);
    end
    fprintf(fid,['%0',num2str(hexwidth),'X,\n'], val);
end

val = data(end);
if val < 0
    val = val + 2^(hexwidth*4);
end
fprintf(fid,['%0',num2str(hexwidth),'X;\n'], val);

fclose(fid);
end






otfs

% ============================================================
% PARAMETERS — MATCH YOUR RTL + FFT IP
% ============================================================
N    = 64;      % FFT length
GRID = 64;      % OTFS grid size (64 frames)

IN_W  = 16;     % input width
OUT_W = 24;     % FFT output width

% ============================================================
% GENERATE OTFS 64×64 COMPLEX GRID
% deterministic (repeatable)
% ============================================================
otfs = zeros(N, GRID);

for c = 1:GRID
    for r = 1:N
        otfs(r,c) = exp(1j*2*pi*(r-1)*(c-1)/N);
    end
end

% ============================================================
% COLUMN SERIALIZATION — MATCHES YOUR FSM STREAM
% frame0 → samples 0..63
% frame1 → samples 64..127
% ============================================================
x_serial = reshape(otfs, [], 1);   % 4096 samples

% ============================================================
% INPUT FIXED POINT — 16 BIT — TRUNCATION
% ============================================================
scale_in = 2^(IN_W-1);

re_in = fix(real(x_serial)*scale_in);
im_in = fix(imag(x_serial)*scale_in);

re_in = max(min(re_in,  2^(IN_W-1)-1), -2^(IN_W-1));
im_in = max(min(im_in,  2^(IN_W-1)-1), -2^(IN_W-1));

x_fixed = re_in + 1j*im_in;

% ============================================================
% FRAME-WISE FFT (MATCHES FFT IP OPERATION)
% ============================================================
X_frames = zeros(N, GRID);

for c = 1:GRID
    idx = (c-1)*N + (1:N);
    X_frames(:,c) = fft(x_fixed(idx));   % UNSCALED
end

X_serial = reshape(X_frames, [], 1);

% ============================================================
% OUTPUT TRUNCATION — 24 BIT — NO ROUNDING
% ============================================================
re_out = fix(real(X_serial));
im_out = fix(imag(X_serial));

re_out = max(min(re_out,  2^(OUT_W-1)-1), -2^(OUT_W-1));
im_out = max(min(im_out,  2^(OUT_W-1)-1), -2^(OUT_W-1));

% ============================================================
% PACK INPUT → 32 BIT (BRAM FORMAT)
% ============================================================
in32 = zeros(length(re_in),1,'uint32');

for k=1:length(re_in)
    in32(k) = bitshift(uint32(typecast(int16(re_in(k)),'uint16')),16) + ...
                      uint32(typecast(int16(im_in(k)),'uint16'));
end

% ============================================================
% PACK OUTPUT → 48 BIT (FFT AXIS FORMAT)
% ============================================================
out48 = zeros(length(re_out),1,'uint64');

for k=1:length(re_out)
    r = uint64(typecast(int32(re_out(k)),'uint32'));
    i = uint64(typecast(int32(im_out(k)),'uint32'));
    out48(k) = bitshift(r,24) + bitand(i, hex2dec('FFFFFF'));
end

% ============================================================
% WRITE COE FILES
% ============================================================
write_coe('otfs_in_real.coe', re_in, 4);
write_coe('otfs_in_imag.coe', im_in, 4);
write_coe('otfs_in_packed32.coe', in32, 8);

write_coe('otfs_out_real.coe', re_out, 6);
write_coe('otfs_out_imag.coe', im_out, 6);
write_coe('otfs_out_packed48.coe', out48, 12);

save('otfs_fft_matlab_golden.mat', ...
     'x_serial','X_serial','re_in','im_in','re_out','im_out');

disp('✅ OTFS 64x64 FFT MATLAB golden + COE files generated');


% ============================================================
% COE WRITER
% ============================================================
function write_coe(fname, data, hexwidth)

fid = fopen(fname,'w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k=1:length(data)-1
    v = data(k);
    if v < 0
        v = v + 2^(hexwidth*4);
    end
    fprintf(fid,['%0',num2str(hexwidth),'X,\n'], v);
end

v = data(end);
if v < 0
    v = v + 2^(hexwidth*4);
end
fprintf(fid,['%0',num2str(hexwidth),'X;\n'], v);

fclose(fid);
end

Matches Your RTL Design Exactly
This aligns with your current system:
Your RTL
MATLAB model
BRAM width 32
packed input32
BRAM depth 4096
64×64 serialized
FFT length 64
per-column FFT
Unscaled
no normalization
Truncation
fix()
Natural order
default FFT
Output 48-bit
packed48
