`timescale 1ns / 1ps 

module Top_fft #( 
        parameter DATAWIDTH = 16,
        parameter ADDRESSWIDTH = 6, 
        parameter CONFIGWIDTH = 16, 
        parameter FFT_LEN = 64, 
        parameter INTER_FRAME_GAP = 20 )
( 
    input wire aclk, 
    input wire aresetn,

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

output wire [22:0] out_real_data, 
output wire [46:24]  out_imag_data,

/* ================= STATUS ================= */

output wire event_frame_started,
output wire event_tlast_unexpected,
output wire event_tlast_missing,
output wire event_status_channel_halt,
output wire event_data_in_channel_halt,
output wire event_data_out_channel_halt
);

////////////////////////////// // FSM STATES //////////////////////////////
 
reg [2:0] state; 
localparam IDLE = 3'b000;
localparam CFG = 3'b001; 
localparam BRAM_WAIT = 3'b010; 
localparam STREAM = 3'b011; 
localparam WAIT_FOR_DONE = 3'b100; 
localparam INTER_FRAME = 3'b101;

////////////////////////////// // INTERNAL REGS ////////////////////////////// 
    reg [CONFIGWIDTH-1:0] config_tdata_reg; 
    reg config_tvalid_reg;

    reg [31:0] data_tdata_reg; 
    reg data_tvalid_reg;

    reg [ADDRESSWIDTH-1:0] bram_addr;
    reg [7:0] sample_cnt; 
    reg [7:0] gap_cnt;

////////////////////////////// // ✅ NEW - 2-cycle pipelines ////////////////////////////// 
    reg [1:0] tvalid_pipe; 
    reg [1:0] tlast_pipe;

////////////////////////////// // ASSIGNS ////////////////////////////// 
assign s_axis_config_tdata = config_tdata_reg; 
assign s_axis_config_tvalid = config_tvalid_reg;

assign s_axis_data_tdata = data_tdata_reg;
assign s_axis_data_tvalid = tvalid_pipe[1]; // delayed 2 cycles 
assign s_axis_data_tlast = tlast_pipe[1]; // delayed 2 cycles

assign m_axis_data_tready = 1'b1;

//    BRAM    // 
wire [31:0] bram_dout;

blk_mem_gen_0 u_bram ( 
.clka(aclk), 
.addra(bram_addr), 
.douta(bram_dout) 
);

////////////////////////////// // FSM ////////////////////////////// 
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
             end else begin
               data_tvalid_reg <= 0;   
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

////////////////////////////// // Data path (unchanged) ////////////////////////////// 
always @(*) begin 
data_tdata_reg = bram_dout;
 end

////////////////////////////// // FFT IP ////////////////////////////// 
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

assign out_real_data = m_axis_data_tdata[22:0]; 
assign out_imag_data = m_axis_data_tdata[46:24];

endmodule

`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
module fft_tb;

parameter DATAWIDTH    = 16;
parameter ADDRESSWIDTH = 6;
parameter CONFIGWIDTH  = 16;
parameter FFT_LEN      = 64;


// CLOCK / RESET////

reg aclk = 0;
reg aresetn = 0;

always #5 aclk = ~aclk;   // 100 MHz

initial begin
    aresetn = 0;
    #100;
    aresetn = 1;
end

//// AXI SIGNALS    ///////

wire [CONFIGWIDTH-1:0] s_axis_config_tdata;
wire                   s_axis_config_tvalid;
wire                   s_axis_config_tready;

wire [7:0]             sample_cnt;
wire [7:0]             bram_addr; 

wire [2*DATAWIDTH-1:0] s_axis_data_tdata;
wire                   s_axis_data_tvalid;
wire                   s_axis_data_tready;
wire                   s_axis_data_tlast;


wire [47:0]            m_axis_data_tdata;
wire                   m_axis_data_tvalid;
wire                   m_axis_data_tready;
wire                   m_axis_data_tlast;

wire [22:0] out_real_data;
wire [46:24]  out_imag_data;

//   EVENTS   //

wire event_frame_started;
wire event_tlast_unexpected;
wire event_tlast_missing;
wire event_status_channel_halt;
wire event_data_in_channel_halt;
wire event_data_out_channel_halt;



////////////////////////////////////////////////////////////
Top_fft uut (
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
    
    .out_real_data(out_real_data),
    .out_imag_data(out_imag_data),

    .m_axis_data_tvalid  (m_axis_data_tvalid),
    .m_axis_data_tready  (m_axis_data_tready),
    .m_axis_data_tlast   (m_axis_data_tlast),
    .event_frame_started(event_frame_started),
    .event_tlast_unexpected(event_tlast_unexpected),
    .event_tlast_missing(event_tlast_missing),
    .event_status_channel_halt(event_status_channel_halt),
    .event_data_in_channel_halt(event_data_in_channel_halt),
    .event_data_out_channel_halt(event_data_out_channel_halt)
);


// OUTPUT MONITOR //
always @(posedge aclk) begin
    if (m_axis_data_tvalid) begin
        
        $display("TIME:%t | FFT OUT: REAL=%d IMAG=%d TLAST=%b",$time, $signed(out_real_data),$signed(out_imag_data),m_axis_data_tlast);
    end

    if (uut.event_tlast_unexpected) $display("ERROR: TLAST Unexpected");

end

endmodule



// MATLAB //

% ============================================================
% XILINX FFT IP TRUE MATCH MODEL (TWIDDLE QUANTIZED)
% Matches: pipelined streaming, fixed, truncation, unscaled
% ============================================================

clear; clc;

N = 64;
IN_W  = 16;
TW_W  = 16;
OUT_W = 23;

% ============================================================
% INPUT SIGNAL (same as your sim)
% ============================================================
n = 0:N-1;
x = exp(1j*2*pi*5*n/N);

SCALE = 2^(IN_W-1);
xr = sat(fix(real(x)*SCALE), IN_W);
xi = sat(fix(imag(x)*SCALE), IN_W);
x_fixed = xr + 1j*xi;

% ============================================================
% BUILD QUANTIZED TWIDDLE MATRIX (16-bit like IP)
% ============================================================
TW_SCALE = 2^(TW_W-1);

W = zeros(N,N);

for k = 0:N-1
    for m = 0:N-1
        w = exp(-1j*2*pi*k*m/N);
        wr = sat(fix(real(w)*TW_SCALE), TW_W);
        wi = sat(fix(imag(w)*TW_SCALE), TW_W);
        W(k+1,m+1) = (wr + 1j*wi) / TW_SCALE;   % back to scaled fixed
    end
end

% ============================================================
% FFT USING QUANTIZED TWIDDLES
% ============================================================
X = W * x_fixed(:);

% ============================================================
% FINAL TRUNCATION ONLY (like IP output stage)
% ============================================================
re_out = sat(fix(real(X)), OUT_W);
im_out = sat(fix(imag(X)), OUT_W);

% ============================================================
% PACK INPUT → 32 BIT (SWAPPED)
% [31:16]=IMAG  [15:0]=REAL
% ============================================================
in_word = zeros(N,1,'uint32');

for k=1:N
    r = typecast(int16(xr(k)),'uint16');
    i = typecast(int16(xi(k)),'uint16');
    in_word(k) = bitshift(uint32(i),16) + uint32(r);
end

% ============================================================
% PACK OUTPUT → 48 BIT (SWAPPED)
% [47:24]=IMAG  [23:0]=REAL
% ============================================================
MASK24 = uint64(hex2dec('FFFFFF'));
out_word = zeros(N,1,'uint64');

for k=1:N
    r = bitand(uint64(typecast(int32(re_out(k)),'uint32')), MASK24);
    i = bitand(uint64(typecast(int32(im_out(k)),'uint32')), MASK24);
    out_word(k) = bitshift(i,24) + r;
end

% ============================================================
% WRITE ALL 6 COE FILES
% ============================================================
write_coe('fft_in_real.coe', xr, 4);
write_coe('fft_in_imag.coe', xi, 4);
write_coe('fft_in_packed32.coe', in_word, 8);

write_coe('fft_out_real.coe', re_out, 6);
write_coe('fft_out_imag.coe', im_out, 6);
write_coe('fft_out_packed48.coe', out_word, 12);

disp('✅ COE files generated — twiddle-quantized FFT model');

% ============================================================
% HELPERS
% ============================================================
function y = sat(x,W)
L = 2^(W-1);
y = min(max(x,-L),L-1);
end

function write_coe(fname,data,hexw)
fid=fopen(fname,'w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k=1:length(data)-1
    v = data(k);
    if v < 0
        v = v + 2^(hexw*4);
    end
    fprintf(fid,['%0',num2str(hexw),'X,\n'], v);
end

v = data(end);
if v < 0
    v = v + 2^(hexw*4);
end
fprintf(fid,['%0',num2str(hexw),'X;\n'], v);
fclose(fid);
end



  // Matlab fft output dump to coe //
  % ============================================================
% FFT IP HEX DUMP → COE GENERATOR (BIT-EXACT)
% ============================================================

clear; clc;

INFILE  = 'fft_ip_dump.txt';
OUT_COE = 'fft_out_packed48.coe';

% ============================================================
% READ HEX LINES
% ============================================================
raw = readlines(INFILE);
raw = strtrim(raw);

% remove commas / semicolons
raw = erase(raw, ',');
raw = erase(raw, ';');

% remove empty lines
raw = raw(raw ~= "");

N = length(raw);
data48 = zeros(N,1,'uint64');

for k = 1:N
    data48(k) = uint64(hex2dec(raw(k)));
end

% ============================================================
% WRITE PACKED 48-BIT COE (EXACT COPY)
% ============================================================
fid = fopen(OUT_COE,'w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');

for k = 1:N-1
    fprintf(fid,'%012X,\n', data48(k));
end
fprintf(fid,'%012X;\n', data48(N));
fclose(fid);

disp('✅ Packed 48-bit COE written (bit-exact)');

% ============================================================
% ALSO GENERATE SPLIT REAL / IMAG COE FILES
% Packing assumed:
% [47:24] = IMAG
% [23:0]  = REAL
% ============================================================

real24 = bitand(data48, hex2dec('FFFFFF'));
imag24 = bitshift(data48, -24);

write24coe('fft_out_real.coe', real24);
write24coe('fft_out_imag.coe', imag24);

disp('✅ Split REAL/IMAG COE files written');

% ============================================================
% helper
% ============================================================
function write24coe(fname,vec)
fid=fopen(fname,'w');
fprintf(fid,'memory_initialization_radix=16;\n');
fprintf(fid,'memory_initialization_vector=\n');
for i=1:length(vec)-1
    fprintf(fid,'%06X,\n', vec(i));
end
fprintf(fid,'%06X;\n', vec(end));
fclose(fid);
end

