`timescale 1ns / 1ps ////////////////////////////////////////////////////////////////////////////////// // Cleaned FINAL Top_fft — only timing & BRAM pipeline fixes applied ////////////////////////////////////////////////////////////////////////////////// module Top_fft #( parameter DATAWIDTH    = 16, parameter ADDRESSWIDTH = 6,   // 64 depth parameter CONFIGWIDTH  = 16, parameter FFT_LEN      = 64, parameter INTER_FRAME_GAP = 8 )( input  wire aclk, input  wire aresetn,

/* ================= FFT CONFIG ================= */
output wire [CONFIGWIDTH-1:0] s_axis_config_tdata,
output wire s_axis_config_tvalid,
input  wire s_axis_config_tready,

/* ================= FFT INPUT ================= */
output wire [2*DATAWIDTH-1:0] s_axis_data_tdata,
output wire s_axis_data_tvalid,
input  wire s_axis_data_tready,
output wire s_axis_data_tlast,

/* ================= FFT OUTPUT ================= */
input  wire [31:0] m_axis_data_tdata,
input  wire m_axis_data_tvalid,
output wire m_axis_data_tready,
input  wire m_axis_data_tlast,

/* ================= OUTPUT SPLIT ================= */
output wire [23:0] out_real_data,
output wire [23:0] out_imag_data,

/* ================= EVENTS ================= */
output wire event_frame_started,
output wire event_tlast_unexpected,
output wire event_tlast_missing,
output wire event_status_channel_halt,
output wire event_data_in_channel_halt,
output wire event_data_out_channel_halt

);

//////////////////////////////////////////////////////////// // FSM STATES //////////////////////////////////////////////////////////// localparam IDLE        = 3'd0; localparam CFG         = 3'd1; localparam BRAM_WAIT   = 3'd2; localparam BRAM_PRIME  = 3'd3; localparam STREAM      = 3'd4; localparam WAIT_DONE   = 3'd5; localparam INTER_FRAME = 3'd6;

reg [2:0] state;

//////////////////////////////////////////////////////////// // REGISTERS //////////////////////////////////////////////////////////// reg [CONFIGWIDTH-1:0]   config_tdata_reg; reg                     config_tvalid_reg;

reg [2*DATAWIDTH-1:0]   data_tdata_reg; reg                     data_tvalid_reg;

reg [ADDRESSWIDTH-1:0]  bram_addr; reg [7:0]               sample_cnt; reg [7:0]               gap_cnt;

//////////////////////////////////////////////////////////// // HANDSHAKE //////////////////////////////////////////////////////////// wire stream_fire; assign stream_fire = s_axis_data_tvalid && s_axis_data_tready;

//////////////////////////////////////////////////////////// // ASSIGNS //////////////////////////////////////////////////////////// assign s_axis_config_tdata  = config_tdata_reg; assign s_axis_config_tvalid = config_tvalid_reg;

assign s_axis_data_tdata  = data_tdata_reg; assign s_axis_data_tvalid = data_tvalid_reg; assign s_axis_data_tlast  = (sample_cnt == FFT_LEN-1) && data_tvalid_reg;

assign m_axis_data_tready = 1'b1;

//////////////////////////////////////////////////////////// // BRAM //////////////////////////////////////////////////////////// wire [31:0] bram_dout; reg  [31:0] bram_dout_r;

blk_mem_gen_0 u_bram ( .clka  (aclk), .addra (bram_addr), .douta (bram_dout) );

// BRAM output pipeline register (required) always @(posedge aclk) bram_dout_r <= bram_dout;

//////////////////////////////////////////////////////////// // FSM //////////////////////////////////////////////////////////// always @(posedge aclk) begin if (!aresetn) begin state <= IDLE; bram_addr <= 0; sample_cnt <= 0; gap_cnt <= 0;

config_tvalid_reg <= 0;
    data_tvalid_reg   <= 0;

    config_tdata_reg <= 16'b0000_0000_1000_0110; // forward FFT
end
else begin
    case (state)

    IDLE: begin
        config_tvalid_reg <= 1;
        state <= CFG;
    end

    CFG: begin
        if (s_axis_config_tready && config_tvalid_reg) begin
            config_tvalid_reg <= 0;
            bram_addr <= 0;
            state <= BRAM_WAIT;
        end
    end

    // allow BRAM one cycle latency
    BRAM_WAIT: begin
        state <= BRAM_PRIME;
    end

    BRAM_PRIME: begin
        data_tvalid_reg <= 1;
        sample_cnt <= 0;
        state <= STREAM;
    end

    STREAM: begin
        if (stream_fire) begin
            bram_addr <= bram_addr + 1;

            if (sample_cnt == FFT_LEN-1) begin
                data_tvalid_reg <= 0;
                sample_cnt <= 0;
                state <= WAIT_DONE;
            end
            else begin
                sample_cnt <= sample_cnt + 1;
            end
        end
    end

    WAIT_DONE: begin
        if (m_axis_data_tvalid && m_axis_data_tlast) begin
            gap_cnt <= 0;
            state <= INTER_FRAME;
        end
    end

    INTER_FRAME: begin
        if (gap_cnt >= INTER_FRAME_GAP) begin
            bram_addr <= 0;
            state <= BRAM_WAIT;
        end
        else begin
            gap_cnt <= gap_cnt + 1;
        end
    end

    endcase
end

end

//////////////////////////////////////////////////////////// // DATA PATH — ALWAYS TRACK BRAM PIPELINE (fix repeat bug) //////////////////////////////////////////////////////////// always @(posedge aclk) data_tdata_reg <= bram_dout_r;

//////////////////////////////////////////////////////////// // FFT IP //////////////////////////////////////////////////////////// xfft_1 fft_i ( .aclk(aclk), .aresetn(aresetn),

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

.event_frame_started(event_frame_started),
.event_tlast_unexpected(event_tlast_unexpected),
.event_tlast_missing(event_tlast_missing),
.event_status_channel_halt(event_status_channel_halt),
.event_data_in_channel_halt(event_data_in_channel_halt),
.event_data_out_channel_halt(event_data_out_channel_halt)

);

//////////////////////////////////////////////////////////// // OUTPUT SPLIT //////////////////////////////////////////////////////////// assign out_real_data = m_axis_data_tdata[15:0]; assign out_imag_data = m_axis_data_tdata[31:16];

endmodule
