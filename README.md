 `timescale 1ns / 1ps
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

module Top_fft #(

    parameter         DATAWIDTH    = 16,
    parameter         ADDRESSWIDTH = 6,
    parameter         CONFIGWIDTH  = 16,
    parameter         FFT_LEN      = 64,
    parameter         INTER_FRAME_GAP = 20
)( 
    input         wire      aclk,
    input         wire      aresetn,

    /* ================= FFT CONFIG ================= */
    
    output wire [CONFIGWIDTH-1:0] s_axis_config_tdata,
    output wire                   s_axis_config_tvalid,
    input  wire                   s_axis_config_tready,

    /* ================= FFT INPUT ================= */
    
    output wire [2*DATAWIDTH-1:0] s_axis_data_tdata,
    output wire                   s_axis_data_tvalid,
    input  wire                   s_axis_data_tready,
    output wire                   s_axis_data_tlast,

    /* ================= FFT OUTPUT ================= */
    
    input  wire [47:0]            m_axis_data_tdata,
    input  wire                   m_axis_data_tvalid,
    output wire                   m_axis_data_tready,
    input  wire                   m_axis_data_tlast,
    
    // * FFT OUTPUT SPLIT *//
    output wire [46:24] out_real_data, 
    output wire [22:0] out_imag_data,
    

    /* ================= STATUS ================= */
    output wire event_frame_started,
    output wire event_tlast_unexpected,
    output wire event_tlast_missing,
    output wire event_status_channel_halt,
    output wire event_data_in_channel_halt,
    output wire event_data_out_channel_halt
   );

    //*  FSM States *// 
  reg [2:0] state;
  localparam       IDLE          =  3'b000;
  localparam       CFG           =  3'b001;
  localparam       BRAM_WAIT     =  3'b010;                    
  localparam       STREAM        =  3'b011;
  localparam       WAIT_FOR_DONE =  3'b100;
  localparam       INTER_FRAME     =  3'b101;

  // * INTERNAL REGS  * //
  reg [CONFIGWIDTH-1:0]      config_tdata_reg;
  reg                        config_tvalid_reg;
  reg [31:0]                 data_tdata_reg;
  reg                        data_tvalid_reg;
  reg  [ADDRESSWIDTH-1:0]    bram_addr;
 
  reg [7:0]                   sample_cnt;
  reg [7:0]                   gap_cnt;
reg [1:0] vld_delay; // 2-bit shift register for 2-cycle BRAM latency


// * Continuous Assignment *//

assign s_axis_config_tdata  =  config_tdata_reg;
assign s_axis_config_tvalid =  config_tvalid_reg;
assign s_axis_data_tdata   =   data_tdata_reg;
assign s_axis_data_tvalid = data_tvalid_reg;


assign s_axis_data_tlast   =   (sample_cnt == FFT_LEN-1)&& data_tvalid_reg;
assign m_axis_data_tready  =  1'b1;   //Always ready to receive output


//* BRAM INSTANTIATION *//

wire [31:0] bram_dout;
blk_mem_gen_0 u_bram (
  .clka(aclk),    // input wire clka
  .addra(bram_addr),  // input wire [5 : 0] addra
  .douta(bram_dout)  // output wire [31 : 0] douta
);

// ============================================================
// FSM  LOGIC
// ============================================================
always @(posedge aclk) begin
    if (!aresetn) begin
        state <= IDLE;
        bram_addr <= 0;
        sample_cnt <= 0;
        gap_cnt <= 0;
        config_tvalid_reg <= 0;
        data_tvalid_reg <= 0;
        config_tdata_reg <= 16'b0000000010000110;
         vld_delay <= 2'b00;


    end else begin
      // Shift the valid signal to match BRAM latency (Adjust length if needed)
        vld_delay[1] <= vld_delay[0];
        case (state)
           IDLE: begin
            config_tvalid_reg <= 1'b1;
            state <= CFG;
        end

        CFG: begin
            if (s_axis_config_tready && s_axis_config_tvalid)begin
              config_tvalid_reg <= 0;
                state <= BRAM_WAIT;
                bram_addr <= 0;
           end
        end

        // wait 1 cycle for BRAM latency
        BRAM_WAIT: begin
            state <= STREAM;
            data_tvalid_reg <= 1'b1;
        end

        STREAM: begin
            if (s_axis_data_tready) begin
                if (sample_cnt == FFT_LEN-1) begin
                    sample_cnt <= 0;
                    // Instead of data_tvalid_reg <= 0, we clear the pipe input
            vld_delay[0] <= 0;
                    state <= WAIT_FOR_DONE;
                end else begin
                    sample_cnt <= sample_cnt + 1;
                    bram_addr <= bram_addr + 1;
                    vld_delay[0] <= 1; // Start the valid pulse into the pipe
                    data_tdata_reg <= bram_dout;
                 end 
                end else begin
                  data_tdata_reg <= 0;
              end
            end
            WAIT_FOR_DONE: begin
            if (m_axis_data_tlast && m_axis_data_tvalid)begin
                state <= INTER_FRAME;
                gap_cnt <= 0;
             end
            end
             INTER_FRAME: begin
               if(gap_cnt >= INTER_FRAME_GAP)begin
                  state <= BRAM_WAIT;
                  bram_addr <= 0;
                end else begin
                  gap_cnt <= gap_cnt +1;
                  end
                  end
             endcase
           end
        end
        
  always@(*) begin
     data_tdata_reg = bram_dout;
  end
// ============================================================
// FFT IP INSTANCE 
// ============================================================
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
  assign s_axis_data_tvalid = vld_delay[1];
 assign out_real_data = m_axis_data_tdata[46:24];
 assign out_imag_data = m_axis_data_tdata[22:0];
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

wire [46:24] out_real_data;
wire [22:0]  out_imag_data;

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
