#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Tue Mar 24 13:31:57 2015
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import time
import datetime
import wx

class top_block(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Top Block")
        _icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
        self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 5e6
        self.center_freq2 = center_freq2 = 10e6
        self.center_freq = center_freq = 1.57542e09
        self.Gain3 = Gain3 = 26
        self.Gain2 = Gain2 = 60
        self.Gain1 = Gain1 = 60
        self.BW = BW = 1.023e6

        ##################################################
        # Blocks
        ##################################################
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
        	self.GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="FFT Plot",
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	"addr0=192.168.10.2",
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0.set_clock_source("gpsdo", 0)
        self.uhd_usrp_source_0.set_time_source("gpsdo", 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_0.set_gain(60, 0)
        self.uhd_usrp_source_0.set_antenna("RX2", 0)
        self.uhd_usrp_source_0.set_bandwidth(BW, 0)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.uhd_usrp_source_0, 0), (self.wxgui_fftsink2_0, 0))


# QT sink close method reimplementation

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)

    def get_center_freq2(self):
        return self.center_freq2

    def set_center_freq2(self, center_freq2):
        self.center_freq2 = center_freq2

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 0)

    def get_Gain3(self):
        return self.Gain3

    def set_Gain3(self, Gain3):
        self.Gain3 = Gain3

    def get_Gain2(self):
        return self.Gain2

    def set_Gain2(self, Gain2):
        self.Gain2 = Gain2

    def get_Gain1(self):
        return self.Gain1

    def set_Gain1(self, Gain1):
        self.Gain1 = Gain1

    def get_BW(self):
        return self.BW

    def set_BW(self, BW):
        self.BW = BW
        self.uhd_usrp_source_0.set_bandwidth(self.BW, 0)

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"
    parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
    (options, args) = parser.parse_args()
    tb = top_block()
    tb.Start(True)
    gpsdo = tb.uhd_usrp_source_0.get_mboard_sensor("gps_time")
    gpsdo = datetime.datetime.fromtimestamp(gpsdo).strftime('%Y-%m-%d %H:%M:%S.%f')
        # check gpsdo status
    print(tb.uhd_usrp_source_0.get_mboard_sensor("lo_locked"))
    print(tb.uhd_usrp_source_0.get_mboard_sensor("gps_locked"))
    print(tb.uhd_usrp_source_0.get_mboard_sensor("gps_gpgga"))
    print(tb.uhd_usrp_source_0.get_mboard_sensor("gps_gprmc"))
    print(tb.uhd_usrp_source_0.get_mboard_sensor("ref_locked"))
    print gpsdo
    print datetime.datetime.now().time()
    tb.Wait()

