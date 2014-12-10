from __future__ import division

__author__ = 'volod_kuzn'
__version__ = '0.1'

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import struct


class Axis():
    def __init__(self, name, units, start, stop):
        self.name = name
        self.units = units
        self.start = start
        self.stop = stop
        self.step = 0

    def full_name(self):
        return "%s, %s" % (self.name, self.units)

    def __str__(self):
        return "%s, %s from %u to %u with step %f" % (self.name, self.units, self.start, self.stop, self.step)


def show():
    plt.show()


def get_pixel_coord(axis, pixel):
    """
    :param axis: This is actual axis parameter
    :param pixel: This is pixel coordinate in data units
    :return: Pixel coordinate in pixel units
    """
    return int((pixel - axis.start) / axis.step)


class Spectrum2D():
    def __init__(self, file_name, y_axis=None):
        if not (file_name.endswith('_.SPE')):
            raise ValueError("Could not create spectrum from not _.SPE file")

        self.name = file_name
        self.x_axis = None
        if y_axis is None:
            self.y_axis = Axis("", "", 0, 0)
        else:
            self.y_axis = y_axis

        self._read_spectrum(file_name)
        self.lum = self.raw_lum.copy()
        self.fig1d = None
        self.fig2d = None
        self.ax2d = None
        self.ax1d = None
        self.image2d = None
        self.colorbar = None
        self.rec_select = None

    def _read_spectrum(self, file_name):
        self._raw_data = read_spe(file_name)
        self.raw_lum = self._raw_data['data'][0]
        p = self._raw_data['XCALIB']["polynom_coeff"]
        self.calib_polynom = np.array([p[2], p[1], p[0]])
        self.wavelength = np.polyval(self.calib_polynom, xrange(1, 1 + len(self.raw_lum[0])))
        self.x_axis = Axis("Wavelength", "nm", self.wavelength[0], self.wavelength[-1])
        self.x_axis.step = (self.x_axis.stop - self.x_axis.start) / (len(self.raw_lum[0]) - 1)
        self.y_axis.step = (self.y_axis.stop - self.y_axis.start) / (len(self.raw_lum) - 1)

    @staticmethod
    def _construct_extent(x_axis, y_axis):
        extent = list()
        extent.append(x_axis.start - x_axis.step / 2)
        extent.append(x_axis.stop + x_axis.step / 2)
        extent.append(y_axis.start - y_axis.step / 2)
        extent.append(y_axis.stop + y_axis.step / 2)
        return extent

    def plot(self, figsize=(10, 10)):
        self.fig2d = plt.figure(num=1, figsize=figsize)
        self.ax2d = self.fig2d.add_subplot(111)
        extent = self._construct_extent(self.x_axis, self.y_axis)
        self.image2d = self.ax2d.imshow(self.lum, aspect='auto', extent=extent, interpolation='nearest')
        self.ax2d.set_xlabel(self.x_axis.full_name())
        self.ax2d.set_ylabel(self.y_axis.full_name())
        self.colorbar = self.fig2d.colorbar(self.image2d)

        def onselect(eclick, erelease):
            """

            eclick and erelease are matplotlib events at press and release

            """
            # print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)
            # print ' startposition : (%u, %u)' % (get_pixel_coord(self.x_axis, eclick.xdata),
            # get_pixel_coord(self.y_axis, eclick.ydata))
            # print self.ax2d.images
            # print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)
            # print ' used button   : ', eclick.button
            # print ' endposition : (%u, %u)' % (get_pixel_coord(self.x_axis, erelease.xdata),
            # get_pixel_coord(self.y_axis, erelease.ydata))

            x = [0, 0]
            x[0] = get_pixel_coord(self.x_axis, eclick.xdata)
            x[1] = get_pixel_coord(self.x_axis, erelease.xdata)
            x.sort()
            y = [0, 0]
            y[0] = get_pixel_coord(self.y_axis, eclick.ydata)
            y[1] = get_pixel_coord(self.y_axis, erelease.ydata)
            y.sort()

            tmp = self.lum[-y[1]: -y[0], x[0]: x[1]]
            # print("Tmp slice: %u: %u, %u: %u" % (x[0], x[1], y[0], y[1]))
            # print("Min: %u, max: %u" % (tmp.min(), tmp.max()))
            # print("All min: %u, max: %u" % (self.lum.min(), self.lum.max()))

            self.image2d.norm.vmin = tmp.min()
            self.image2d.norm.vmax = tmp.max()
            self.image2d.changed()
            # self.ax2d.imshow(self.lum, aspect='auto', extent=extent, interpolation='nearest', vmin=tmp.min(),
            # vmax=tmp.max())
            self.fig2d.canvas.draw()

        self.rec_select = RectangleSelector(self.ax2d, onselect=onselect, drawtype='box', minspanx=5, minspany=5,
                                            spancoords=u'pixels')

    def plot_slice(self, spec_num=0, figsize=(10, 10)):
        lum = self.lum[spec_num]
        self.fig1d = plt.figure(num=2, figsize=figsize)
        self.ax1d = self.fig1d.add_subplot(111)
        self.ax1d.plot(self.wavelength, lum)
        self.ax1d.set_xlabel(self.x_axis.full_name())
        self.ax1d.set_ylabel("Intensity, a. u.")


# from Kasey
# http://people.seas.harvard.edu/~krussell/html-tutorial/_modules/winspec.html
# noinspection PyPep8
def read_spe(spefilename, verbose=False):
    """
    Read a binary PI SPE file into a python dictionary

    Inputs:

        spefilename --  string specifying the name of the SPE file to be read
        verbose     --  boolean print debug statements (True) or not (False)

        Outputs
        spedict

            python dictionary containing header and data information
            from the SPE file
            Content of the dictionary is:
            spedict = {'data':[],    # a list of 2D numpy arrays, one per image
            'IGAIN':pimax_gain,
            'EXPOSURE':exp_sec,
            'SPEFNAME':spefilename,
            'OBSDATE':date,
            'CHIPTEMP':detectorTemperature
            }

    I use the struct module to unpack the binary SPE data.
    Some useful formats for struct.unpack_from() include:
    fmt   c type          python
    c     char            string of length 1
    s     char[]          string (Ns is a string N characters long)
    h     short           integer
    H     unsigned short  integer
    l     long            integer
    f     float           float
    d     double          float

    The SPE file defines new c types including:
        BYTE  = unsigned char
        WORD  = unsigned short
        DWORD = unsigned long


    Example usage:
    Given an SPE file named test.SPE, you can read the SPE data into
    a python dictionary named spedict with the following:
    """

    # open SPE file as binary input
    spe = open(spefilename, "rb")

    # Header length is a fixed number
    n_bytes_in_header = 4100

    # Read the entire header
    header = spe.read(n_bytes_in_header)

    # version of WinView used
    sw_version = struct.unpack_from("16s", header, offset=688)[0]

    # version of header used
    # Eventually, need to adjust the header unpacking
    # based on the header_version.
    header_version = struct.unpack_from("f", header, offset=1992)[0]

    # which camera controller was used?
    controller_version = struct.unpack_from("h", header, offset=0)[0]
    if verbose:
        print "sw_version         = ", sw_version
        print "header_version     = ", header_version
        print "controller_version = ", controller_version

    # Date of the observation
    # (format is DDMONYYYY  e.g. 27Jan2009)
    date = struct.unpack_from("9s", header, offset=20)[0]

    # Exposure time (float)
    exp_sec = struct.unpack_from("f", header, offset=10)[0]

    # Intensifier gain
    pimax_gain = struct.unpack_from("h", header, offset=148)[0]

    # Not sure which "gain" this is
    gain = struct.unpack_from("H", header, offset=198)[0]

    # Data type (0=float, 1=long integer, 2=integer, 3=unsigned int)
    data_type = struct.unpack_from("h", header, offset=108)[0]

    comments = struct.unpack_from("400s", header, offset=200)[0]

    # CCD Chip Temperature (Degrees C)
    detector_temperature = struct.unpack_from("f", header, offset=36)[0]

    # The following get read but are not used
    # (this part is only lightly tested...)
    analog_gain = struct.unpack_from("h", header, offset=4092)[0]
    noscan = struct.unpack_from("h", header, offset=34)[0]
    pimax_used = struct.unpack_from("h", header, offset=144)[0]
    # pimaxMode = struct.unpack_from("h", header, offset=146)[0]

    # here's from Kasey
    # int avgexp 2 number of accumulations per scan (why don't they call this "accumulations"?)
    # this isn't actually accumulations, so fix it...
    accumulations = struct.unpack_from("h", header, offset=668)[0]
    if accumulations == -1:
        # if > 32767, set to -1 and
        # see lavgexp below (668)
        # accumulations = struct.unpack_from("l", header, offset=668)[0]
        # or should it be DWORD, NumExpAccums (1422): Number of Time experiment accumulated
        accumulations = struct.unpack_from("l", header, offset=1422)[0]

    """Start of X Calibration Structure (although I added things to it that I thought were relevant,
       like the center wavelength..."""
    xcalib = {'SpecAutoSpectroMode': bool(struct.unpack_from("h", header, offset=70)[0]),
              'SpecCenterWlNm': struct.unpack_from("f", header, offset=72)[0],
              'SpecGlueFlag': bool(struct.unpack_from("h", header, offset=76)[0]),
              'SpecGlueStartWlNm': struct.unpack_from("f", header, offset=78)[0],
              'SpecGlueEndWlNm': struct.unpack_from("f", header, offset=82)[0],
              'SpecGlueMinOvrlpNm': struct.unpack_from("f", header, offset=86)[0],
              'SpecGlueFinalResNm': struct.unpack_from("f", header, offset=90)[0],
              'background_applied': struct.unpack_from("h", header, offset=150)[0]}

    # short   BackGrndApplied              150  1 if background subtraction done
    background_applied = False
    if xcalib['background_applied'] == 1:
        background_applied = True

    # float   SpecGrooves                  650  Spectrograph Grating Grooves
    xcalib['SpecGrooves'] = struct.unpack_from("f", header, offset=650)[0]

    # short   flat_field_applied             706  1 if flat field was applied.
    xcalib['flat_field_applied'] = struct.unpack_from("h", header, offset=706)[0]
    flat_field_applied = False
    if xcalib['flat_field_applied'] == 1:
        flat_field_applied = True

    # double offset # 3000 offset for absolute data scaling */
    xcalib['offset'] = struct.unpack_from("d", header, offset=3000)[0]

    # double factor # 3008 factor for absolute data scaling */
    xcalib['factor'] = struct.unpack_from("d", header, offset=3008)[0]

    # char current_unit # 3016 selected scaling unit */
    xcalib['current_unit'] = struct.unpack_from("c", header, offset=3016)[0]

    # char reserved1 # 3017 reserved */
    xcalib['reserved1'] = struct.unpack_from("c", header, offset=3017)[0]

    # char string[40] # 3018 special string for scaling */
    xcalib['string'] = struct.unpack_from("40s", header, offset=3018)

    # char reserved2[40] # 3058 reserved */
    xcalib['reserved2'] = struct.unpack_from("40s", header, offset=3058)

    # char calib_valid # 3098 flag if calibration is valid */
    xcalib['calib_valid'] = struct.unpack_from("c", header, offset=3098)[0]

    # char input_unit # 3099 current input units for */
    xcalib['input_unit'] = struct.unpack_from("c", header, offset=3099)[0]
    """/* "calib_value" */"""

    # char polynom_unit # 3100 linear UNIT and used */
    xcalib['polynom_unit'] = struct.unpack_from("c", header, offset=3100)[0]
    """/* in the "polynom_coeff" */"""

    # char polynom_order # 3101 ORDER of calibration POLYNOM */
    xcalib['polynom_order'] = struct.unpack_from("c", header, offset=3101)[0]

    # char calib_count # 3102 valid calibration data pairs */
    xcalib['calib_count'] = struct.unpack_from("c", header, offset=3102)[0]

    # double pixel_position[10];/* 3103 pixel pos. of calibration data */
    xcalib['pixel_position'] = struct.unpack_from("10d", header, offset=3103)

    # double calib_value[10] # 3183 calibration VALUE at above pos */
    xcalib['calib_value'] = struct.unpack_from("10d", header, offset=3183)

    # double polynom_coeff[6] # 3263 polynom COEFFICIENTS */
    xcalib['polynom_coeff'] = struct.unpack_from("6d", header, offset=3263)

    # double laser_position # 3311 laser wavenumber for relativ WN */
    xcalib['laser_position'] = struct.unpack_from("d", header, offset=3311)[0]

    # char reserved3 # 3319 reserved */
    xcalib['reserved3'] = struct.unpack_from("c", header, offset=3319)[0]

    # unsigned char new_calib_flag # 3320 If set to 200, valid label below */
    # xcalib['calib_value'] = struct.unpack_from("BYTE", header, offset=3320)[0] # how to do this?

    # char calib_label[81] # 3321 Calibration label (NULL term'd) */
    xcalib['calib_label'] = struct.unpack_from("81s", header, offset=3321)

    # char expansion[87] # 3402 Calibration Expansion area */
    xcalib['expansion'] = struct.unpack_from("87s", header, offset=3402)
    # end of Kasey's addition

    if verbose:
        print "date      = [" + date + "]"
        print "exp_sec   = ", exp_sec
        print "pimax_gain = ", pimax_gain
        print "gain (?)  = ", gain
        print "data_type = ", data_type
        print "comments  = [" + comments + "]"
        print "analog_gain = ", analog_gain
        print "noscan = ", noscan
        print "detectorTemperature [C] = ", detector_temperature
        print "pimax_used = ", pimax_used

    # Determine the data type format string for
    # upcoming struct.unpack_from() calls
    if data_type == 0:
        # float (4 bytes)
        data_type_str = "f"  # untested
        bytes_per_pixel = 4
        dtype = "float32"
    elif data_type == 1:
        # long (4 bytes)
        data_type_str = "l"  # untested
        bytes_per_pixel = 4
        dtype = "int32"
    elif data_type == 2:
        # short (2 bytes)
        data_type_str = "h"  # untested
        bytes_per_pixel = 2
        dtype = "int32"
    elif data_type == 3:
        # unsigned short (2 bytes)
        data_type_str = "H"  # 16 bits in python on intel mac
        bytes_per_pixel = 2
        dtype = "int32"  # for numpy.array().
        # other options include:
        # IntN, UintN, where N = 8,16,32 or 64
        # and Float32, Float64, Complex64, Complex128
        # but need to verify that pyfits._ImageBaseHDU.ImgCode cna handle it
        # right now, ImgCode must be float32, float64, int16, int32, int64 or uint8
    else:
        print "unknown data type"
        print "returning..."
        sys.exit()

    # Number of pixels on x-axis and y-axis
    nx = struct.unpack_from("H", header, offset=42)[0]
    ny = struct.unpack_from("H", header, offset=656)[0]

    # Number of image frames in this SPE file
    nframes = struct.unpack_from("l", header, offset=1446)[0]

    if verbose:
        print "nx, ny, nframes = ", nx, ", ", ny, ", ", nframes

    npixels = nx * ny
    npixstr = str(npixels)
    fmt_str = npixstr + data_type_str
    if verbose:
        print "fmt_str = ", fmt_str

    # How many bytes per image?
    n_bytes_per_frame = npixels * bytes_per_pixel
    if verbose:
        print "n_bytes_per_frame = ", n_bytes_per_frame

    # Create a dictionary that holds some header information
    # and contains a placeholder for the image data
    spedict = {'data': [],  # can have more than one image frame per SPE file
               'IGAIN': pimax_gain,
               'EXPOSURE': exp_sec,
               'SPEFNAME': spefilename,
               'OBSDATE': date,
               'CHIPTEMP': detector_temperature,
               'COMMENTS': comments,
               'XCALIB': xcalib,
               'ACCUMULATIONS': accumulations,
               'FLATFIELD': flat_field_applied,
               'BACKGROUND': background_applied
               }

    # Now read in the image data
    # Loop over each image frame in the image
    if verbose:
        print "Reading image frames number ",
    for ii in range(nframes):
        data = spe.read(n_bytes_per_frame)
        if verbose:
            print ii, " ",

        # read pixel values into a 1-D numpy array. the "=" forces it to use
        # standard python datatype size (4bytes for 'l') rather than native
        # (which on 64bit is 8bytes for 'l', for example).
        # See http://docs.python.org/library/struct.html
        data_arr = np.array(struct.unpack_from("=" + fmt_str, data, offset=0),
                            dtype=dtype)

        # Resize array to nx by ny pixels
        # notice order... (y,x)
        data_arr.resize((ny, nx))
        # print data_arr.shape

        # Push this image frame data onto the end of the list of images
        # but first cast the datatype to float (if it's not already)
        # this isn't necessary, but shouldn't hurt and could save me
        # from doing integer math when i really meant floating-point...
        # noinspection PyTypeChecker
        spedict['data'].append(data_arr.astype(float))

    if verbose:
        print ""

    return spedict

    ###############################################################################
    ###############################################################################
    #           Description of the header structure used to create piUtils      ###
    ###############################################################################
    ###############################################################################
    #
    # WINHEAD.TXT
    #
    # $Date: 3/23/04 11:36 $
    #
    # Header Structure For WinView/WinSpec (WINX) Files
    #
    # The current data file used for WINX files consists of a 4100 (1004 Hex)
    # byte header followed by the data.
    #
    # Beginning with Version 2.5, many more items were added to the header to
    # make it a complete as possible record of the data collection.  This includes
    # spectrograph and pulser information.  Much of these additions were accomplished
    # by recycling old information which had not been used in many versions.
    # All data files created under previous 2.x versions of WinView/WinSpec CAN
    # still be read correctly.  HOWEVER, files created under the new versions
    # (2.5 and higher) CANNOT be read by previous versions of WinView/WinSpec
    # OR by the CSMA software package.
    #
    #
    # ***************************************************
    #
    # Decimal Byte
    # Offset
    # -----------
    # short   ControllerVersion              0  Hardware Version
    # short   LogicOutput                    2  Definition of Output BNC
    # WORD    AmpHiCapLowNoise               4  Amp Switching Mode
    # WORD    xDimDet                        6  Detector x dimension of chip.
    #  short   mode                           8  timing mode
    #  float   exp_sec                       10  alternitive exposure, in sec.
    #  short   VChipXdim                     14  Virtual Chip X dim
    #  short   VChipYdim                     16  Virtual Chip Y dim
    #  WORD    yDimDet                       18  y dimension of CCD or detector.
    #  char    date[DATEMAX]                 20  date
    #  short   VirtualChipFlag               30  On/Off
    #  char    Spare_1[2]                    32
    #  short   noscan                        34  Old number of scans - should always be -1
    #  float   DetTemperature                36  Detector Temperature Set
    #  short   DetType                       40  CCD/DiodeArray type
    #  WORD    xdim                          42  actual # of pixels on x axis
    #  short   stdiode                       44  trigger diode
    #  float   DelayTime                     46  Used with Async Mode
    #  WORD    ShutterControl                50  Normal, Disabled Open, Disabled Closed
    #  short   AbsorbLive                    52  On/Off
    #  WORD    AbsorbMode                    54  Reference Strip or File
    #  short   CanDoVirtualChipFlag          56  T/F Cont/Chip able to do Virtual Chip
    #  short   ThresholdMinLive              58  On/Off
    #  float   ThresholdMinVal               60  Threshold Minimum Value
    #  short   ThresholdMaxLive              64  On/Off
    #  float   ThresholdMaxVal               66  Threshold Maximum Value
    #  short   SpecAutoSpectroMode           70  T/F Spectrograph Used
    #  float   SpecCenterWlNm                72  Center Wavelength in Nm
    #  short   SpecGlueFlag                  76  T/F File is Glued
    #  float   SpecGlueStartWlNm             78  Starting Wavelength in Nm
    #  float   SpecGlueEndWlNm               82  Starting Wavelength in Nm
    #  float   SpecGlueMinOvrlpNm            86  Minimum Overlap in Nm
    #  float   SpecGlueFinalResNm            90  Final Resolution in Nm
    #  short   PulserType                    94  0=None, PG200=1, PTG=2, DG535=3
    #  short   CustomChipFlag                96  T/F Custom Chip Used
    #  short   XPrePixels                    98  Pre Pixels in X direction
    #  short   XPostPixels                  100  Post Pixels in X direction
    #  short   YPrePixels                   102  Pre Pixels in Y direction
    #  short   YPostPixels                  104  Post Pixels in Y direction
    #  short   asynen                       106  asynchron enable flag  0 = off
    #  short   datatype                     108  experiment datatype
    #                                             0 =   float (4 bytes)
    #                                             1 =   long (4 bytes)
    #                                             2 =   short (2 bytes)
    #                                             3 =   unsigned short (2 bytes)
    #  short   PulserMode                   110  Repetitive/Sequential
    #  WORD    PulserOnChipAccums           112  Num PTG On-Chip Accums
    #  DWORD   PulserRepeatExp              114  Num Exp Repeats (Pulser SW Accum)
    #  float   PulseRepWidth                118  Width Value for Repetitive pulse (usec)
    #  float   PulseRepDelay                122  Width Value for Repetitive pulse (usec)
    #  float   PulseSeqStartWidth           126  Start Width for Sequential pulse (usec)
    #  float   PulseSeqEndWidth             130  End Width for Sequential pulse (usec)
    #  float   PulseSeqStartDelay           134  Start Delay for Sequential pulse (usec)
    #  float   PulseSeqEndDelay             138  End Delay for Sequential pulse (usec)
    #  short   PulseSeqIncMode              142  Increments: 1=Fixed, 2=Exponential
    #  short   PImaxUsed                    144  PI-Max type controller flag
    #  short   PImaxMode                    146  PI-Max mode
    #  short   PImaxGain                    148  PI-Max Gain
    #  short   BackGrndApplied              150  1 if background subtraction done
    #  short   PImax2nsBrdUsed              152  T/F PI-Max 2ns Board Used
    #  WORD    minblk                       154  min. # of strips per skips
    #  WORD    numminblk                    156  # of min-blocks before geo skps
    #  short   SpecMirrorLocation[2]        158  Spectro Mirror Location, 0=Not Present
    #  short   SpecSlitLocation[4]          162  Spectro Slit Location, 0=Not Present
    #  short   CustomTimingFlag             170  T/F Custom Timing Used
    #  char    ExperimentTimeLocal[TIMEMAX] 172  Experiment Local Time as hhmmss\0
    #  char    ExperimentTimeUTC[TIMEMAX]   179  Experiment UTC Time as hhmmss\0
    #  short   ExposUnits                   186  User Units for Exposure
    #  WORD    ADCoffset                    188  ADC offset
    #  WORD    ADCrate                      190  ADC rate
    #  WORD    ADCtype                      192  ADC type
    #  WORD    ADCresolution                194  ADC resolution
    #  WORD    ADCbitAdjust                 196  ADC bit adjust
    #  WORD    gain                         198  gain
    #  char    Comments[5][COMMENTMAX]      200  File Comments
    #  WORD    geometric                    600  geometric ops: rotate 0x01,
    #                                             reverse 0x02, flip 0x04
    #  char    xlabel[LABELMAX]             602  intensity display string
    #  WORD    cleans                       618  cleans
    #  WORD    NumSkpPerCln                 620  number of skips per clean.
    #  short   SpecMirrorPos[2]             622  Spectrograph Mirror Positions
    #  float   SpecSlitPos[4]               626  Spectrograph Slit Positions
    #  short   AutoCleansActive             642  T/F
    #  short   UseContCleansInst            644  T/F
    #  short   AbsorbStripNum               646  Absorbance Strip Number
    #  short   SpecSlitPosUnits             648  Spectrograph Slit Position Units
    #  float   SpecGrooves                  650  Spectrograph Grating Grooves
    #  short   srccmp                       654  number of source comp. diodes
    #  WORD    ydim                         656  y dimension of raw data.
    #  short   scramble                     658  0=scrambled,1=unscrambled
    #  short   ContinuousCleansFlag         660  T/F Continuous Cleans Timing Option
    #  short   ExternalTriggerFlag          662  T/F External Trigger Timing Option
    #  long    lnoscan                      664  Number of scans (Early WinX)
    #  long    lavgexp                      668  Number of Accumulations
    #  float   ReadoutTime                  672  Experiment readout time
    #  short   TriggeredModeFlag            676  T/F Triggered Timing Option
    #  char    Spare_2[10]                  678
    #  char    sw_version[FILEVERMAX]       688  Version of SW creating this file
    #  short   type                         704   1 = new120 (Type II)
    #                                             2 = old120 (Type I )
    #                                             3 = ST130
    #                                             4 = ST121
    #                                             5 = ST138
    #                                             6 = DC131 (PentaMax)
    #                                             7 = ST133 (MicroMax/SpectroMax)
    #                                             8 = ST135 (GPIB)
    #                                             9 = VICCD
    #                                            10 = ST116 (GPIB)
    #                                            11 = OMA3 (GPIB)
    #                                            12 = OMA4
    # short   flat_field_applied             706  1 if flat field was applied.
    #  char    Spare_3[16]                  708
    #  short   kin_trig_mode                724  Kinetics Trigger Mode
    #  char    dlabel[LABELMAX]             726  Data label.
    #  char    Spare_4[436]                 742
    #  char    PulseFileName[HDRNAMEMAX]   1178  Name of Pulser File with
    #                                             Pulse Widths/Delays (for Z-Slice)
    #  char    AbsorbFileName[HDRNAMEMAX]  1298 Name of Absorbance File (if File Mode)
    #  DWORD   NumExpRepeats               1418  Number of Times experiment repeated
    #  DWORD   NumExpAccums                1422  Number of Time experiment accumulated
    #  short   YT_Flag                     1426  Set to 1 if this file contains YT data
    #  float   clkspd_us                   1428  Vert Clock Speed in micro-sec
    #  short   HWaccumFlag                 1432  set to 1 if accum done by Hardware.
    #  short   StoreSync                   1434  set to 1 if store sync used
    #  short   BlemishApplied              1436  set to 1 if blemish removal applied
    #  short   CosmicApplied               1438  set to 1 if cosmic ray removal applied
    #  short   CosmicType                  1440  if cosmic ray applied, this is type
    #  float   CosmicThreshold             1442  Threshold of cosmic ray removal.
    #  long    NumFrames                   1446  number of frames in file.
    #  float   MaxIntensity                1450  max intensity of data (future)
    #  float   MinIntensity                1454  min intensity of data (future)
    #  char    ylabel[LABELMAX]            1458  y axis label.
    #  WORD    ShutterType                 1474  shutter type.
    #  float   shutterComp                 1476  shutter compensation time.
    #  WORD    readoutMode                 1480  readout mode, full,kinetics, etc
    #  WORD    WindowSize                  1482  window size for kinetics only.
    #  WORD    clkspd                      1484  clock speed for kinetics & frame transfer
    #  WORD    interface_type              1486  computer interface
    #                                             (isa-taxi, pci, eisa, etc.)
    #  short   NumROIsInExperiment         1488  May be more than the 10 allowed in
    #                                             this header (if 0, assume 1)
    #  char    Spare_5[16]                 1490
    #  WORD    controllerNum               1506  if multiple controller system will
    #                                             have controller number data came from.
    #                                             this is a future item.
    #  WORD    SWmade                      1508  Which software package created this file
    #  short   NumROI                      1510  number of ROIs used. if 0 assume 1.
    #
    #
    # ------------------------------------------------------------------------------
    #
    #                        ROI entries   (1512 - 1631)
    #
    #  struct ROIinfo
    #  {
    #    WORD  startx                            left x start value.
    #    WORD  endx                              right x value.
    #    WORD  groupx                            amount x is binned/grouped in hw.
    #    WORD  starty                            top y start value.
    #    WORD  endy                              bottom y value.
    #    WORD  groupy                            amount y is binned/grouped in hw.
    #  } ROIinfoblk[ROIMAX]
    #                                            ROI Starting Offsets:
    #
    #                                              ROI  1  =  1512
    #                                              ROI  2  =  1524
    #                                              ROI  3  =  1536
    #                                              ROI  4  =  1548
    #                                              ROI  5  =  1560
    #                                              ROI  6  =  1572
    #                                              ROI  7  =  1584
    #                                              ROI  8  =  1596
    #                                              ROI  9  =  1608
    #                                              ROI 10  =  1620
    #
    # ------------------------------------------------------------------------------
    #
    #  char    FlatField[HDRNAMEMAX]       1632  Flat field file name.
    #  char    background[HDRNAMEMAX]      1752  background sub. file name.
    #  char    blemish[HDRNAMEMAX]         1872  blemish file name.
    #  float   file_header_ver             1992  version of this file header
    #  char    YT_Info[1000]               1996-2995  Reserved for YT information
    #  long    WinView_id                  2996  == 0x01234567L if file created by WinX
    #
    # ------------------------------------------------------------------------------
    #
    #                        START OF X CALIBRATION STRUCTURE (3000 - 3488)
    #
    #  double  offset                      3000  offset for absolute data scaling
    #  double  factor                      3008  factor for absolute data scaling
    #  char    current_unit                3016  selected scaling unit
    #  char    reserved1                   3017  reserved
    #  char    string[40]                  3018  special string for scaling
    #  char    reserved2[40]               3058  reserved
    #  char    calib_valid                 3098  flag if calibration is valid
    #  char    input_unit                  3099  current input units for
    #                                            "calib_value"
    #  char    polynom_unit                3100  linear UNIT and used
    #                                            in the "polynom_coeff"
    #  char    polynom_order               3101  ORDER of calibration POLYNOM
    #  char    calib_count                 3102  valid calibration data pairs
    #  double  pixel_position[10]          3103  pixel pos. of calibration data
    #  double  calib_value[10]             3183  calibration VALUE at above pos
    #  double  polynom_coeff[6]            3263  polynom COEFFICIENTS
    #  double  laser_position              3311  laser wavenumber for relativ WN
    #  char    reserved3                   3319  reserved
    #  BYTE    new_calib_flag              3320  If set to 200, valid label below
    #  char    calib_label[81]             3321  Calibration label (NULL term'd)
    #  char    expansion[87]               3402  Calibration Expansion area
    #
    # ------------------------------------------------------------------------------
    #
    #                        START OF Y CALIBRATION STRUCTURE (3489 - 3977)
    #
    #  double  offset                      3489  offset for absolute data scaling
    #  double  factor                      3497  factor for absolute data scaling
    #  char    current_unit                3505  selected scaling unit
    #  char    reserved1                   3506  reserved
    #  char    string[40]                  3507  special string for scaling
    #  char    reserved2[40]               3547  reserved
    #  char    calib_valid                 3587  flag if calibration is valid
    #  char    input_unit                  3588  current input units for
    #                                            "calib_value"
    #  char    polynom_unit                3589  linear UNIT and used
    #                                            in the "polynom_coeff"
    #  char    polynom_order               3590  ORDER of calibration POLYNOM
    #  char    calib_count                 3591  valid calibration data pairs
    #  double  pixel_position[10]          3592  pixel pos. of calibration data
    #  double  calib_value[10]             3672  calibration VALUE at above pos
    #  double  polynom_coeff[6]            3752  polynom COEFFICIENTS
    #  double  laser_position              3800  laser wavenumber for relativ WN
    #  char    reserved3                   3808  reserved
    #  BYTE    new_calib_flag              3809  If set to 200, valid label below
    #  char    calib_label[81]             3810  Calibration label (NULL term'd)
    #  char    expansion[87]               3891  Calibration Expansion area
    #
    #                         END OF CALIBRATION STRUCTURES
    #
    # ------------------------------------------------------------------------------
    #
    #  char    Istring[40]                 3978  special intensity scaling string
    #  char    Spare_6[25]                 4018
    #  BYTE    SpecType                    4043  spectrometer type (acton, spex, etc.)
    #  BYTE    SpecModel                   4044  spectrometer model (type dependent)
    #  BYTE    PulseBurstUsed              4045  pulser burst mode on/off
    #  DWORD   PulseBurstCount             4046  pulser triggers per burst
    #  double  ulseBurstPeriod             4050  pulser burst period (in usec)
    #  BYTE    PulseBracketUsed            4058  pulser bracket pulsing on/off
    #  BYTE    PulseBracketType            4059  pulser bracket pulsing type
    #  double  PulseTimeConstFast          4060  pulser slow exponential time constant (in usec)
    #  double  PulseAmplitudeFast          4068  pulser fast exponential amplitude constant
    #  double  PulseTimeConstSlow          4076  pulser slow exponential time constant (in usec)
    #  double  PulseAmplitudeSlow          4084  pulser slow exponential amplitude constant
    #  short   AnalogGain;                 4092  analog gain
    #  short   AvGainUsed                  4094  avalanche gain was used
    #  short   AvGain                      4096  avalanche gain value
    #  short   lastvalue                   4098  Always the LAST value in the header
    #
    #                         END OF HEADER
    #
    # ------------------------------------------------------------------------------
    #
    #                                      4100  Start of Data
    #
    #
    #
    #        ***************************** E.O.F.*****************************
    #
    #
    #
    #
    #  Definitions of array sizes:
    #  ---------------------------
    #
    #    HDRNAMEMAX  = 120     Max char str length for file name
    #    USERINFOMAX = 1000    User information space
    #    COMMENTMAX  = 80      User comment string max length (5 comments)
    #    LABELMAX    = 16      Label string max length
    #    FILEVERMAX  = 16      File version string max length
    #    DATEMAX     = 10      String length of file creation date string as ddmmmyyyy\0
    #    ROIMAX      = 10      Max size of roi array of structures
    #    TIMEMAX     = 7       Max time store as hhmmss\0
    #
    #
    #
    #  Custom Data Types used in the structure:
    #  ----------------------------------------
    #
    #    BYTE = unsigned char
    #    WORD = unsigned short
    #    DWORD = unsigned long
    #
    #
    #
    #
    #  READING DATA:
    #  -------------
    #
    #    The data follows the header beginning at offset 4100.
    #
    #    Data is stored as sequential points.
    #
    #    The X, Y and Frame dimensions are determined by the header.
    #
    #      The X dimension of the stored data is in "xdim" ( Offset 42  ).
    #      The Y dimension of the stored data is in "ydim" ( Offset 656 ).
    #      The number of frames of data stored is in "NumFrames" ( Offset 1446 ).
    #
    #    The size of a frame (in bytes) is:
    #
    #      One frame size = xdim x ydim x (datatype Offset 108)
    #
    ###############################################################################
    ###############################################################################