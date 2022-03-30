Examples for BRUKERIO
=====================

Reading Bruker headers and data shouldn't be all that complicated. Here are
some examples on how to do this:

>>> import sarpy
>>> fname = 'sarpy/test/test-data/Moosvi.ii1/2/acqp'
>>> header = sarpy.readJCAMP(fname)

The header is a dictionary with parameter names as the keys:

>>> type(header)
<type 'dict'>

All of the BRUKER header files are in this JCAMP format. So the *method* parameter file can be loaded as well:

>>> fname = 'sarpy/test/test-data/Moosvi.ii1/2/method'
>>> header = sarpy.readJCAMP(fname)
>>> header['ExcPulse']
'(1, 5400, 30, 9.02336512447152, 100, 0, 100, LIB_EXCITATION, < hermite.exc>, 5400, 0.1794, 50, 0.1024, conventional)'

Loading k-space data (BRUKER calls the file fid) is also doable:

>>> fname = 'sarpy/test/test-data/Moosvi.ii1/2/fid'
>>> data = sarpy.readfid(fname)
/tmp/sarpy/sarpy/test/test-data/Moosvi.ii1/2
opening /tmp/sarpy/sarpy/test/test-data/Moosvi.ii1/2/acqp
opening /tmp/sarpy/sarpy/test/test-data/Moosvi.ii1/2/method
ACQ_size = [128, 128, 15], NR=1, ACQ_obj_order=[0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13], EncSteps=[  
   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53
  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71
  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89
  90  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107
 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125
 126 127]
size(fid) = (245760,) ?=? ACQ_size*NR=245760

Function readfid returns another dictionary that bundles the actual data and header information:

>>> data.keys()
['isImage', 'header', 'data']
>>> data['data'].shape
(128, 128, 15, 1)
>>> data['isImage']
False

Because we are now dealing with these dictionaries nested in dictionaries, there is no need to load header data separately before the fid data (unless you don't need the actual data and would rather avoid filling your memory).

>>> data['header']['method']['ExcPulse']
'(1, 5400, 30, 9.02336512447152, 100, 0, 100, LIB_EXCITATION, < hermite.exc>, 5400, 0.1794, 50, 0.1024, conventional)'

Here is a little python string kung-fu to turn this structure into a list of strings for further processing:

>>> data['header']['method']['ExcPulse'].strip('()').split(',')
['1', ' 5400', ' 30', ' 9.02336512447152', ' 100', ' 0', ' 100', ' LIB_EXCITATION', ' < hermite.exc>', ' 5400', ' 0.1794', ' 50', ' 0.1024', ' conventional']

------
The standard way of interacting with BRUKER data is to read in a scan:

>>> scn = sarpy.Scan('OES2P4s13.cZ1/7')
>>> scn
sarpy.BRUKERIO.BRUKER_classes.Scan("OES2P4s13.cZ1/7")
>>> print(scn)
Scan object: Scan("OES2P4s13.cZ1/7")
   protocol: 07_FLASH2D.6sl-OEMRI
   TE = [2.6696]
   TR = [132.9832]
   FA = 40
   FoV = [2.88, 1.92]
   matrix = [192, 64]

All the parameters are accessible through direct attribute access (as well as dictionary access).

>>> scn.method.DATE
'Tue Mar 15 12:51:20 2022 PDT (UT-7h)'
>>> scn.acqp.ACQ_echo_time
[2.6696]




