#!/usr/bin/env python

try:
    import PyQt4
except ImportError:
    try:
        import PySide
    except ImportError:
        print('neither PyQt4 nor PySide found')
    else:
        print('pyside-uic')
else:
    print('pyuic4-2.7') 
