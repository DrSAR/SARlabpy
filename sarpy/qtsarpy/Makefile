# create the automated UI code view PySide or Qt uic converter

qtsarpy_mainwindow_auto.py:: UIC=$(shell python test_pyside.py)

qtsarpy_mainwindow_auto.py: ui/qtsarpy_mainwindow.ui
	$(UIC) ui/qtsarpy_mainwindow.ui -o qtsarpy_mainwindow_auto.py
