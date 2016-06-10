# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 22:28:37 2013

@author: stefan
"""
from collections import OrderedDict
from numpy import sort

def obj2dict(obj):
    copied_dict = OrderedDict()
    try:
        candidate_dict = obj.__dict__
    except AttributeError:
        return obj
    for attr in sort(candidate_dict.keys()):
        vals = candidate_dict[attr]
        if isinstance (vals, list):
            # iterate over list
            collection = []
            for elem in vals:
                collection.append(obj2dict(elem))
            copied_dict[attr] = collection
        else:
            copied_dict[attr] = obj2dict(vals)
    return copied_dict
    
def show_tree(obj):
    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui
    
    app = QtGui.QApplication([])
    tr = obj2dict(obj)
    tree = pg.DataTreeWidget(data=tr)
    tree.show()
    QtGui.QApplication.instance().exec_()


if __name__ == '__main__':
    import BRUKER_classes
    exp_readfid = BRUKER_classes.Experiment('readfidTest', lazy=False)
    show_tree(exp_readfid)