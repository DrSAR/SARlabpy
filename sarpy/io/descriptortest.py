# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 13:19:35 2013

@author: stefan
"""

class D(object):
    def __get__(self, instance, value):
        print('getting variable')
        return 'forty two'
    def __set__(*args): raise AttributeError
    
    
class C(object):
    a = D()
    
x = C()
print(x.a)
#x.a = 42

print('-'*40)

    

class Dnew(object):
    def __init__(self):
        self.__yet_loaded = False
        print('running Dnew __init__')
    def __get__(self, instance, value):
        print('getting variable')
        if not self.__yet_loaded:
            self._load()
        return self.data
    def __set__(*args): raise AttributeError
    def _load(self):
        print('loading now')
        self.data=[42,42,42]
        self.__yet_loaded = True
    
    
class Cnew(object):
    data = Dnew()
    def __init__(self):
        print('this is Cnew __init__')
        
print('Cnew not yet instantiated')
y = Cnew()
print('instance exists, data not yet loaded')
print(y.data)
print('let us do this again')
print(y.data)