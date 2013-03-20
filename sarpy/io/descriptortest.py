# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 13:19:35 2013

@author: stefan
"""

class dict_plus(dict):
    def __key__(self, key, value):
        print('write access to dict')
        self.x=value

class lazy_property_desc(object):

    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__
        self.__computed_val = None

    def __get__(self,obj,cls):
        ''' this is called on attribute read access
        it will be called only the first time. 
        '''
        if obj is None:
            return None
        if self.__computed_val is None:
            self.__computed_val = self.fget(obj)
        return self.__computed_val
#    def __set__(self, obj, val):
#        print('someone is trying to set this to {0}'.format(val))



class tryout(object):
    def __init__(self, someval):
        self._name = someval
    @lazy_property_desc
    def adata(self):
        print('someheavy liftin to compute the value')
        tempd = {}
        tempd[self._name]=self._name*3
        return tempd
        
    def getname(self):
        return self._name
    def setname(self, value):
        self._name = value
    name = property(getname, setname)

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

    
if __name__ == '__main__':
    y = Cnew()
    print('instance exists, data not yet loaded')
    print(y.data)
    print('let us do this again')
    print(y.data)