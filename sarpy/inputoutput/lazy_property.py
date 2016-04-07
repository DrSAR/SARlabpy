# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:06:35 2013

@author: stefan
"""

class lazy_property(object):
    '''
    Meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    
    An example illustrates the use. First define a class that contains
    a lazy_property. That call it and observe how the calculation is 
    only performed once.
    
    >>> class Test(object):
    ...     @lazy_property
    ...     def results(self):
    ...         calcs = 42 # do a lot of calculation here
    ...         print('phew, this took a lot of effort - glad I do this so rarely...')
    ...         return calcs
    >>> A=Test()
    >>> A.__dict__
    {}
    >>> A.results
    phew, this took a lot of effort - glad I do this so rarely...
    42
    >>> A.results  # on a second call the calculation is not needed !
    42
    >>> A.__dict__  # see how results is no a proper attribute
    {'results': 42}

    '''

    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__

    def __get__(self,obj,cls):
        ''' this is called on attribute read access
        it will be called only the first time. 
        '''
        if obj is None:
            return None
        value = self.fget(obj)
        # the following overwrites the dict in the calling object
        # thereby erasing any access to this function (which)
        # calls the expensive calculation once.
        setattr(obj,self.func_name,value)
        return value
