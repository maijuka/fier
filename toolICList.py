class ICList(list):
    """
    A list-like object to keep track of the changes.
    
    To initialize from any iterable i, use iclist = ICList(i). By default, iclist.isChanged = False
    just after initialization; to override, use iclist = ICList(i, isChanged=True). Following
    operations set iclist.isChanged to true:
    iclist.isChanged = True
    iclist[i] = x
    del iclist[i]
    iclist[i:j] = aniterator
    del iclist[i:j]
    iclist += aniterator
    iclist *= i
    iclist.append(x)
    iclist.extend(alist)
    iclist.insert(i, x)
    iclist.remove(x)
    iclist.pop()
    iclist.sort()
    iclist.reverse()
    
    The following operations will return a new ICList with isChanged = True
    iclist + alist
    iclist * i
    """
    def __init__(self, data=[], isChanged = False):
        list.__init__(self, data)
        self._isChanged = isChanged

    @property
    def isChanged(self):
        """The property tracking wether the iclist has changed. Accepts only boolean values."""
        return self._isChanged

    @isChanged.setter
    def isChanged(self, value):
        if isinstance(value, bool):
            self._isChanged = value
        else:
            raise TypeError('The value of isChanged must be Boolean, is '+str(type(value)))

    # Container-type methods
    def __setitem__(self, key, value):
        list.__setitem__(self, key, value)
        self._isChanged = True

    def __delitem__(self, key):
        list.__delitem__(self, key)
        self._isChanged = True
        
    def __setslice__(self, i, j, sequence):
        list.__setslice__(self, i, j, sequence)
        self._isChanged = True

    def __delslice__(self, i, j):
        list.__delslice__(self, i, j)
        self._isChanged = True

    # "Numeric" operations for concatenation and repetition
    def __add__(self, other):
        return ICList(list.__add__(self, other), True)

    def __radd__(self, other):
        return ICList(list.__radd__(self, other), True)
    
    def __iadd__(self, other):
        list.__iadd__(self, other)
        self._isChanged = True
        return self

    def __mul__(self, other):
        return ICList(list.__mul__(self, other), True)

    def __rmull__(self, other):
        return ICList(list.__rmul__(self, other), True)

    def __imul__(self, other):
        list.__imul__(self, other)
        self._isChanged = True
        return self

    # Public methods for list
    def append(self, val):
        list.append(self, val)
        self._isChanged = True

    def extend(self, L):
        list.extend(self, L)
        self._isChanged = True

    def insert(self, i, x):
        list.insert(self, i, x)
        self._isChanged = True

    def remove(self, x):
        list.remove(self, x)
        self._isChanged = True

    def pop(self, i = None):
        self._isChanged = True
        if i is None:
            i = len(self)-1
        return list.pop(self, i)
        
    def sort(self, key=None, reverse=False):
        list.sort(self, key=key, reverse=reverse)
        self._isChanged = True

    def reverse(self):
        list.reverse(self)
        self._isChanged = True

    def reset(self):
        del self[:]
        self._isChanged = True

    def getElement(self, i):
        return self[i]
    def getIds(self):
        return range(len(self))

    def naturalOrder(self):
        return sorted(self)
    
    # Printing
    def __repr__(self):
        return 'ICList('+list.__repr__(self)+', isChanged = '+str(self.isChanged)+')'

    def __str__(self):
        return 'ICList: ('+", ".join(["%s" for v in self])+'), isChanged = ' + str(self.isChanged)
 
# def unitTest():
#     d = ICList()
#     assert(d.isChanged == False)
#     d.append(3)
#     assert(d.isChanged == True)
#     d.extend([4, 6])
#     assert(3 in d)
#     assert(5 not in d)
#     print(d)
    
# if __name__ == '__main__':
#     unitTest()
