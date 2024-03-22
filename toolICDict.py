
def first_value(d):
    try:
        k, *_ = d
        return d[k]
    except ValueError:
        return None

class ICDict(dict):
    """
    A dict like object that keeps track of the changes.

    Initialization and use like normal dict, except that doesn't support initialization from keywords
    (as in 'dict(one=1, two=2)'). Has property isChanged, which is initially set to False, unless other
    value is given to __init__() as second argument.
    The following actions set d.isChanged to True for ICDict d:
    d[key] = value
    del d[key]
    d.clear()
    d.pop(key)
    d.popitem()
    d.setdefault(key)
    d.update()

    Method d.copy() will return a new ICDict with isChanged property set to that of d.

    Method d.fromkeys(seq) will return a new ICDict with isChanged = False.
    """

    def __init__(self, iterable=[], isChanged = False):
        dict.__init__(self, iterable)
        self._isChanged = bool(isChanged)

    @property
    def isChanged(self):
        "The property for tracking if the ICDict has changed. Accepts only Boolean values"
        return self._isChanged

    @isChanged.setter
    def isChanged(self, value):
        if isinstance(value, bool):
            self._isChanged = value
        else:
            raise TypeError('The value of isChanged must be Boolean, is '+str(type(value)))

    def __setitem__(self, key, value):
        try:
            dict.__setitem__(self, key, value)
        except:
            raise
        else:
            self._isChanged = True

    def __delitem__(self, key):
        try:
            dict.__delitem__(self, key)
        except:
            raise
        else:
            self._isChanged = True

    def copy(self):
        new = ICDict(self.data)
        new.isChanged = self._isChanged
        return new

    def reset(self):
        self.clear()
        self._isChanged = True

    def getElement(self, i):
        return self.get(i, None)
    def getIds(self):
        return self.keys()


    @classmethod
    def fromkeys(cls, seq, value=None):
        return ICDict(dict.fromkeys(seq, value), False)        

    def __repr__(self):
        return 'ICDict('+ dict.__repr__(self) + ', isChanged=' + str(self._isChanged) + ')'

    def __str__(self):
        return 'ICDict: {'+ ", ".join(["%s: %s" % (k,v) for k,v in self.items()]) +'}, isChanged = ' + str(self._isChanged)

# def unitTest():
#     t = ICDict.fromkeys(["x", "y", "z"], 4)
#     print(t)
#     d = ICDict()
#     assert(d.isChanged == False)
#     d[3] = "three"
#     assert(d.isChanged == True)
#     assert(3 in d)
#     assert(4 not in d)
#     print(d)
    
# if __name__ == '__main__':
#     unitTest()
