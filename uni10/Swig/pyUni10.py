# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _pyUni10
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class PySwigIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PySwigIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PySwigIterator, name)
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _pyUni10.delete_PySwigIterator
    __del__ = lambda self : None;
    def value(*args): return _pyUni10.PySwigIterator_value(*args)
    def incr(*args): return _pyUni10.PySwigIterator_incr(*args)
    def decr(*args): return _pyUni10.PySwigIterator_decr(*args)
    def distance(*args): return _pyUni10.PySwigIterator_distance(*args)
    def equal(*args): return _pyUni10.PySwigIterator_equal(*args)
    def copy(*args): return _pyUni10.PySwigIterator_copy(*args)
    def next(*args): return _pyUni10.PySwigIterator_next(*args)
    def previous(*args): return _pyUni10.PySwigIterator_previous(*args)
    def advance(*args): return _pyUni10.PySwigIterator_advance(*args)
    def __eq__(*args): return _pyUni10.PySwigIterator___eq__(*args)
    def __ne__(*args): return _pyUni10.PySwigIterator___ne__(*args)
    def __iadd__(*args): return _pyUni10.PySwigIterator___iadd__(*args)
    def __isub__(*args): return _pyUni10.PySwigIterator___isub__(*args)
    def __add__(*args): return _pyUni10.PySwigIterator___add__(*args)
    def __sub__(*args): return _pyUni10.PySwigIterator___sub__(*args)
    def __iter__(self): return self
PySwigIterator_swigregister = _pyUni10.PySwigIterator_swigregister
PySwigIterator_swigregister(PySwigIterator)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pyUni10.IntVector_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pyUni10.IntVector___nonzero__(*args)
    def __len__(*args): return _pyUni10.IntVector___len__(*args)
    def pop(*args): return _pyUni10.IntVector_pop(*args)
    def __getslice__(*args): return _pyUni10.IntVector___getslice__(*args)
    def __setslice__(*args): return _pyUni10.IntVector___setslice__(*args)
    def __delslice__(*args): return _pyUni10.IntVector___delslice__(*args)
    def __delitem__(*args): return _pyUni10.IntVector___delitem__(*args)
    def __getitem__(*args): return _pyUni10.IntVector___getitem__(*args)
    def __setitem__(*args): return _pyUni10.IntVector___setitem__(*args)
    def append(*args): return _pyUni10.IntVector_append(*args)
    def empty(*args): return _pyUni10.IntVector_empty(*args)
    def size(*args): return _pyUni10.IntVector_size(*args)
    def clear(*args): return _pyUni10.IntVector_clear(*args)
    def swap(*args): return _pyUni10.IntVector_swap(*args)
    def get_allocator(*args): return _pyUni10.IntVector_get_allocator(*args)
    def begin(*args): return _pyUni10.IntVector_begin(*args)
    def end(*args): return _pyUni10.IntVector_end(*args)
    def rbegin(*args): return _pyUni10.IntVector_rbegin(*args)
    def rend(*args): return _pyUni10.IntVector_rend(*args)
    def pop_back(*args): return _pyUni10.IntVector_pop_back(*args)
    def erase(*args): return _pyUni10.IntVector_erase(*args)
    def __init__(self, *args): 
        this = _pyUni10.new_IntVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pyUni10.IntVector_push_back(*args)
    def front(*args): return _pyUni10.IntVector_front(*args)
    def back(*args): return _pyUni10.IntVector_back(*args)
    def assign(*args): return _pyUni10.IntVector_assign(*args)
    def resize(*args): return _pyUni10.IntVector_resize(*args)
    def insert(*args): return _pyUni10.IntVector_insert(*args)
    def reserve(*args): return _pyUni10.IntVector_reserve(*args)
    def capacity(*args): return _pyUni10.IntVector_capacity(*args)
    __swig_destroy__ = _pyUni10.delete_IntVector
    __del__ = lambda self : None;
IntVector_swigregister = _pyUni10.IntVector_swigregister
IntVector_swigregister(IntVector)

class Qnum_arr(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Qnum_arr, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Qnum_arr, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pyUni10.Qnum_arr_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pyUni10.Qnum_arr___nonzero__(*args)
    def __len__(*args): return _pyUni10.Qnum_arr___len__(*args)
    def pop(*args): return _pyUni10.Qnum_arr_pop(*args)
    def __getslice__(*args): return _pyUni10.Qnum_arr___getslice__(*args)
    def __setslice__(*args): return _pyUni10.Qnum_arr___setslice__(*args)
    def __delslice__(*args): return _pyUni10.Qnum_arr___delslice__(*args)
    def __delitem__(*args): return _pyUni10.Qnum_arr___delitem__(*args)
    def __getitem__(*args): return _pyUni10.Qnum_arr___getitem__(*args)
    def __setitem__(*args): return _pyUni10.Qnum_arr___setitem__(*args)
    def append(*args): return _pyUni10.Qnum_arr_append(*args)
    def empty(*args): return _pyUni10.Qnum_arr_empty(*args)
    def size(*args): return _pyUni10.Qnum_arr_size(*args)
    def clear(*args): return _pyUni10.Qnum_arr_clear(*args)
    def swap(*args): return _pyUni10.Qnum_arr_swap(*args)
    def get_allocator(*args): return _pyUni10.Qnum_arr_get_allocator(*args)
    def begin(*args): return _pyUni10.Qnum_arr_begin(*args)
    def end(*args): return _pyUni10.Qnum_arr_end(*args)
    def rbegin(*args): return _pyUni10.Qnum_arr_rbegin(*args)
    def rend(*args): return _pyUni10.Qnum_arr_rend(*args)
    def pop_back(*args): return _pyUni10.Qnum_arr_pop_back(*args)
    def erase(*args): return _pyUni10.Qnum_arr_erase(*args)
    def __init__(self, *args): 
        this = _pyUni10.new_Qnum_arr(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pyUni10.Qnum_arr_push_back(*args)
    def front(*args): return _pyUni10.Qnum_arr_front(*args)
    def back(*args): return _pyUni10.Qnum_arr_back(*args)
    def assign(*args): return _pyUni10.Qnum_arr_assign(*args)
    def resize(*args): return _pyUni10.Qnum_arr_resize(*args)
    def insert(*args): return _pyUni10.Qnum_arr_insert(*args)
    def reserve(*args): return _pyUni10.Qnum_arr_reserve(*args)
    def capacity(*args): return _pyUni10.Qnum_arr_capacity(*args)
    __swig_destroy__ = _pyUni10.delete_Qnum_arr
    __del__ = lambda self : None;
Qnum_arr_swigregister = _pyUni10.Qnum_arr_swigregister
Qnum_arr_swigregister(Qnum_arr)

class Bond_arr(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Bond_arr, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Bond_arr, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pyUni10.Bond_arr_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pyUni10.Bond_arr___nonzero__(*args)
    def __len__(*args): return _pyUni10.Bond_arr___len__(*args)
    def pop(*args): return _pyUni10.Bond_arr_pop(*args)
    def __getslice__(*args): return _pyUni10.Bond_arr___getslice__(*args)
    def __setslice__(*args): return _pyUni10.Bond_arr___setslice__(*args)
    def __delslice__(*args): return _pyUni10.Bond_arr___delslice__(*args)
    def __delitem__(*args): return _pyUni10.Bond_arr___delitem__(*args)
    def __getitem__(*args): return _pyUni10.Bond_arr___getitem__(*args)
    def __setitem__(*args): return _pyUni10.Bond_arr___setitem__(*args)
    def append(*args): return _pyUni10.Bond_arr_append(*args)
    def empty(*args): return _pyUni10.Bond_arr_empty(*args)
    def size(*args): return _pyUni10.Bond_arr_size(*args)
    def clear(*args): return _pyUni10.Bond_arr_clear(*args)
    def swap(*args): return _pyUni10.Bond_arr_swap(*args)
    def get_allocator(*args): return _pyUni10.Bond_arr_get_allocator(*args)
    def begin(*args): return _pyUni10.Bond_arr_begin(*args)
    def end(*args): return _pyUni10.Bond_arr_end(*args)
    def rbegin(*args): return _pyUni10.Bond_arr_rbegin(*args)
    def rend(*args): return _pyUni10.Bond_arr_rend(*args)
    def pop_back(*args): return _pyUni10.Bond_arr_pop_back(*args)
    def erase(*args): return _pyUni10.Bond_arr_erase(*args)
    def __init__(self, *args): 
        this = _pyUni10.new_Bond_arr(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pyUni10.Bond_arr_push_back(*args)
    def front(*args): return _pyUni10.Bond_arr_front(*args)
    def back(*args): return _pyUni10.Bond_arr_back(*args)
    def assign(*args): return _pyUni10.Bond_arr_assign(*args)
    def resize(*args): return _pyUni10.Bond_arr_resize(*args)
    def insert(*args): return _pyUni10.Bond_arr_insert(*args)
    def reserve(*args): return _pyUni10.Bond_arr_reserve(*args)
    def capacity(*args): return _pyUni10.Bond_arr_capacity(*args)
    __swig_destroy__ = _pyUni10.delete_Bond_arr
    __del__ = lambda self : None;
Bond_arr_swigregister = _pyUni10.Bond_arr_swigregister
Bond_arr_swigregister(Bond_arr)

PRT_EVEN = _pyUni10.PRT_EVEN
PRT_ODD = _pyUni10.PRT_ODD
PRTF_EVEN = _pyUni10.PRTF_EVEN
PRTF_ODD = _pyUni10.PRTF_ODD
class Qnum(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Qnum, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Qnum, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _pyUni10.new_Qnum(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _pyUni10.delete_Qnum
    __del__ = lambda self : None;
    def U1(*args): return _pyUni10.Qnum_U1(*args)
    def prt(*args): return _pyUni10.Qnum_prt(*args)
    def prtF(*args): return _pyUni10.Qnum_prtF(*args)
    def assign(*args): return _pyUni10.Qnum_assign(*args)
    __swig_getmethods__["isFermionic"] = lambda x: _pyUni10.Qnum_isFermionic
    if _newclass:isFermionic = staticmethod(_pyUni10.Qnum_isFermionic)
    def __eq__(*args): return _pyUni10.Qnum___eq__(*args)
    def __mul__(*args): return _pyUni10.Qnum___mul__(*args)
    def __neg__(*args): return _pyUni10.Qnum___neg__(*args)
    def __copy__(*args): return _pyUni10.Qnum___copy__(*args)
    def cp(*args): return _pyUni10.Qnum_cp(*args)
    def __str__(*args): return _pyUni10.Qnum___str__(*args)
    def assignF(*args): return _pyUni10.Qnum_assignF(*args)
    U1_UPB = _pyUni10.Qnum_U1_UPB
    U1_LOB = _pyUni10.Qnum_U1_LOB
Qnum_swigregister = _pyUni10.Qnum_swigregister
Qnum_swigregister(Qnum)
QnumF = _pyUni10.QnumF
Qnum_isFermionic = _pyUni10.Qnum_isFermionic

BD_IN = _pyUni10.BD_IN
BD_OUT = _pyUni10.BD_OUT
class Bond(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Bond, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Bond, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _pyUni10.new_Bond(*args)
        try: self.this.append(this)
        except: self.this = this
    def assign(*args): return _pyUni10.Bond_assign(*args)
    def type(*args): return _pyUni10.Bond_type(*args)
    def dim(*args): return _pyUni10.Bond_dim(*args)
    def change(*args): return _pyUni10.Bond_change(*args)
    def combine(*args): return _pyUni10.Bond_combine(*args)
    def degeneracy(*args): return _pyUni10.Bond_degeneracy(*args)
    def Qlist(*args): return _pyUni10.Bond_Qlist(*args)
    __swig_destroy__ = _pyUni10.delete_Bond
    __del__ = lambda self : None;
Bond_swigregister = _pyUni10.Bond_swigregister
Bond_swigregister(Bond)


combine = _pyUni10.combine

