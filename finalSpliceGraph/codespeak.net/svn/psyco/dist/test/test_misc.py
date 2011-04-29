import py, sys, psyco
import math


def test_index():
    if sys.version_info < (2, 5):
        py.test.skip("for Python 2.5")

    class X(object):
        def __index__(self):
            return -3

    def f(x):
        return (12, 23, 34, 45)[x]

    res = psyco.proxy(f)(X())
    assert res == 23


def test_index_slice():
    if sys.version_info < (2, 5):
        py.test.skip("for Python 2.5")

    class Observer(object):
        def __getslice__(self, i, j):
            return "slice", i, j
        def __getitem__(self, x):
            raise Exception("in __getitem__")
    class X(object):
        def __index__(self):
            return 5
    class Y(object):
        def __index__(self):
            return -2

    def f(o, x, y):
        return o[x:y]

    res = psyco.proxy(f)(Observer(), X(), Y())
    assert res == ("slice", 5, -2)

def test_index_repeat():
    if sys.version_info < (2, 5):
        py.test.skip("for Python 2.5")

    class X(object):
        def __index__(self):
            return 3

    def f(x):
        return (12, 23) * x

    res = psyco.proxy(f)(X())
    assert res == (12, 23, 12, 23, 12, 23)

def test_slice_overflow():
    class X(object):
        def __len__(self):
            return sys.maxint
        def __getslice__(self, i, j):
            return i, j
    def f(X):
        return X()[-(2**100):2**100]

    res = psyco.proxy(f)(X)
    assert res == f(X)

def test_id():
    def f(x):
        return id(x)
    obj = []
    res = psyco.proxy(f)(obj)
    assert res == id(obj)

def test_math_module():
    import math
    global _operation
    for name in ["acos", "asin", "atan", "ceil", "cos", "cosh", "exp",
                 "fabs", "floor", "sin", "sinh", "sqrt", "tan", "tanh"]:
        _operation = getattr(math, name)

        def tester(a):
            return _operation(a)
        res = psyco.proxy(tester)(0.17)
        assert abs(res - _operation(0.17)) < 1e6

        def tester(a):
            return _operation(*a)
        res = psyco.proxy(tester)((0.71,))
        assert abs(res - _operation(0.71)) < 1e6

    for name in ["atan2", "fmod", "hypot", "pow"]:
        _operation = getattr(math, name)

        def tester(a, b):
            return _operation(a, b)
        res = psyco.proxy(tester)(0.17, 0.081)
        assert abs(res - _operation(0.17, 0.081)) < 1e6

        def tester(a):
            return _operation(*a)
        res = psyco.proxy(tester)((0.71, 0.2643))
        assert abs(res - _operation(0.71, 0.2643)) < 1e6

def test_circular_list():
    def f():
        l = [None]
        l[0] = l
        return l
    lst = psyco.proxy(f)()
    assert type(lst) is list
    assert len(lst) == 1
    assert lst[0] is lst

    def g():
        a = [None]
        b = [None]
        a[0] = b
        b[0] = a
        return b
    lst = psyco.proxy(g)()
    assert type(lst) is list
    assert len(lst) == 1
    assert lst[0] is not lst
    assert type(lst[0]) is list
    assert len(lst[0]) == 1
    assert lst[0][0] is lst

def test_big_dict_constructor():
    def f(x):
        return {1: 2, 3: 4, x: 6, 7: 8, 9: x,
                'a': 'a', 'b': 'c', 'd': 'e', 'f': 'g'}
    assert psyco.proxy(f)(5) == f(5)
    assert psyco.proxy(f)('b') == f('b')

def test_nonfloat_in_math():
    class X(object):
        def __float__(self):
            return -123.45
    def f(x):
        return math.ceil(x), math.floor(x)
    assert psyco.proxy(f)(X()) == (-123.0, -124.0)

def test_power_error_case():
    def f(x, y):
        return x ** y
    py.test.raises(ValueError, f, -4.5, 0.9)
    py.test.raises(ValueError, psyco.proxy(f), -4.5, 0.9)
