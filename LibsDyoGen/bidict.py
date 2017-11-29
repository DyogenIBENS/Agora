#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Josh Bronson and contributors
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
Overview
--------

``bidict`` provides a bidirectional mapping data structure and related
functionality to naturally work with one-to-one relations in Python.

Unlike alternative implementations, ``bidict`` builds on top of the dict
API and supports the familiar ``__getitem__`` syntax. It also supports a
convenient slice syntax to express an inverse mapping::

    >>> element_by_symbol = bidict({'H': 'hydrogen'})
    >>> element_by_symbol['H']  # forward mapping works just like with dict
    'hydrogen'
    >>> element_by_symbol[:'hydrogen']  # use slice for the inverse mapping
    'H'

Syntax hacks ftw.

Motivation & More Examples
--------------------------

Python's built-in dict lets us associate unique keys with arbitrary values.
Because keys must be hashable, values can be looked up by key in constant time.
Different keys can map to the same value, but a single key cannot map to two
different values. For instance, {-1: 1, 0: 0, 1: 1} is a dict with
three unique keys and two unique values, because the keys -1 and 1 both map to
1. If you try to write its inverse {1: -1, 0: 0, 1: 1}, the dict that
comes out has only two mappings, one for key 1 and one for key 0; since key 1
is not allowed to map to both -1 and 1, one of these mappings is discarded.

Sometimes the relation we're modeling will only ever have a single key mapping
to a single value, as in the relation of chemical elements and their symbols.
This is called a one-to-one (or injective) mapping (see
https://en.wikipedia.org/wiki/Injective_mapping).

In this case we can be sure that the inverse mapping has the same number of
items as the forward mapping, and moreover that if key k maps to value v in the
forward mapping, value v maps to key k in the inverse. It would be useful then
to be able to look up keys by value in constant time in addition to being able
to look up values by key. With the additional constraint that values must be
hashable as well as keys, we can get constant-time forward and inverse lookups
via convenient syntax with ``bidict``.

Expanding on the previous example, anywhere the ``__getitem__`` syntax can be
used to reference a forward mapping, slice syntax can be used too::

    >>> element_by_symbol['H'] = 'Hydrogen'
    >>> element_by_symbol['H':]
    'Hydrogen'

Including setting and deleting items in either direction::

    >>> element_by_symbol['He':] = 'helium'
    >>> element_by_symbol[:'lithium'] = 'Li'
    >>> del element_by_symbol['H':]
    >>> del element_by_symbol[:'lithium']
    >>> element_by_symbol
    bidict({'He': 'helium'})

The rest of the ``MutableMapping`` interface is supported too::

    >>> 'C' in element_by_symbol
    False
    >>> element_by_symbol.get('C', 'carbon')
    'carbon'
    >>> element_by_symbol.pop('He')
    'helium'
    >>> element_by_symbol
    bidict({})
    >>> element_by_symbol.update(Hg='mercury')
    >>> element_by_symbol
    bidict({'Hg': 'mercury'})

You can also use the unary inverse operator ~ on a ``bidict`` to get the
inverse mapping in constant time::

    >>> ~element_by_symbol
    bidict({'mercury': 'Hg'})

Inverse can be composed with other ``MutableMapping`` APIs at no extra cost::

    >>> 'mercury' in ~element_by_symbol  # no more expensive than ``in element_by_symbol``
    True
    >>> (~element_by_symbol).pop('mercury')  # no more expensive than ``element_by_symbol.pop``
    'Hg'

See the ``bidict`` class for more examples.

The ``inverted`` iterator is also provided in the spirit of the built-in
function ``reversed``. Pass in a mapping to get the inverse mapping, an
iterable of pairs to get the pairs' inverses, or any object implementing an
``__inverted__`` method. See the ``inverted`` class for examples.

Note: It is intentional that the term "inverse" is used rather than "reverse".
Consider a collection of (k, v) pairs. Taking the reverse of the collection
can only be done if it is ordered (i.e. a sequence), and reverses the order of
the pairs in the collection, but each original (k, v) pair remains in the
resulting collection. By contrast, taking the inverse of such a collection does
not require an original ordering or say anything about the resulting ordering,
but rather just replaces every (k, v) pair with the inverse pair (v, k).

The ``namedbidict`` class factory can be used to create a bidirectional mapping
with customized names for the forward and the inverse mappings accessible via
attributes. See the ``namedbidict`` function for examples.

The built-in ``htmlentitydefs`` module provides an example of where ``bidict``
could be used in the Python standard library instead of maintaining the two
``name2codepoint`` and ``codepoint2name`` dictionaries separately.


Caveats
-------

Because ``bidict`` is a bidirectional dict, values as well as keys must be
hashable. Attempting to insert an unhashable value will result in an error::

    >>> anagrams_by_alphagram = bidict({'opt': ['opt', 'pot', 'top']})
    ... # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    TypeError:...unhashable...
    >>> bidict({'opt': ('opt', 'pot', 'top')})
    bidict({'opt': ('opt', 'pot', 'top')})

When instantiating or updating a ``bidict``, remember that mappings for
like values with differing keys will be silently dropped (just as the dict
literal {1: 'one', 1: 'uno'} silently drops a mapping), to maintain
bidirectionality::

    >>> nils = bidict({'zero': 0, 'zilch': 0, 'zip': 0})
    >>> len(nils)
    1
    >>> nils.update(nix=0, nada=0)
    >>> len(nils)
    1

When mapping the key of one existing mapping to the value of another (or
vice versa), the two mappings silently collapse into one::

    >>> b = bidict({1: 'one', 2: 'two'})
    >>> b[1] = 'two'
    >>> b
    bidict({1: 'two'})
    >>> b = bidict({1: 'one', 2: 'two'})
    >>> b[:'two'] = 1
    >>> b
    bidict({1: 'two'})

Links
-----

* Documentation: https://bidict.readthedocs.org

* Development: https://bitbucket.org/jab/bidict

Credits
-------

Thanks to Terry Reedy for the idea for the slice syntax.
Thanks to Raymond Hettinger for the idea for namedbidict and pointing out
various caveats.
Thanks to Francis Carr for the idea of storing the inverse bidict.

See the bidict module for further documentation.
------------------------------------------------
'''

from collections import MutableMapping
from functools import wraps
from re import compile
from sys import version_info

# compatibility
PY2 = version_info[0] == 2
if PY2:
    iteritems = lambda x: x.iteritems()
else:
    iteritems = lambda x: iter(x.items())


def iteritems_multi(*map_or_it, **kw):
    '''
    Generator yielding the mappings provided mirroring dict's ``__init__``.
    Useful so that an entire set of mappings as expressed by this signature
    does not have to be stored inside a new dict instance just to iterate
    over them::

        >>> list(iteritems_multi())
        []
        >>> list(iteritems_multi({1: 1}))
        [(1, 1)]
        >>> list(sorted(iteritems_multi({1: 1, 2: 4})))
        [(1, 1), (2, 4)]
        >>> list(iteritems_multi([(1, 1)]))
        [(1, 1)]
        >>> list(sorted(iteritems_multi([(1, 1), (2, 4)])))
        [(1, 1), (2, 4)]
        >>> list(sorted(iteritems_multi(one=1, two=2)))
        [('one', 1), ('two', 2)]
        >>> list(sorted(iteritems_multi([('one', 1)], two=2)))
        [('one', 1), ('two', 2)]
        >>> list(sorted(iteritems_multi({'one': 1}, two=2, three=3)))
        [('one', 1), ('three', 3), ('two', 2)]

    Of course, if both arguments are given and there are repeat entries, they
    will all be yielded, whereas with dict the last repeat clobbers earlier
    ones::

        >>> dict({'one': 1}, one=2)
        {'one': 2}
        >>> list(iteritems_multi({'one': 1}, one=2))
        [('one', 1), ('one', 2)]
    '''
    if map_or_it:
        l = len(map_or_it)
        if l != 1:
            raise TypeError('expected at most 1 argument, got %d' % l)
        map_or_it = map_or_it[0]
        if PY2 and hasattr(map_or_it, 'iteritems') or \
           (not PY2) and hasattr(map_or_it, 'items'):
               for (k, v) in iteritems(map_or_it):
                   yield (k, v)
        else:
            for (k, v) in map_or_it:
                yield (k, v)
    for (k, v) in iteritems(kw):
        yield (k, v)


class inverted(object):
    '''
    An iterator in the spirit of ``reversed``. Useful for inverting a mapping::

        >>> keys = (1, 2, 3)
        >>> vals = ('one', 'two', 'three')
        >>> fwd = dict(zip(keys, vals))
        >>> inv = dict(inverted(fwd))
        >>> inv == dict(zip(vals, keys))
        True

    Passing an iterable of pairs produces an iterable of the pairs' inverses::

        >>> seq = [(1, 'one'), (2, 'two'), (3, 'three')]
        >>> list(inverted(seq))
        [('one', 1), ('two', 2), ('three', 3)]

    Under the covers, ``inverted`` first tries to call ``__inverted__`` on the
    wrapped object and returns an iterator over the result if the call
    succeeds. If the call fails, ``inverted`` next tries to call ``items`` on
    the wrapped object, returning the inverses of the resulting pairs if the
    call succeeds. Finally, if the ``items`` call fails, ``inverted`` falls
    back on unpacking pairs from the wrapped object directly.

    This allows for passing an ``inverted`` object back into ``inverted`` to
    to get the original sequence of pairs back out::

        >>> seq == list(inverted(inverted(seq)))
        True

    Be careful with passing the inverse of a non-injective mapping into
    ``dict``; mappings for like values with differing keys will be dropped
    silently, just as ``{1: 'one', 1: 'uno'}`` silently drops a mapping::

        >>> squares = {-2: 4, -1: 1, 0: 0, 1: 1, 2: 4}
        >>> len(squares)
        5
        >>> len(dict(inverted(squares)))
        3
    '''

    def __init__(self, data):
        self._data = data

    def __iter__(self):
        data = self._data
        try:
            it = data.__inverted__
        except AttributeError:
            return self.next()
        else:
            return it()

    def __next__(self):
        data = self._data
        try:
            for k, v in iteritems(data):  # mapping?
                yield v, k
        except AttributeError:  # no, assume sequence
            for k, v in data:
                yield v, k

    next = __next__


class bidict(object):
    '''
    Bidirectional mapping implementing the ``MutableMapping`` interface, with
    additional facilities for retrieving inverse mappings. The API is a
    superset of the ``dict`` API (minus the ``fromkeys`` method, which doesn't
    make sense for a bidirectional mapping because keys *and* values must be
    unique).

    Examples::

        >>> keys = (1, 2, 3)
        >>> vals = ('one', 'two', 'three')
        >>> bi = bidict(zip(keys, vals))
        >>> bi == bidict({1: 'one', 2: 'two', 3: 'three'})
        True
        >>> bidict(inverted(bi)) == bidict(zip(vals, keys))
        True

    You can use standard subscripting syntax with a key to get or set a forward
    mapping::

        >>> bi[2]
        'two'
        >>> bi[2] = 'twain'
        >>> bi[2]
        'twain'
        >>> bi[4]
        Traceback (most recent call last):
            ...
        KeyError: 4

    Or use a slice with only a ``start``::

        >>> bi[2:]
        'twain'
        >>> bi[0:] = 'naught'
        >>> bi[0:]
        'naught'

    Use a slice with only a ``stop`` to get or set an inverse mapping::

        >>> bi[:'one']
        1
        >>> bi[:'aught'] = 1
        >>> bi[:'aught']
        1
        >>> bi[1]
        'aught'
        >>> bi[:'one']
        Traceback (most recent call last):
            ...
        KeyError: 'one'

    Deleting items from the bidict works the same way::

        >>> del bi[0]
        >>> del bi[2:]
        >>> del bi[:'three']
        >>> bi
        bidict({1: 'aught'})

    bidicts maintain references to their inverses via the ``inv`` property,
    which can also be used to access or modify them::

        >>> bi.inv
        bidict({'aught': 1})
        >>> bi.inv['aught']
        1
        >>> bi.inv[:1]
        'aught'
        >>> bi.inv[:1] = 'one'
        >>> bi.inv
        bidict({'one': 1})
        >>> bi
        bidict({1: 'one'})
        >>> bi.inv.inv is bi
        True
        >>> bi.inv.inv.inv is bi.inv
        True

    A ``bidict``’s inverse can also be accessed via the unary ~ operator, by
    analogy to the unary bitwise inverse operator::

        >>> ~bi
        bidict({'one': 1})
        >>> ~bi is bi.inv
        True

    Because ~ binds less tightly than brackets, parentheses are necessary for
    something like::

        >>> (~bi)['one']
        1

    bidicts work with ``inverted`` as expected::

        >>> biinv = bidict(inverted(bi))
        >>> biinv
        bidict({'one': 1})

    This of course creates a new object (equivalent but not identical)::

        >>> biinv == bi.inv
        True
        >>> biinv is bi.inv
        False

    This just demonstrated that ``__eq__`` has been implemented to work as
    expected. ``__neq__`` has too::

        >>> bi != biinv
        True

    bidicts should compare as expected to instances of other mapping
    types::

        >>> bi == dict([(1, 'one')])
        True

    Inverting the inverse should round-trip::

        >>> bi == bidict(inverted(inverted(bi)))
        True

    Use ``invert`` to invert the mapping in place::

        >>> bi.invert()
        >>> bi
        bidict({'one': 1})

    The rest of the ``MutableMapping`` interface is supported too::

        >>> bi.get('one')
        1
        >>> bi.setdefault('one', 2)
        1
        >>> bi.setdefault('two', 2)
        2
        >>> len(bi)  # calls __len__
        2
        >>> bi.pop('one')
        1
        >>> bi.popitem()
        ('two', 2)
        >>> bi.inv.setdefault(3, 'three')
        'three'
        >>> bi
        bidict({'three': 3})
        >>> [key for key in bi]  # calls __iter__, returns keys like dict
        ['three']
        >>> 'three' in bi  # calls __contains__
        True
        >>> list(bi.keys())
        ['three']
        >>> list(bi.values())
        [3]
        >>> bi.update([('four', 4)])
        >>> bi.update({'five': 5}, six=6, seven=7)
        >>> sorted(bi.items(), key=lambda x: x[1])
        [('three', 3), ('four', 4), ('five', 5), ('six', 6), ('seven', 7)]

    When instantiating or updating a ``bidict``, remember that mappings for
    like values with differing keys will be silently dropped (just as the
    literal ``{1: 'one', 1: 'uno'}`` silently drops a mapping), to maintain
    bidirectionality::

        >>> nils = bidict({'zero': 0, 'zilch': 0, 'zip': 0})
        >>> len(nils)
        1
        >>> nils.update(nix=0, nada=0)
        >>> len(nils)
        1

    Another caveat: when mapping the key of one existing mapping to the value
    of another (or vice versa), the two mappings collapse into one::

        >>> b = bidict({1: 'one', 2: 'two'})
        >>> b[1] = 'two'
        >>> b
        bidict({1: 'two'})
        >>> b = bidict({1: 'one', 2: 'two'})
        >>> b[:'two'] = 1
        >>> b
        bidict({1: 'two'})
    '''

    def __init__(self, *args, **kw):
        self._fwd = {}
        self._bwd = {}
        for (k, v) in iteritems_multi(*args, **kw):
            self[k] = v
        inv = object.__new__(self.__class__)
        inv._fwd = self._bwd
        inv._bwd = self._fwd
        inv._inv = self
        self._inv = inv

    def __invert__(self):
        '''
        Called when unary ~ operator is applied.
        '''
        return self._inv

    inv = property(__invert__)

    def __inverted__(self):
        return iteritems(self._bwd)

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self._fwd)
    __str__ = __repr__

    def __eq__(self, other):
        if isinstance(other, bidict):
            return self._fwd == other._fwd
        return self._fwd == other

    def __neq__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def _fwd_slice(slice):
        '''
        Raises ``TypeError`` if the given slice does not have either only
        its start or only its stop set to a non-None value.

        Returns True if only its start is not None and False if only its stop
        is not None.
        '''
        if slice.step is not None or \
            (not ((slice.start is None) ^ (slice.stop is None))):
            raise TypeError('Slice must specify only either start or stop')
        return slice.start is not None

    def __getitem__(self, keyorslice):
        if isinstance(keyorslice, slice):
            # forward lookup (by key): b[key:]
            if self._fwd_slice(keyorslice):
                return self._fwd[keyorslice.start]
            else:  # inverse lookup (by val): b[:val]
                return self._bwd[keyorslice.stop]
        else:  # keyorslice is a key: b[key]
            return self._fwd[keyorslice]

    def __del(self, key):
        val = self._fwd[key]
        del self._fwd[key]
        del self._bwd[val]

    def __delitem__(self, keyorslice):
        if isinstance(keyorslice, slice):
            # delete by key: del b[key:]
            if self._fwd_slice(keyorslice):
                self.__del(keyorslice.start)
            else:  # delete by value: del b[:val]
                self.__del(self._bwd[keyorslice.stop])
        else:  # keyorslice is a key: del b[key]
            self.__del(keyorslice)

    def __set(self, key, val):
        try:
            oldkey = self._bwd[val]
        except KeyError:
            pass
        else:
            del self._fwd[oldkey]
        try:
            oldval = self._fwd[key]
        except KeyError:
            pass
        else:
            del self._bwd[oldval]
        self._fwd[key] = val
        self._bwd[val] = key

    def __setitem__(self, keyorslice, keyorval):
        if isinstance(keyorslice, slice):
            # keyorslice.start is key, keyorval is val: b[key:] = val
            if self._fwd_slice(keyorslice):
                self.__set(keyorslice.start, keyorval)
            else:  # keyorval is key, keyorslice.stop is val: b[:val] = key
                self.__set(keyorval, keyorslice.stop)
        else:  # keyorslice is a key, keyorval is a val: b[key] = val
            self.__set(keyorslice, keyorval)

    def clear(self):
        self._fwd.clear()
        self._bwd.clear()

    def copy(self):
        return self.__class__(self._fwd)

    def invert(self):
        self._fwd, self._bwd = self._bwd, self._fwd
        self._inv._fwd, self._inv._bwd = self._inv._bwd, self._inv._fwd

    def pop(self, key, *args):
        val = self._fwd.pop(key, *args)
        del self._bwd[val]
        return val

    def popitem(self):
        if not self._fwd:
            raise KeyError
        key, val = self._fwd.popitem()
        del self._bwd[val]
        return key, val

    def setdefault(self, key, default=None):
        val = self._fwd.setdefault(key, default)
        self._bwd[val] = key
        return val

    def update(self, *args, **kw):
        for k, v in iteritems_multi(*args, **kw):
            self[k] = v

    def _proxied_to_fwd(method):
        '''
        Decorator which proxies calls to the given bidict method on to the
        self._fwd dict.
        '''
        @wraps(method, ('__name__', '__doc__'))
        def wrapper(self, *args, **kw):
            return method(self._fwd, *args, **kw)
        return wrapper

    _methodnames = [
        '__contains__', '__iter__', '__len__',
        'get', 'keys', 'items', 'values',
    ]
    if PY2:
        _methodnames.extend(('iteritems', 'iterkeys', 'itervalues'))
    for _methodname in _methodnames:
        locals()[_methodname] = _proxied_to_fwd(getattr(dict, _methodname))
    del _methodname, _methodnames

MutableMapping.register(bidict)


_LEGALNAMEPAT = '^[a-zA-Z][a-zA-Z0-9_]*$'
_LEGALNAMERE = compile(_LEGALNAMEPAT)

def empty_namedbidict(mapname, fwdname, invname):
    '''
    Create an empty instance of a custom bidict (namedbidict). This method is
    used to make 'namedbidict' instances picklable.
    '''
    return namedbidict(mapname, fwdname, invname)()

def namedbidict(mapname, fwdname, invname):
    '''
    Generate a custom bidict class in the spirit of ``namedtuple`` with
    custom attribute-based access to forward and inverse mappings::

        >>> ElementMap = namedbidict('ElementMap', 'symbol', 'element')
        >>> noble_gases = ElementMap({'He': 'helium'})
        >>> noble_gases.symbol['He']
        'helium'
        >>> noble_gases.element['helium']
        'He'
        >>> noble_gases.symbol['Ne'] = 'neon'
        >>> del noble_gases.element['helium']
        >>> noble_gases
        ElementMap({'Ne': 'neon'})

    Pass to ``bidict`` to get back a regular ``bidict``::

        >>> bidict(noble_gases)
        bidict({'Ne': 'neon'})

    Comparison works as expected::

        >>> noble_gases2 = ElementMap({'Ne': 'neon'})
        >>> noble_gases2 == noble_gases
        True
        >>> noble_gases2 == bidict(noble_gases)
        True
        >>> noble_gases2 == dict(noble_gases)
        True
        >>> noble_gases2['Rn'] = 'radon'
        >>> noble_gases2 == noble_gases
        False
        >>> noble_gases2 == bidict(noble_gases)
        False
        >>> noble_gases2 == dict(noble_gases)
        False
    '''
    for name in mapname, fwdname, invname:
        if _LEGALNAMERE.match(name) is None:
            raise ValueError('"%s" does not match pattern %s' %
                             (name, _LEGALNAMEPAT))
                             
    __dict__ = {fwdname: property(lambda self: self),
                invname: bidict.inv}

    custombidict = type(mapname, (bidict,), __dict__)

    # tell pickle how to handle this dynamically generated custom bidict
    custombidict.__reduce__ = lambda self: \
        (empty_namedbidict, (mapname, fwdname, invname), self.__dict__)

    return custombidict


if __name__ == '__main__':
    from doctest import testmod, ELLIPSIS
    testmod(optionflags=ELLIPSIS)
