# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of ngs_crumbs.
# ngs_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# ngs_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ngs_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
from crumbs.collectionz import OrderedSet, KeyedSet, RecentlyAddedCache

# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111


class TestCollections(unittest.TestCase):

    def test_ordered_set(self):
        in_list = [1, 2, 3, 4, 5, 8, 10]
        not_in_list = [6, 9, 11, 13]
        ordered_set = OrderedSet(in_list)
        for item in in_list:
            assert item in ordered_set
        assert ordered_set.check_add(7)
        assert ordered_set._items == [1, 2, 3, 4, 5, 7, 8, 10]
        assert not ordered_set.check_add(2)
        assert ordered_set._items == [1, 2, 3, 4, 5, 7, 8, 10]
        assert ordered_set.check_add(0)
        assert ordered_set._items == [0, 1, 2, 3, 4, 5, 7, 8, 10]
        for item in not_in_list:
            assert item not in ordered_set

    def test_unordered_set(self):
        in_set = [1, 2, 3, 4, 5, 8, 10]
        not_in_set = [6, 9, 11, 13]
        keyed_set = KeyedSet(in_set)
        for item in in_set:
            assert item in keyed_set
        assert keyed_set.check_add(7)
        assert keyed_set._items == set([1, 2, 3, 4, 5, 7, 8, 10])
        assert not keyed_set.check_add(2)
        assert keyed_set._items == set([1, 2, 3, 4, 5, 7, 8, 10])
        assert keyed_set.check_add(0)
        assert keyed_set._items == set([0, 1, 2, 3, 4, 5, 7, 8, 10])

        for item in not_in_set:
            assert item not in keyed_set

        #with key
        a = 'a'
        in_set = [(1, a), (2, a), (3, a), (4, a), (5, a),  (8, a), (10, a)]
        not_in_set = [(6, a), (9, a), (11, a), (13, a)]

        def key(item):
            return item[0]
        keyed_set = KeyedSet(in_set, key=key)
        for item in in_set:
            assert item in keyed_set
        assert keyed_set.check_add((7, a))
        assert keyed_set._items == set([1, 2, 3, 4, 5, 7, 8, 10])
        assert not keyed_set.check_add((2, a))
        assert keyed_set._items == set([1, 2, 3, 4, 5, 7, 8, 10])
        assert keyed_set.check_add((0, a))
        assert keyed_set._items == set([0, 1, 2, 3, 4, 5, 7, 8, 10])
        for item in not_in_set:
            assert item not in keyed_set


class RecentlyAddedCacheTest(unittest.TestCase):
    def test_recently_added_cache(self):
        cache = RecentlyAddedCache(3)
        cache[1] = 1
        assert cache[1]
        cache[2] = 2
        cache[3] = 3
        assert cache == {1: 1, 2: 2, 3: 3}
        cache[4] = 4
        assert cache == {2: 2, 3: 3, 4: 4}


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'TestCollections']
    unittest.main()
