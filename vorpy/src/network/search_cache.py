"""Bounded caches for recomputable search data, never solved vertex state."""
from collections import OrderedDict
import sys

import numpy as np


CACHE_BYTES = 32 * 1024 * 1024  # Per bucket; five buckets use at most ~160 MiB.
CACHE_ENTRIES = 8192


def retained_size(value, seen=None):
    """Conservatively count containers and native NumPy buffers per entry."""
    if seen is None:
        seen = set()
    identity = id(value)
    if identity in seen:
        return 0
    seen.add(identity)
    size = sys.getsizeof(value)
    if isinstance(value, np.ndarray):
        # Owning arrays include their buffer in getsizeof; views can keep a
        # much larger base allocation alive, so include that allocation too.
        if value.base is not None:
            size += retained_size(value.base, seen)
    elif isinstance(value, dict):
        size += sum(retained_size(k, seen) + retained_size(v, seen)
                    for k, v in value.items())
    elif isinstance(value, (list, tuple)):
        size += sum(retained_size(item, seen) for item in value)
    return size


class BoundedCache:
    """Evict least recently used entries; oversized values are not retained."""
    def __init__(self, max_bytes=CACHE_BYTES, max_entries=CACHE_ENTRIES, size_of=retained_size):
        self.max_bytes = max_bytes
        self.max_entries = max_entries
        self.bytes = 0
        self.evictions = 0
        self.skipped = 0
        self.hits = 0
        self.misses = 0
        self.size_of = size_of
        self._values = OrderedDict()

    def __len__(self):
        return len(self._values)

    def __contains__(self, key):
        return key in self._values

    def __getitem__(self, key):
        try:
            value, _ = self._values[key]
        except KeyError:
            self.misses += 1
            raise
        self._values.move_to_end(key)
        self.hits += 1
        return value

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __setitem__(self, key, value):
        size = self.size_of((key, value)) + 128  # Mapping/bookkeeping overhead.
        old = self._values.pop(key, None)
        if old is not None:
            self.bytes -= old[1]
        if size > self.max_bytes or self.max_entries <= 0:
            self.skipped += 1
            return
        while self._values and (self.bytes + size > self.max_bytes
                                or len(self) >= self.max_entries):
            _, (_, removed_size) = self._values.popitem(last=False)
            self.bytes -= removed_size
            self.evictions += 1
        self._values[key] = (value, size)
        self.bytes += size


def _index_list_size(indices):
    # Ball/lookup indices fit in platform integers. Count each separately
    # (including duplicates) without walking thousands of scalar objects.
    return sys.getsizeof(indices) + (len(indices) * 36 if indices is not None else 0)


def _candidate_size(entry):
    key, value = entry
    return sys.getsizeof(entry) + retained_size(key) + _index_list_size(value)


def _neighborhood_size(entry):
    key, value = entry
    size = sys.getsizeof(entry) + retained_size(key) + sys.getsizeof(value)
    for item in value:
        if isinstance(item, np.ndarray):
            size += retained_size(item)
        elif isinstance(item, dict):
            size += sys.getsizeof(item) + len(item) * 72
        else:
            size += _index_list_size(item)
    return size


def cache_bucket(cache, name):
    if name not in cache:
        size_of = (_candidate_size if name == 'candidates' else
                   _neighborhood_size if name in ('surrounding', 'aw_verification') else
                   retained_size)
        cache[name] = BoundedCache(size_of=size_of)
    return cache[name]


def cache_summary(cache):
    buckets = [bucket for bucket in cache.values() if isinstance(bucket, BoundedCache)]
    return (f'cache={sum(bucket.bytes for bucket in buckets) / 1024 ** 2:.1f} MiB'
            f'; cache entries={sum(len(bucket) for bucket in buckets):,}'
            f'; evictions={sum(bucket.evictions for bucket in buckets):,}'
            f'; oversized skipped={sum(bucket.skipped for bucket in buckets):,}'
            + '; cache hits/misses=' + ','.join(
                f'{name}:{bucket.hits}/{bucket.misses}' for name, bucket in cache.items()
                if isinstance(bucket, BoundedCache)))
