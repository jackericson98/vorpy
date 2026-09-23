"""Restricted JSON/NumPy codec for the versioned network archive schema.

No dynamic imports, arbitrary constructors, pickle, or filesystem extraction.
Numeric table cells are packed by column, so millions of triangles do not
create millions of JSON values or ZIP members.
"""
import json
import math
import tempfile

import numpy as np
import pandas as pd


class ArchiveError(ValueError):
    pass


class Writer:
    def __init__(self, archive, references):
        self.archive = archive
        self.references = references
        self.array_count = 0
        self.tables = {}

    def json(self, path, value):
        with self.archive.open(path, 'w', force_zip64=True) as stream:
            for chunk in json.JSONEncoder(allow_nan=False, separators=(',', ':')).iterencode(value):
                stream.write(chunk.encode('utf-8'))

    def array(self, value):
        value = np.asarray(value)
        if value.dtype.kind not in 'biufc':
            raise ArchiveError(f'Unsupported numerical dtype: {value.dtype}')
        path = f'arrays/{self.array_count:06d}.npy'
        self.array_count += 1
        with self.archive.open(path, 'w', force_zip64=True) as stream:
            np.lib.format.write_array(stream, value, allow_pickle=False)
        return {'array': path, 'dtype': value.dtype.str, 'shape': list(value.shape)}

    def encode(self, value):
        if id(value) in self.references:
            return {'ref': self.references[id(value)]}
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, np.generic):
            return self.encode(value.item())
        if isinstance(value, float):
            return value if math.isfinite(value) else {'float': str(value)}
        if isinstance(value, np.ndarray):
            return self.array(value)
        if isinstance(value, pd.DataFrame):
            return self.table(value)
        if isinstance(value, dict):
            return {'map': [[self.encode(k), self.encode(v)] for k, v in value.items()]}
        if isinstance(value, (list, tuple, set)):
            kind = 'list' if isinstance(value, list) else 'tuple' if isinstance(value, tuple) else 'set'
            if len(value) > 32:
                try:
                    array = np.asarray(value)
                    if array.dtype.kind in 'biuf' and array.ndim <= 3:
                        return {'sequence': kind, 'values': self.array(array)}
                except (ValueError, TypeError):
                    pass
            return {kind: [self.encode(item) for item in value]}
        raise ArchiveError(f'Unsupported scientific value type: {type(value).__name__}')

    def column(self, values):
        # Scalar numerical columns, including ordinary nullable topology flags.
        array = np.asarray(values)
        if array.dtype.kind in 'biufc':
            return {'kind': 'numeric', 'values': self.array(array)}
        # Ragged numeric rows: coordinates, triangles, adjacency, tensors.
        dtype = None
        shapes = []
        sizes = []
        containers = []
        for value in values:
            if not isinstance(value, (list, tuple, np.ndarray)):
                break
            try:
                row = np.asarray(value)
            except (ValueError, TypeError):
                break
            if row.dtype.kind not in 'biufc' or row.ndim > 3:
                break
            if row.size:
                dtype = row.dtype if dtype is None else np.result_type(dtype, row.dtype)
            shapes.append(list(row.shape) + [-1] * (3 - row.ndim))
            sizes.append(row.size)
            containers.append(0 if isinstance(value, np.ndarray) else 1 if isinstance(value, list) else 2)
        else:
            if len(values):
                offsets = np.concatenate(([0], np.cumsum(sizes, dtype=np.int64)))
                # Disk-backed packing bounds peak RAM to one row plus the input.
                with tempfile.TemporaryDirectory(prefix='vorpy_archive_') as directory:
                    if offsets[-1]:
                        packed = np.lib.format.open_memmap(directory + '/packed.npy', mode='w+', dtype=dtype or np.int64,
                                                          shape=(int(offsets[-1]),))
                        for index, value in enumerate(values):
                            packed[offsets[index]:offsets[index + 1]] = np.asarray(value).reshape(-1)
                        data = self.array(packed)
                        del packed
                    else:
                        data = self.array(np.empty(0, dtype=dtype or np.int64))
                return {'kind': 'ragged', 'values': data, 'offsets': self.array(offsets),
                        'shapes': self.array(np.asarray(shapes, dtype=np.int64)),
                        'containers': self.array(np.asarray(containers, dtype=np.int8))}
        return {'kind': 'records', 'values': [self.encode(value) for value in values]}

    def table(self, table):
        if id(table) in self.tables:
            return {'table': self.tables[id(table)]}
        path = f'tables/{len(self.tables):04d}.json'
        self.tables[id(table)] = path
        if not table.index.is_unique:
            raise ArchiveError('Scientific table indexes must be unique')
        columns = {}
        for key in table.columns:
            if key in {'draw_points', 'draw_tris', 'tri_colors'}:
                continue  # disposable style-dependent drawing caches
            if not isinstance(key, str):
                raise ArchiveError('Scientific table column names must be strings')
            columns[key] = self.column(table[key].to_numpy())
        self.json(path, {'count': len(table), 'index': self.encode(table.index.to_numpy()), 'columns': columns})
        return {'table': path}


class Reader:
    def __init__(self, archive):
        self.archive = archive
        self.tables = {}
        self.references = {}

    def json(self, path):
        try:
            entry = self.archive.getinfo(path)
            if entry.file_size > 256 * 1024 * 1024:
                raise ArchiveError(f'JSON member too large: {path}')
            return json.loads(self.archive.read(path), parse_constant=lambda value: (_ for _ in ()).throw(ArchiveError('Nonstandard JSON constant')))
        except (KeyError, UnicodeError, json.JSONDecodeError) as error:
            raise ArchiveError(f'Invalid or missing JSON member {path}: {error}') from error

    def array(self, node):
        path = node['array']
        if not isinstance(path, str) or not path.startswith('arrays/') or not path.endswith('.npy'):
            raise ArchiveError('Invalid array member name')
        try:
            with self.archive.open(path) as stream:
                version = np.lib.format.read_magic(stream)
                shape, fortran, dtype = np.lib.format._read_array_header(stream, version)
                size = math.prod(shape) * dtype.itemsize
                if dtype.kind not in 'biufc' or size > self.archive.getinfo(path).file_size - stream.tell():
                    raise ArchiveError(f'Invalid array dtype/size: {path}')
            with self.archive.open(path) as stream:
                value = np.load(stream, allow_pickle=False)
        except (KeyError, ValueError, EOFError) as error:
            raise ArchiveError(f'Invalid array {path}: {error}') from error
        if list(value.shape) != node['shape'] or value.dtype.str != node['dtype']:
            raise ArchiveError(f'Array shape/dtype mismatch: {path}')
        return value

    def decode(self, node):
        if node is None or isinstance(node, (str, bool, int, float)):
            return node
        if not isinstance(node, dict):
            raise ArchiveError('Invalid encoded value')
        if 'ref' in node:
            if node['ref'] not in self.references:
                raise ArchiveError(f'Dangling object reference: {node["ref"]}')
            return self.references[node['ref']]
        if 'array' in node:
            return self.array(node)
        if 'table' in node:
            return self.table(node['table'])
        if 'map' in node:
            return {self.decode(key): self.decode(value) for key, value in node['map']}
        if 'float' in node and node['float'] in {'nan', 'inf', '-inf'}:
            return float(node['float'])
        for kind, constructor in (('list', list), ('tuple', tuple), ('set', set)):
            if kind in node:
                return constructor(self.decode(value) for value in node[kind])
            if node.get('sequence') == kind:
                return constructor(self.decode(node['values']).tolist())
        raise ArchiveError('Unknown encoded value')

    def column(self, node, count):
        kind = node['kind']
        if kind == 'numeric':
            result = self.array(node['values'])
            if result.shape != (count,):
                raise ArchiveError('Invalid scalar column shape')
            return result
        if kind == 'records':
            result = [self.decode(value) for value in node['values']]
        elif kind == 'ragged':
            data, offsets, shapes, containers = (self.array(node[key]) for key in ('values', 'offsets', 'shapes', 'containers'))
            if (data.ndim != 1 or offsets.shape != (count + 1,) or shapes.shape != (count, 3)
                    or containers.shape != (count,) or offsets.dtype.kind not in 'iu'
                    or shapes.dtype.kind not in 'iu' or offsets[0] != 0 or offsets[-1] != len(data)
                    or np.any(np.diff(offsets) < 0) or np.any(shapes < -1)):
                raise ArchiveError('Invalid packed column offsets/shapes')
            result = []
            for index in range(count):
                shape = tuple(int(v) for v in shapes[index] if v != -1)
                if math.prod(shape) != offsets[index + 1] - offsets[index]:
                    raise ArchiveError('Invalid packed row shape')
                value = data[offsets[index]:offsets[index + 1]].reshape(shape)
                container = containers[index]
                if container not in (0, 1, 2):
                    raise ArchiveError('Invalid packed row container')
                result.append(value if container == 0 else list(value) if container == 1 else tuple(value))
        else:
            raise ArchiveError('Unknown column encoding')
        if len(result) != count:
            raise ArchiveError('Column count mismatch')
        # Force object Series to keep arrays/lists as single DataFrame cells.
        return pd.Series(result, dtype=object).to_numpy()

    def table(self, path):
        if path in self.tables:
            return self.tables[path]
        if not path.startswith('tables/') or not path.endswith('.json'):
            raise ArchiveError('Invalid table member name')
        node = self.json(path)
        count = node['count']
        if not isinstance(count, int) or count < 0:
            raise ArchiveError('Invalid table count')
        index = self.decode(node['index'])
        if index.shape != (count,) or index.dtype.kind not in 'iu' or len(set(index)) != count:
            raise ArchiveError('Invalid table index')
        table = pd.DataFrame({key: self.column(value, count) for key, value in node['columns'].items()}, index=index)
        self.tables[path] = table
        return table
