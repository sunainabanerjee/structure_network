
__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['CacheData']


class CacheData:
    def __init__(self, **kwargs):
        self._cache = {}
        self._allowed_key = kwargs.get('allowed_keys', None)

    @property
    def restricted_access(self):
        return self._allowed_key is not None

    def add(self, ordered_keys, value):
        assert isinstance(ordered_keys, (list, tuple))
        curr_dict = self._cache
        depth = 0
        for x in ordered_keys[:-1]:
            if (self._allowed_key is not None) and (len(self._allowed_key) > depth):
                if x not in self._allowed_key[depth]:
                    raise RuntimeError("Error: Key: {} can not "
                                       "be inserted @ depth: {}".format(x, depth))
            if x not in curr_dict:
                curr_dict[x] = {}
            curr_dict = curr_dict[x]
            depth += 1
        if (ordered_keys[-1] in curr_dict) and \
                isinstance(curr_dict[ordered_keys[-1]], dict) and \
                isinstance(value, dict):
            curr_dict[ordered_keys[-1]].update(value)
        else:
            curr_dict[ordered_keys[-1]] = value
        return self

    def get(self, ordered_keys):
        assert isinstance(ordered_keys, (list, tuple))

        curr_dict = self._cache
        for x in ordered_keys:
            if x not in curr_dict:
                return {}
            curr_dict = curr_dict[x]
        return curr_dict

    def clear(self, ordered_keys=[]):
        assert isinstance(ordered_keys, (list, tuple))
        if len(ordered_keys) == 0:
            self._cache.clear()
        else:
            curr_dict = self._cache
            for x in ordered_keys[:-1]:
                if x not in curr_dict:
                    return self
                curr_dict = curr_dict[x]
            curr_dict.pop(ordered_keys[-1], {})
        return self

