from collections import namedtuple

Node = namedtuple('Node', 'name start end left right')


def partition(features):
    """Partition `features` into a binary tree where `features` is an list
    of items (name string, start int, end int) ordered by 'start'. The
    tree is constructed of Nodes, each of which is a namedtuple with
    elements 'name start end left right'.

    * name - a string or other object identifying the node
    * start, end - the starting and ending coordinates of a range
    * left - a subtree in which node.start < start for all nodes or None
    * right - a subtree in which node.end > end for all nodes or None

    """

    if features:
        i = len(features) // 2
        left, (name, start, end), right = features[:i], features[i], features[i + 1:]
        return Node(name, start, end, partition(left), partition(right))


def assign(tree, val):
    """Assign `val` a name if it falls within an interval corresponding to
    a Node in `tree` or None otherwise.

    """

    if tree:
        if val < tree.start:
            return assign(tree.left, val)
        elif val > tree.end:
            return assign(tree.right, val)
        else:
            return tree.name


def check_features(features):
    """Check for elements of `features` with overlapping or mis-ordered
    ranges. In the case of overlap, print an informative error message
    and return names and positions of overlapping features.

    """

    features = features[:]
    errors = []
    for i in range(len(features) - 1):
        prev_name, prev_start, prev_end = features[i]
        name, start, end = features[i + 1]
        pair = ((prev_name, prev_start, prev_end), (name, start, end))

        if prev_end >= start:
            errors.append(pair)

    return errors


class RangeTree(object):
    """
    >>> features = [('foo', 1, 5), ('bar', 7, 10), ('baz', 12, 15)]
    >>> rt = RangeTree(features)
    >>> assert rt.assign(2) == 'foo'
    """

    def __init__(self, features):
        errors = check_features(features)
        if errors:
            raise ValueError('invalid features: {}'.format(errors))
        self.features = features
        self.tree = partition(features)

    def assign(self, val):
        return assign(self.tree, val)




# features = [('foo', 1, 5), ('bar', 7, 10), ('baz', 12, 15)]
# overlapping = [('foo', 1, 7), ('bar', 7, 10), ('baz', 12, 15)]
# misordered = [('bar', 7, 10), ('foo', 1, 7), ('baz', 12, 15)]

# import sys
# male_ages = [(-1, 0, 34), (0, 35, 39), (1, 40, 44), (2, 45, 49), (3, 50, 54),
#              (4, 55, 59), (5, 60, 64), (6, 65, 69), (7, 70, sys.maxint)]
# assert check_features(male_ages) == []

# female_ages = [(-9, 0, 34), (-4, 35, 39), (0, 40, 44), (3, 45, 49), (6, 50, 54),
#                (7, 55, 59), (7, 60, sys.maxint)]
# assert check_features(female_ages) == []


# male_tree = partition(male_ages)
# female_tree = partition(female_ages)
# assign(male_tree, 55)
# assign(female_tree, 55)


if __name__ == '__main__':
    print('testing...')
    import doctest
    doctest.testmod()
