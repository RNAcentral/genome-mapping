import inspect

import attr


class NonUniqueName(Exception):
    pass


def get_children(module, parent, names, ignore=set()):
    found = []
    for cls in children_of(module, parent, ignore=ignore):
        if cls.name in names:
            found.append(cls)
    if found:
        if len(found) == len(names):
            return found
        raise ValueError("Could not find all children %s" % names)
    raise ValueError("Unknown children %s" % names)


def get_child(module, parent, name, ignore=set()):
    return get_children(module, parent, set([name]))[0]


def children_of(module, parent, ignore=set()):
    ignore = ignore or set()

    def is_child(member):
        return inspect.isclass(member) and issubclass(member, parent) and \
            member != parent and member not in ignore

    pairs = inspect.getmembers(module, is_child)
    classes = {p[1] for p in pairs}
    names = set()
    for cls in classes:
        if cls.name in names:
            raise NonUniqueName("Non unique name %s" % cls.name)
        names.add(cls.name)
    return classes


def names_of_children(module, parent, ignore=set()):
    return {cls.name for cls in children_of(module, parent, ignore=set())}


def properities_of(klass):
    return {a.name for a in attr.fields(klass)}
