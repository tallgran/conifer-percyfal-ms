import copy


def get_envmodules(config, envmodules):
    retmodules = []
    if not isinstance(envmodules, list):
        envmodules = [envmodules]
    try:
        retmodules = copy.deepcopy(config["envmodules"]["__site__"])
    except KeyError:
        retmodules = []
    for mod in envmodules:
        try:
            retmodules.extend(config["envmodules"][mod])
        except KeyError:
            pass
    return retmodules
