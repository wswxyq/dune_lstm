""" Import Counter from collections module """
from collections import Counter

PDGLIST = {
    "p": 2212,
    "n": 2112,
    "pi+": 211,
    "pi-": -211,
    "pi0": 111,
}

FSIdic = {
    1: "one proton and nothing else",
    2: "one pion and nothing else",
    3: "at least one proton and any other particles",
    4: "at least one pion and any other particles but no protons",
    5: "no particles",
    6: "one neutron and any other particles (other than neutrons)",
    7: "multiple neutrons and any other particles",
    8: "only neutron/s",
    9: "no neutrons but there are other particles",
}


def get_fsi_state(input_s: str) -> list:
    # pylint: disable-msg=too-many-branches
    """
    takes a string and returns a list of the FSI labels.
    """
    _temp = input_s.split(",")[1:]
    _slist = []
    for x in _temp:
        if float(x) < 2000000000:
            _slist.append(int(float(x)))
    scouter = Counter(_slist)
    statuscode = []

    # one proton and nothing else
    if scouter == {PDGLIST["p"]: 1}:
        statuscode.append(1)

    # one pion and nothing else
    # if scouter == {PDGLIST["pi+"]: 1}:
    if _slist == [PDGLIST["pi+"]]:
        statuscode.append(2)
    # elif scouter == {PDGLIST["pi-"]: 1}:
    elif _slist == [PDGLIST["pi-"]]:
        statuscode.append(2)
    # elif scouter == {PDGLIST["pi0"]: 1}:
    elif _slist == [PDGLIST["pi0"]]:
        statuscode.append(2)

    # at least one proton and any other particles
    if (
        PDGLIST["p"] in scouter.keys()
        and scouter[PDGLIST["p"]] >= 1
        and len(scouter) > 1
    ):
        statuscode.append(3)

    # at least one pion and any other particles but no protons
    if (
        PDGLIST["pi+"] in scouter.keys()
        and scouter[PDGLIST["pi+"]] >= 1
        and len(scouter) > 1
        and PDGLIST["p"] not in scouter.keys()
    ):
        statuscode.append(4)
    elif (
        PDGLIST["pi-"] in scouter.keys()
        and scouter[PDGLIST["pi-"]] >= 1
        and len(scouter) > 1
        and PDGLIST["p"] not in scouter.keys()
    ):
        statuscode.append(4)
    elif (
        PDGLIST["pi0"] in scouter.keys()
        and scouter[PDGLIST["pi0"]] >= 1
        and len(scouter) > 1
        and PDGLIST["p"] not in scouter.keys()
    ):
        statuscode.append(4)

    # no particles
    if not scouter:
        statuscode.append(5)

    # one neutron and any other particles (other than neutrons)
    if (
        PDGLIST["n"] in scouter.keys()
        and scouter[PDGLIST["n"]] == 1
        and len(scouter) > 1
    ):
        statuscode.append(6)

    # multiple neutrons and any other particles
    if (
        PDGLIST["n"] in scouter.keys()
        and scouter[PDGLIST["n"]] > 1
        and len(scouter) > 1
    ):
        statuscode.append(7)

    # only neutron/s
    if PDGLIST["n"] in scouter.keys() and len(scouter) == 1:
        statuscode.append(8)

    # no neutrons but there are other particles
    if PDGLIST["n"] not in scouter.keys() and len(scouter) >= 1:
        statuscode.append(9)
    return statuscode


def get_fsi_bool(input_s: str, label: int) -> bool:
    """
    takes a string and returns a boolean if the label is present.
    label is an integer.
    """
    return label in get_fsi_state(input_s)


def get_str_max(input_s: str):
    """
    takes a string of array and returns the maximum element.
    """
    if isinstance(input_s, str):
        _temp = input_s.split(",")
        _slist = [float(x) for x in _temp]
        return max(_slist)
    return 0

def track_shower_01e8(input_s: str):
    """
    make a custom cut of track/shower energy > 0.1e8.
    """
    return get_str_max(input_s) > 0.1e8
    