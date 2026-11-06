"""Default COSMIC BSE flags, read from COSMIC's own settings file."""
import ast
import json
from importlib.resources import files


def default_bse_dict(**overrides):
    """Return COSMIC's documented default BSE flags, with optional overrides.

    Array-valued flags stored as strings in cosmic-settings.json are parsed.
    """
    path = files("cosmic.data").joinpath("cosmic-settings.json")
    with open(path) as f:
        settings = json.load(f)
    bse = {}
    for cat in settings:
        if cat["category"] != "bse":
            continue
        for s in cat["settings"]:
            dv = [o["name"] for o in s.get("options", []) if o.get("default")]
            val = dv[0] if dv else None
            if isinstance(val, str) and val.strip().startswith("["):
                val = eval(val, {"__builtins__": {}})  # e.g. "[2.0/21.0, ...]"
            bse[s["name"]] = val
    bse.update(overrides)
    return bse
