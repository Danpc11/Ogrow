"""Single source of parameters and seeds."""
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]

def load(path=None):
    p = Path(path) if path else ROOT / "config" / "params.yaml"
    with open(p) as f:
        return yaml.safe_load(f)

CFG = load()

def seeds(name):
    s = CFG["seeds"][name]
    return range(*s) if isinstance(s, list) else [s]
