import sys, json, numpy as np, warnings
warnings.filterwarnings("ignore")
from scipy import optimize
import importlib.util
spec = importlib.util.spec_from_file_location("rl", "ref_llbern.py"); rl = importlib.util.module_from_spec(spec); spec.loader.exec_module(rl)
what = sys.argv[1]
if what == "S2":
    eng = rl.make(2.0, 6.0, 0.8, 200.0, 2.0, 52.0)
    for n in (40, 48):
        print("S2", n, eng(n, [1.5])["probs"][0], flush=True)
else:
    lev, lo, hi, n = {"median": (0.5, 1.0, 2.0, 32), "upper": (0.975, 2.0, 60.0, 32)}[what]
    eng = rl.make(1.0, 6.0, 0.8, 200.0, 1.0, 104.0)
    f = lambda lq: eng(n, [np.exp(lq)])["probs"][0] - (1 - lev)
    q = float(np.exp(optimize.brentq(f, np.log(lo), np.log(hi), xtol=1e-5)))
    print(what, n, q, flush=True)
