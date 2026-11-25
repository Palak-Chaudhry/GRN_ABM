"""Example rules plugin: new KINDS of signal and behaviour, usable from rules.yaml once this file is listed
under `plugins:`. Copy this file to add your own; nothing else in the code needs to change.

  signal     distance_from_centre      um from the tissue centre (0, 0)
  signal     dying_neighbours          number of contacts that are dying (dying cells do not signal,
                                       so this uses distance, not the signalling adjacency)
  behaviour  module:<MODULE>           raises the transcription on-rate of every pro gene of a module and
                                       lowers every anti gene by the same amount (log on-rate change)
"""
import numpy as np

from grn_abm.rules import register_behaviour, register_signal


@register_signal("distance_from_centre")
def distance_from_centre(rules, sim, arg, cache):
    return np.linalg.norm(sim.tissue.pos, axis=1)


@register_signal("dying_neighbours")
def dying_neighbours(rules, sim, arg, cache):
    tis = sim.tissue
    R, f = tis.par("radius"), tis.par("max_adhesion_distance")
    pr = tis._pairs(2 * (R * f).max())
    out = np.zeros(len(tis.pos))
    if len(pr):
        i, j = pr[:, 0], pr[:, 1]
        close = np.linalg.norm(tis.pos[i] - tis.pos[j], axis=1) < f[i] * R[i] + f[j] * R[j]
        i, j = i[close], j[close]
        out = np.bincount(i, tis.dying[j], len(out)) + np.bincount(j, tis.dying[i], len(out))
    return out


@register_behaviour("module:", check=lambda rules, sim, arg: arg in rules.modules)
def module_programme(rules, sim, arg, value, out):
    if out.get("gene_mod") is None:
        out["gene_mod"] = np.zeros((len(value), len(sim.engine.genes)))
    m = rules.modules[arg]
    out["gene_mod"][:, m["pro"]] += value[:, None]
    out["gene_mod"][:, m["anti"]] -= value[:, None]
