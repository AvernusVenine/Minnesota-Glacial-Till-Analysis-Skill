"""M1 weights, M5-M7: the unit model, the probability of a sample under a unit, and typicality.

Grains: every view (age, lithology, colour, individual grains) is its own Dirichlet-multinomial
with two spreads per unit, between boreholes (lateral, A_L) and within a borehole (A_W).
Larger A = tighter. The overall number is a weighted sum of the views' log probabilities;
each view's weight discounts the counts it shares with other views (M1).
Texture: bivariate normal on sand % and clay % with a between-borehole covariance (S_L)
and a within-borehole covariance (S_W).
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import chi2

LOG_A_MIN, LOG_A_MAX = np.log(1.0), np.log(5000.0)
DEFAULT_AW, DEFAULT_AL = 60.0, 20.0
PI_FLOOR = 1e-6
TEX_FLOOR = 0.25      # %^2: texture is recorded to the nearest percent


def dm_logpmf(x, alpha):
    """Log probability of counts x under a Dirichlet-multinomial with parameters alpha."""
    x = np.asarray(x, float)
    alpha = np.asarray(alpha, float)
    n, A = x.sum(-1), alpha.sum(-1)
    return (gammaln(n + 1) - gammaln(x + 1).sum(-1) + gammaln(A) - gammaln(n + A)
            + (gammaln(x + alpha) - gammaln(alpha)).sum(-1))


def neff(n, AW):
    """Grains a sample is worth for locating its borehole's composition (design effect)."""
    return n * (AW + 1.0) / (AW + n)


def dm_pred_alpha(pi, AW, AL, S, E):
    """M6: parameters for a new sample, given same-borehole evidence S (effective counts, total E)."""
    a = AL * pi + S
    atot = AL + E
    r_post, r_w = 1.0 / (atot + 1.0), 1.0 / (AW + 1.0)
    r = r_post + r_w - r_post * r_w
    return (1.0 / r - 1.0) * a / atot


def norm_pred(mu, SW, SL, ysum, m):
    """M6 for texture: mean and covariance for a new sample given m same-borehole samples."""
    if m <= 0:
        return mu, SL + SW
    iW = np.linalg.inv(SW)
    Vb = np.linalg.inv(np.linalg.inv(SL) + m * iW)
    return mu + Vb @ (iW @ (ysum - m * mu)), Vb + SW


def mvn_logpdf(y, mean, cov):
    d = y - mean
    _, logdet = np.linalg.slogdet(cov)
    return -0.5 * (d @ np.linalg.solve(cov, d) + logdet + len(y) * np.log(2 * np.pi))


def _floor(pi):
    pi = np.maximum(np.asarray(pi, float), PI_FLOOR)
    return pi / pi.sum()


def _psd(M, floor):
    w, U = np.linalg.eigh((M + M.T) / 2)
    return (U * np.maximum(w, floor)) @ U.T


def _loo_loglik(logAW, logAL, x, n, g, pi_ref):
    """Sum of leave-one-out predictive log probabilities (the quantity the spreads maximise)."""
    AW, AL = np.exp(logAW), np.exp(logAL)
    ne = neff(n, AW)
    p = x / n[:, None]
    E = np.bincount(g, ne)
    Sg = np.column_stack([np.bincount(g, ne * p[:, k]) for k in range(x.shape[1])])
    Eo = np.maximum(E[g] - ne, 0.0)
    So = np.maximum(Sg[g] - ne[:, None] * p, 0.0)
    a = AL * pi_ref + So
    atot = AL + Eo
    r_post, r_w = 1.0 / (atot + 1.0), 1.0 / (AW + 1.0)
    r = r_post + r_w - r_post * r_w
    alpha = (1.0 / r - 1.0)[:, None] * a / atot[:, None]
    return dm_logpmf(x, alpha).sum()


class Model:
    """Fitted on the samples in `fit` (clean, coded till samples of the area).

    unit: array of unit labels (None where a sample has none). Spreads use every fit
    sample in the area (county + buffer); proportions and means use county samples only.
    withheld_views: {view name: boolean array} - measurements of a fit sample that are kept out
    of every spread, proportion and borehole evidence (M7: only the atypical views of an outlier).
    """

    def __init__(self, data, unit, fit, nodes, cfg, withheld_views=None):
        self.cfg = cfg
        self.nodes = nodes
        self.node = {nd.name: nd for nd in nodes}
        self.overall = [nd.name for nd in nodes if nd.overall]
        self.grain_views = [nd.name for nd in nodes if nd.overall and nd.view != "individual"]
        self.hole = data.S.hole.to_numpy()
        self.county = data.S.in_county.to_numpy()
        self.unit = np.asarray(unit, dtype=object)
        self.fit = np.asarray(fit, bool)
        self.tex = data.tex
        self.has_tex = np.isfinite(data.tex).all(axis=1)
        self.excl = {nd.name: np.zeros(len(self.unit), bool) for nd in nodes}
        self.excl["texture"] = np.zeros(len(self.unit), bool)
        for j, m in (withheld_views or {}).items():
            self.excl[j] = self.excl[j] | np.asarray(m, bool)
        from .variables import node_counts
        self.cnt = {nd.name: node_counts(data.raw, nd) for nd in nodes}
        self.n = {j: c.sum(axis=1) for j, c in self.cnt.items()}
        idx = np.flatnonzero(self.fit)
        self.members = {u: idx[self.unit[idx] == u] for u in sorted(set(self.unit[idx]))}
        self.members_county = {u: m[self.county[m]] for u, m in self.members.items()}
        self.units = [u for u, m in self.members_county.items() if len(m)]
        self.AW, self.AL, self.pooled = {}, {}, {}
        self.area_pi = {}
        for nd in nodes:
            self._fit_grain(nd.name)
        self._fit_texture()
        self._cache = {}
        self._fit_weights()

    # ---------------------------------------------------------------- M5: spreads
    def _fit_grain(self, j):
        tau = self.cfg["shrink_tau"]
        prior = self.cfg["prior_grains"]
        idx = np.flatnonzero(self.fit & (self.n[j] > 0) & ~self.excl[j])
        x, n = self.cnt[j][idx], self.n[j][idx]
        self.area_pi[j] = _floor(x.sum(axis=0) + 1e-9 if len(idx) else np.ones(self.cnt[j].shape[1]))
        self.AW[j], self.AL[j] = {}, {}
        if len(idx) == 0:
            self.pooled[j] = (DEFAULT_AW, DEFAULT_AL)
            return
        unit, hole = self.unit[idx], self.hole[idx]
        g = np.unique(np.array([f"{h}|{u}" for h, u in zip(hole, unit)]), return_inverse=True)[1]
        pi_ref = np.zeros_like(x)
        multi = np.zeros(len(idx), bool)
        for u in set(unit):
            m = unit == u
            tot = x[m].sum(axis=0)
            holes_u = set(hole[m])
            multi[m] = len(holes_u) > 1
            for h in holes_u:
                mh = m & (hole == h)
                rest = tot - x[mh].sum(axis=0) if len(holes_u) > 1 else tot
                pi_ref[mh] = _floor(rest + prior * self.area_pi[j])

        def fit(mask, both, prior_log=None):
            if mask.sum() < 2:
                return None
            xm, nm, pm = x[mask], n[mask], pi_ref[mask]
            gm = np.unique(g[mask], return_inverse=True)[1]

            def obj(t):
                lw = np.clip(t[0], LOG_A_MIN, LOG_A_MAX)
                ll_ = np.clip(t[1], LOG_A_MIN, LOG_A_MAX) if both else prior_log[1]
                pen = 0.0
                if prior_log is not None:
                    pen = ((lw - prior_log[0]) ** 2 + (ll_ - prior_log[1]) ** 2) / (2 * tau ** 2)
                return -_loo_loglik(lw, ll_, xm, nm, gm, pm) + pen
            start = list(prior_log) if prior_log is not None else [np.log(DEFAULT_AW), np.log(DEFAULT_AL)]
            res = minimize(obj, start if both else start[:1], method="Nelder-Mead",
                           options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 400})
            lw = np.clip(res.x[0], LOG_A_MIN, LOG_A_MAX)
            ll_ = np.clip(res.x[1], LOG_A_MIN, LOG_A_MAX) if both else prior_log[1]
            return float(lw), float(ll_)

        pooled = fit(multi, True) if multi.sum() >= 10 else None
        pooled = pooled or (np.log(DEFAULT_AW), np.log(DEFAULT_AL))
        self.pooled[j] = (float(np.exp(pooled[0])), float(np.exp(pooled[1])))
        for u in set(unit):
            m = unit == u
            est = fit(m, bool(multi[m].any()), pooled) or pooled
            self.AW[j][u], self.AL[j][u] = float(np.exp(est[0])), float(np.exp(est[1]))

    def _fit_texture(self):
        k = self.cfg["shrink_df"]
        idx = np.flatnonzero(self.fit & self.has_tex & ~self.excl["texture"])
        I = np.eye(2)
        stats = {}
        Wsum, Wdf, Bsum, Bdf, inv_m = np.zeros((2, 2)), 0, np.zeros((2, 2)), 0, []
        for u in set(self.unit[idx]):
            iu = idx[self.unit[idx] == u]
            W, df, means, ms = np.zeros((2, 2)), 0, [], []
            for h in set(self.hole[iu]):
                y = self.tex[iu[self.hole[iu] == h]]
                d = y - y.mean(axis=0)
                W += d.T @ d
                df += len(y) - 1
                means.append(y.mean(axis=0))
                ms.append(len(y))
            means = np.array(means)
            B = (means - means.mean(axis=0)).T @ (means - means.mean(axis=0))
            stats[u] = (W, df, B, len(ms) - 1, float(np.mean(1.0 / np.array(ms))))
            Wsum += W; Wdf += df; Bsum += B; Bdf += len(ms) - 1
            inv_m += [1.0 / m for m in ms]
        SW0 = _psd(Wsum / Wdf, TEX_FLOOR) if Wdf > 0 else 25.0 * I
        SL0 = _psd(Bsum / Bdf - SW0 * np.mean(inv_m), TEX_FLOOR) if Bdf > 0 else 25.0 * I
        self.SW0, self.SL0 = SW0, SL0
        self.SW, self.SL = {}, {}
        for u, (W, df, B, dfb, im) in stats.items():
            self.SW[u] = _psd((W + k * SW0) / (df + k), TEX_FLOOR)
            SL_u = _psd(B / dfb - self.SW[u] * im, TEX_FLOOR) if dfb > 0 else SL0
            self.SL[u] = _psd((dfb * SL_u + k * SL0) / (dfb + k), TEX_FLOOR)
        self.area_mu = self.tex[idx].mean(axis=0) if len(idx) else np.array([40.0, 25.0])

    # ---------------------------------------------------------------- M1: view weights
    def _fit_weights(self, R=2000):
        """Weight of each view in the overall number: 1 / (sum of its positive overlaps).

        Overlap = rank correlation between two views' log ratios when the true unit is fixed:
        samples are simulated from each unit (eight broad groups, a new borehole, the unit's own
        grain totals) and compared with every other unit and with the area average. Only
        counting shared between views can correlate them. Texture and individual-grain views
        share no counts with any other view, so their weight is 1.
        """
        from .variables import leaf_map
        self.weight = {j: 1.0 for j in self.overall}
        self.weight["texture"] = 1.0
        g = self.grain_views
        self.overlap = np.eye(len(g))
        if len(g) < 2 or not self.units:
            return
        rng = np.random.default_rng(self.cfg["seed"])
        maps = {j: leaf_map(self.node[j]) for j in g}
        names = g + ["leaf"]

        def alphas(pi, spread):
            return {j: dm_pred_alpha(pi[j], *spread[j], np.zeros(len(pi[j])), 0.0) for j in names}
        alts = {u: alphas(self.unit_params(u)["pi"], {j: self.spread(j, u) for j in names}) for u in self.units}
        alts["_area"] = alphas(self.area_pi, {j: self.pooled[j] for j in names})
        blocks = []
        for u in self.units:
            n_real = self.n["leaf"][self.members[u]]
            n_real = n_real[n_real > 0]
            n = rng.choice(n_real, R) if len(n_real) else np.full(R, 300.0)
            theta = rng.dirichlet(np.maximum(alts[u]["leaf"], 1e-12), size=R)
            x = rng.multinomial(n.astype(np.int64), theta)
            X = {j: x @ maps[j] for j in g}
            own = {j: dm_logpmf(X[j], alts[u][j]) for j in g}
            for v, av in alts.items():
                if v == u:
                    continue
                L = np.column_stack([dm_logpmf(X[j], av[j]) - own[j] for j in g])
                blocks.append(L - L.mean(axis=0))
        Z = np.vstack(blocks)
        ranks = np.argsort(np.argsort(Z, axis=0), axis=0).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            rho = np.corrcoef(ranks, rowvar=False)
        rho = np.where(np.isfinite(rho), rho, 0.0)
        np.fill_diagonal(rho, 1.0)
        self.overlap = rho
        w = 1.0 / np.clip(rho, 0.0, None).sum(axis=1)
        self.weight.update({j: float(wj) for j, wj in zip(g, w)})

    def weight_table(self):
        """Rows for the report: each view, its weight and its overlap with the grain views."""
        rows = []
        for j in self.overall + ["texture"]:
            row = {"view": j, "weight": self.weight[j]}
            for k, other in enumerate(self.grain_views):
                row[f"overlap_with_{other}"] = (float(self.overlap[self.grain_views.index(j), k])
                                                if j in self.grain_views else 0.0)
            rows.append(row)
        return rows

    def combine(self, parts):
        """Weighted sum of per-view log probabilities (numbers or arrays): the overall."""
        return sum(self.weight[j] * parts[j] for j in self.overall) + parts["texture"]

    def spread(self, j, u):
        return self.AW[j].get(u, self.pooled[j][0]), self.AL[j].get(u, self.pooled[j][1])

    # ---------------------------------------------------------------- means
    def pool_pi(self, j, idxs, AW, AL, weights=None):
        """Borehole-weighted proportions of a branch over samples idxs (plus a small prior)."""
        idxs = np.asarray(idxs, int)
        w = np.ones(len(idxs)) if weights is None else np.asarray(weights, float)
        keep = (self.n[j][idxs] > 0) & (w > 0) & ~self.excl[j][idxs]
        idxs, w = idxs[keep], w[keep]
        prior = self.cfg["prior_grains"]
        if len(idxs) == 0:
            return self.area_pi[j]
        x, n = self.cnt[j][idxs], self.n[j][idxs]
        ne = neff(n, AW) * w
        hid = np.unique(self.hole[idxs], return_inverse=True)[1]
        E = np.bincount(hid, ne)
        P = np.column_stack([np.bincount(hid, ne * x[:, k] / n) for k in range(x.shape[1])]) / E[:, None]
        wh = 1.0 / (1.0 / (AL + 1.0) + 1.0 / E)
        pw = (wh[:, None] * P).sum(axis=0) / wh.sum()
        return _floor((E.sum() * pw + prior * self.area_pi[j]) / (E.sum() + prior))

    def pool_mu(self, idxs, SW, SL, weights=None):
        idxs = np.asarray(idxs, int)
        w = np.ones(len(idxs)) if weights is None else np.asarray(weights, float)
        keep = self.has_tex[idxs] & (w > 0) & ~self.excl["texture"][idxs]
        idxs, w = idxs[keep], w[keep]
        if len(idxs) == 0:
            return self.area_mu
        num, den = np.zeros(2), np.zeros((2, 2))
        for h in set(self.hole[idxs]):
            m = self.hole[idxs] == h
            mh = w[m].sum()
            Wh = np.linalg.inv(SL + SW / mh)
            num += Wh @ ((w[m][:, None] * self.tex[idxs[m]]).sum(axis=0) / mh)
            den += Wh
        return np.linalg.solve(den, num)

    def unit_params(self, u, exclude=()):
        """Proportions, means and spreads of unit u, leaving out the samples in `exclude`."""
        key = (u, tuple(sorted(exclude)))
        if key in self._cache:
            return self._cache[key]
        idxs = self.members_county.get(u, np.array([], int))
        if len(exclude):
            idxs = idxs[~np.isin(idxs, list(exclude))]
        P = {"pi": {}, "AW": {}, "AL": {}}
        for nd in self.nodes:
            j = nd.name
            P["AW"][j], P["AL"][j] = self.spread(j, u)
            P["pi"][j] = self.pool_pi(j, idxs, P["AW"][j], P["AL"][j])
        P["SW"], P["SL"] = self.SW.get(u, self.SW0), self.SL.get(u, self.SL0)
        P["mu"] = self.pool_mu(idxs, P["SW"], P["SL"])
        if not exclude:
            self._cache[key] = P
        return P

    # ---------------------------------------------------------------- M6
    def chain(self, idxs, P, weights=None):
        """Joint log probability of samples idxs as one borehole of a unit with params P.

        Sequential (each sample given those before it). Returns {branch: per-sample terms},
        including 'texture'. weights scale a sample's term and what it teaches the borehole.
        """
        idxs = list(idxs)
        w = np.ones(len(idxs)) if weights is None else np.asarray(weights, float)
        out = {}
        for nd in self.nodes:
            j = nd.name
            pi, AW, AL = P["pi"][j], P["AW"][j], P["AL"][j]
            S, E, terms = np.zeros(len(pi)), 0.0, np.zeros(len(idxs))
            for t, i in enumerate(idxs):
                x, n = self.cnt[j][i], self.n[j][i]
                if n > 0 and not self.excl[j][i]:
                    terms[t] = w[t] * dm_logpmf(x, dm_pred_alpha(pi, AW, AL, S, E))
                    ne = neff(n, AW) * w[t]
                    S = S + ne * x / n
                    E += ne
            out[j] = terms
        ysum, m, terms = np.zeros(2), 0.0, np.zeros(len(idxs))
        for t, i in enumerate(idxs):
            if self.has_tex[i] and not self.excl["texture"][i]:
                mean, cov = norm_pred(P["mu"], P["SW"], P["SL"], ysum, m)
                terms[t] = w[t] * mvn_logpdf(self.tex[i], mean, cov)
                ysum = ysum + w[t] * self.tex[i]
                m += w[t]
        out["texture"] = terms
        return out

    def sample_params(self, i, u, others=None):
        """Predictive parameters for sample i under unit u in its own borehole (never uses i)."""
        P = self.unit_params(u, exclude=(i,) if (self.unit[i] == u and self.fit[i]) else ())
        if others is None:
            mem = self.members.get(u, np.array([], int))
            others = mem[(self.hole[mem] == self.hole[i]) & (mem != i)]
        out = {"alpha": {}, "others": len(others)}
        for nd in self.nodes:
            j = nd.name
            AW, AL = P["AW"][j], P["AL"][j]
            S, E = np.zeros(len(P["pi"][j])), 0.0
            for o in others:
                if self.n[j][o] > 0 and not self.excl[j][o]:
                    ne = neff(self.n[j][o], AW)
                    S = S + ne * self.cnt[j][o] / self.n[j][o]
                    E += ne
            out["alpha"][j] = dm_pred_alpha(P["pi"][j], AW, AL, S, E)
        ot = [o for o in others if self.has_tex[o] and not self.excl["texture"][o]]
        ysum = self.tex[ot].sum(axis=0) if ot else np.zeros(2)
        out["mean"], out["cov"] = norm_pred(P["mu"], P["SW"], P["SL"], ysum, len(ot))
        return out

    def logp(self, i, pred):
        """{branch: log probability of sample i} under predictive parameters `pred`."""
        out = {j: (dm_logpmf(self.cnt[j][i], a) if self.n[j][i] > 0 else 0.0) for j, a in pred["alpha"].items()}
        out["texture"] = mvn_logpdf(self.tex[i], pred["mean"], pred["cov"]) if self.has_tex[i] else 0.0
        return out

    def total(self, lp):
        return self.combine(lp)

    # ---------------------------------------------------------------- M7
    def typicality(self, i, pred, rng, R=999):
        """Share of simulated samples (same grain total) no more probable than sample i.

        Simulated samples are drawn as eight broad groups (`leaf`) so the views keep the counts
        they share, then scored view by view and combined with the M1 weights. Each view's own p
        uses simulations with that view's observed total. Returns p overall, p for grains, p for
        texture and the view with the smallest p.
        """
        from .variables import LEAVES, leaf_map
        obs = self.logp(i, pred)
        branch_p = {}
        n_leaf = int(self.n["leaf"][i])
        has_g = n_leaf > 0
        sims = np.zeros(R)
        if has_g:
            x = rng.multinomial(n_leaf, rng.dirichlet(np.maximum(pred["alpha"]["leaf"], 1e-12), size=R))
            for j in self.overall:
                nd = self.node[j]
                if nd.view == "individual":
                    n = x[:, LEAVES.index(nd.within)].astype(np.int64)
                    xj = rng.multinomial(n, rng.dirichlet(np.maximum(pred["alpha"][j], 1e-12), size=R))
                else:
                    xj = x @ leaf_map(nd)
                sims += self.weight[j] * dm_logpmf(xj, pred["alpha"][j])
        for j in self.overall:
            if self.n[j][i] > 0:
                theta = rng.dirichlet(np.maximum(pred["alpha"][j], 1e-12), size=R)
                xs = rng.multinomial(int(self.n[j][i]), theta)
                branch_p[j] = (1 + np.sum(dm_logpmf(xs, pred["alpha"][j]) <= obs[j] + 1e-9)) / (R + 1)
        obs_g = sum(self.weight[j] * obs[j] for j in self.overall)
        p_g = (1 + np.sum(sims <= obs_g + 1e-9)) / (R + 1) if has_g else np.nan
        p_t = np.nan
        if self.has_tex[i]:
            d = self.tex[i] - pred["mean"]
            p_t = float(chi2.sf(d @ np.linalg.solve(pred["cov"], d), 2))
            branch_p["texture"] = p_t
            ysim = rng.multivariate_normal(pred["mean"], pred["cov"], size=R)
            dd = ysim - pred["mean"]
            _, logdet = np.linalg.slogdet(pred["cov"])
            sims_t = -0.5 * (np.einsum("ij,ij->i", dd @ np.linalg.inv(pred["cov"]), dd) + logdet + 2 * np.log(2 * np.pi))
            p = (1 + np.sum((sims if has_g else 0) + sims_t <= self.total(obs) + 1e-9)) / (R + 1)
        else:
            p = p_g
        worst = min(branch_p, key=branch_p.get) if branch_p else ""
        return {"p": float(p), "p_grains": float(p_g), "p_texture": p_t, "worst_branch": worst,
                "branch_p": branch_p}
