from __future__ import annotations

import itertools
import math
from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gamma, weibull_min


FACTORS = ("G", "T", "S", "C")
RELEASE_TYPES = {
    "compatibility_release",
    "tj_compatibility_release",
    "climb_quota_completion",
}


def mechanism_set(regime: str) -> frozenset[str]:
    """Return the active G/T/S/C letters for one factorial regime."""
    if regime == "B0":
        return frozenset()
    return frozenset(letter for letter in FACTORS if letter in regime)


def release_events(events: pd.DataFrame) -> pd.DataFrame:
    if events.empty or "event_type" not in events:
        return events.iloc[0:0].copy()
    return events[events["event_type"].isin(RELEASE_TYPES)].copy()


def parse_grain_ids(value: object) -> tuple[int, ...]:
    if value is None or pd.isna(value):
        return ()
    result = []
    for item in str(value).replace("[", "").replace("]", "").split(";"):
        try:
            result.append(int(item.strip()))
        except ValueError:
            continue
    return tuple(result)


def _motion_intervals(tracks: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for grain_id, group in tracks.groupby("grain_id", sort=False):
        group = group.sort_values("time")
        time = group["time"].to_numpy(float)
        radius = group["radius"].to_numpy(float)
        if len(time) < 2:
            continue
        dt = np.diff(time)
        valid = dt > 0
        if not np.any(valid):
            continue
        rows.append(pd.DataFrame({
            "grain_id": int(grain_id),
            "time": time[1:][valid],
            "dt": dt[valid],
            "motion": np.abs(np.diff(radius)[valid]) / dt[valid],
        }))
    if not rows:
        return pd.DataFrame(columns=["grain_id", "time", "dt", "motion"])
    return pd.concat(rows, ignore_index=True)


def event_burst_coupling(
    tracks: pd.DataFrame,
    events: pd.DataFrame,
    *,
    shuffle_samples: int = 200,
    seed: int = 20260821,
) -> dict[str, float]:
    """Measure grain-scale release/motion association against a stratified null.

    A release is attached to each grain named in its ledger row.  Motion is the
    absolute radius rate over one recorded interval.  The null preserves, for
    every grain, the number of release-associated intervals while randomly
    assigning those intervals across that grain's observed trajectory.
    """
    intervals = _motion_intervals(tracks)
    releases = release_events(events)
    result = {
        "motion_interval_count": int(len(intervals)),
        "release_count": int(len(releases)),
        "release_grain_links": 0,
        "event_window": np.nan,
        "p_large_given_release": np.nan,
        "p_large_unconditional": np.nan,
        "large_burst_risk_ratio": np.nan,
        "large_burst_risk_ratio_shuffled": np.nan,
        "large_burst_risk_ratio_excess": 0.0,
        "event_growth_xcorr_max": np.nan,
        "event_growth_xcorr_max_shuffled": np.nan,
        "event_growth_xcorr_excess": 0.0,
    }
    for pct in (1, 5, 10):
        result[f"top_{pct}pct_preceded_by_release"] = 0.0
        result[f"top_{pct}pct_preceded_by_release_shuffled"] = 0.0
        result[f"top_{pct}pct_release_excess"] = 0.0
    if intervals.empty:
        return result

    cadence = float(np.median(intervals["dt"].to_numpy(float)))
    result["event_window"] = cadence
    motion = intervals["motion"].to_numpy(float)
    positive_count = int(np.count_nonzero(motion > 0.0))
    top_masks = {}
    for pct in (1, 5, 10):
        count = min(positive_count, max(1, int(math.ceil(len(motion) * pct / 100.0))))
        mask = np.zeros(len(motion), dtype=bool)
        if count:
            indices = np.argpartition(motion, -count)[-count:]
            mask[indices] = motion[indices] > 0.0
        top_masks[pct] = mask
    large = top_masks[5]
    result["p_large_unconditional"] = float(np.mean(large))

    event_times: dict[int, list[float]] = {}
    for row in releases.itertuples(index=False):
        for grain_id in parse_grain_ids(getattr(row, "grain_ids", None)):
            event_times.setdefault(grain_id, []).append(float(row.time))
    result["release_grain_links"] = int(sum(map(len, event_times.values())))
    preceded = np.zeros(len(intervals), dtype=bool)
    group_indices: list[np.ndarray] = []
    group_counts: list[int] = []
    for grain_id, raw_indices in intervals.groupby("grain_id", sort=False).indices.items():
        indices = np.asarray(raw_indices, dtype=int)
        times = intervals.iloc[indices]["time"].to_numpy(float)
        grain_event_times = np.sort(np.asarray(event_times.get(int(grain_id), []), float))
        if len(grain_event_times):
            right = np.searchsorted(grain_event_times, times, side="right")
            last = np.full(len(times), -np.inf)
            valid = right > 0
            last[valid] = grain_event_times[right[valid] - 1]
            local = valid & (last > times - cadence - 1e-12)
            preceded[indices] = local
        count = int(np.count_nonzero(preceded[indices]))
        if count:
            group_indices.append(indices)
            group_counts.append(count)

    if not np.any(preceded):
        return result

    actual_risk = float(np.mean(large[preceded]))
    baseline_risk = max(float(np.mean(large)), np.finfo(float).tiny)
    actual_ratio = actual_risk / baseline_risk
    result["p_large_given_release"] = actual_risk
    result["large_burst_risk_ratio"] = actual_ratio

    for pct, mask in top_masks.items():
        result[f"top_{pct}pct_preceded_by_release"] = (
            float(np.mean(preceded[mask])) if np.any(mask) else 0.0
        )

    rng = np.random.default_rng(seed)
    shuffled_risk = []
    shuffled_top = {pct: [] for pct in top_masks}
    for _ in range(shuffle_samples):
        null = np.zeros(len(intervals), dtype=bool)
        for indices, count in zip(group_indices, group_counts):
            null[rng.choice(indices, size=count, replace=False)] = True
        shuffled_risk.append(float(np.mean(large[null])) / baseline_risk)
        for pct, mask in top_masks.items():
            shuffled_top[pct].append(float(np.mean(null[mask])) if np.any(mask) else 0.0)
    null_ratio = float(np.mean(shuffled_risk))
    result["large_burst_risk_ratio_shuffled"] = null_ratio
    result["large_burst_risk_ratio_excess"] = actual_ratio - null_ratio
    for pct, values in shuffled_top.items():
        null_value = float(np.mean(values))
        actual = result[f"top_{pct}pct_preceded_by_release"]
        result[f"top_{pct}pct_preceded_by_release_shuffled"] = null_value
        result[f"top_{pct}pct_release_excess"] = actual - null_value

    # System-level event-rate/growth-rate cross-correlation.  Positive lag means
    # release counts lead motion by that many recorded frames.
    frame_times = np.sort(intervals["time"].unique())
    motion_by_time = intervals.groupby("time")["motion"].mean().reindex(frame_times).to_numpy(float)
    release_times = releases["time"].to_numpy(float)
    edges = np.r_[frame_times[0] - cadence, frame_times]
    event_counts = np.histogram(release_times, bins=edges)[0].astype(float)

    def max_leading_correlation(counts: np.ndarray) -> float:
        correlations = []
        for lag in range(4):
            x = counts[:len(counts) - lag or None]
            y = motion_by_time[lag:]
            if len(x) > 2 and np.std(x) > 0 and np.std(y) > 0:
                correlations.append(float(np.corrcoef(x, y)[0, 1]))
        return max(correlations) if correlations else np.nan

    actual_xcorr = max_leading_correlation(event_counts)
    null_xcorr = []
    for _ in range(shuffle_samples):
        null_xcorr.append(max_leading_correlation(rng.permutation(event_counts)))
    finite_null = np.asarray(null_xcorr, float)
    finite_null = finite_null[np.isfinite(finite_null)]
    mean_null = float(np.mean(finite_null)) if len(finite_null) else np.nan
    result["event_growth_xcorr_max"] = actual_xcorr
    result["event_growth_xcorr_max_shuffled"] = mean_null
    result["event_growth_xcorr_excess"] = (
        actual_xcorr - mean_null if np.isfinite(actual_xcorr) and np.isfinite(mean_null) else 0.0
    )
    return result


def factorial_effects(table: pd.DataFrame, metrics: Iterable[str]) -> pd.DataFrame:
    """Calculate conventional high-minus-low effects for a complete 2^4 design."""
    lookup = {frozenset(mechanism_set(str(row.regime))): row for row in table.itertuples()}
    expected = {frozenset(combo) for size in range(5) for combo in itertools.combinations(FACTORS, size)}
    missing = expected - set(lookup)
    if missing:
        raise ValueError(f"factorial design is incomplete; missing {sorted(map(sorted, missing))}")
    rows = []
    for metric in metrics:
        for order in range(1, len(FACTORS) + 1):
            for interaction in itertools.combinations(FACTORS, order):
                contrast = 0.0
                values = []
                for active, row in lookup.items():
                    value = float(getattr(row, metric))
                    values.append(value)
                    sign = math.prod(1.0 if factor in active else -1.0 for factor in interaction)
                    contrast += sign * value
                effect = contrast / (2 ** (len(FACTORS) - 1))
                rows.append({
                    "metric": metric,
                    "effect": "x".join(interaction),
                    "order": order,
                    "value": effect,
                    "absolute_value": abs(effect),
                    "response_range": max(values) - min(values),
                })
    result = pd.DataFrame(rows)
    result["absolute_rank_within_metric"] = result.groupby("metric")["absolute_value"].rank(
        method="min", ascending=False
    ).astype(int)
    return result.sort_values(["metric", "absolute_rank_within_metric", "order", "effect"])


def pin_episodes(boundaries: pd.DataFrame) -> pd.DataFrame:
    """Extract contiguous blocked-domain episodes with right-censoring flags."""
    rows = []
    for entity_id, group in boundaries.groupby("entity_id", sort=False):
        group = group.sort_values("time")
        time = group["time"].to_numpy(float)
        blocked = group["blocked"].to_numpy(bool)
        if not len(time) or not np.any(blocked):
            continue
        gaps = np.diff(time)
        cadence = float(np.median(gaps[gaps > 0])) if np.any(gaps > 0) else 0.0
        start = None
        for index, value in enumerate(blocked):
            if value and start is None:
                start = index
            terminal = start is not None and (not value or index == len(blocked) - 1)
            if terminal:
                censored = bool(value and index == len(blocked) - 1)
                last = index if censored else index - 1
                end = time[last] + cadence
                rows.append({
                    "entity_id": entity_id,
                    "start": float(time[start]),
                    "end": float(end),
                    "duration": float(max(0.0, end - time[start])),
                    "censored": censored,
                })
                start = None
    return pd.DataFrame(rows, columns=["entity_id", "start", "end", "duration", "censored"])


def activation_wait_episodes(
    events: pd.DataFrame,
    *,
    hit_type: str,
    release_type: str,
    end_time: float,
) -> pd.DataFrame:
    """Return lower-bound waits from first persisted activation hit to release."""
    selected = events[events["event_type"].isin({hit_type, release_type})]
    rows = []
    for entity_id, group in selected.groupby("entity_id", sort=False):
        start = None
        for row in group.sort_values("time").itertuples(index=False):
            if row.event_type == hit_type and start is None:
                start = float(row.time)
            elif row.event_type == release_type and start is not None:
                rows.append({
                    "entity_id": entity_id, "start": start, "end": float(row.time),
                    "duration": max(0.0, float(row.time) - start), "censored": False,
                })
                start = None
        if start is not None:
            rows.append({
                "entity_id": entity_id, "start": start, "end": float(end_time),
                "duration": max(0.0, float(end_time) - start), "censored": True,
            })
    return pd.DataFrame(rows, columns=["entity_id", "start", "end", "duration", "censored"])


def kaplan_meier(episodes: pd.DataFrame) -> pd.DataFrame:
    if episodes.empty:
        return pd.DataFrame(columns=["duration", "survival", "at_risk", "events", "censored"])
    frame = episodes.sort_values("duration")
    survival = 1.0
    rows = []
    for duration in np.sort(frame["duration"].unique()):
        at_risk = int(np.count_nonzero(frame["duration"].to_numpy(float) >= duration))
        at_time = frame[frame["duration"] == duration]
        failures = int(np.count_nonzero(~at_time["censored"].to_numpy(bool)))
        censored = int(np.count_nonzero(at_time["censored"].to_numpy(bool)))
        if at_risk:
            survival *= 1.0 - failures / at_risk
        rows.append({
            "duration": float(duration), "survival": survival, "at_risk": at_risk,
            "events": failures, "censored": censored,
        })
    return pd.DataFrame(rows)


def fit_pin_lifetime_models(episodes: pd.DataFrame) -> pd.DataFrame:
    """Fit censored exponential, gamma, and Weibull waiting-time models."""
    if episodes.empty:
        return pd.DataFrame()
    duration = np.maximum(episodes["duration"].to_numpy(float), np.finfo(float).tiny)
    observed = ~episodes["censored"].to_numpy(bool)
    failures = int(np.count_nonzero(observed))
    if failures == 0:
        return pd.DataFrame()
    rows = []
    rate = failures / float(np.sum(duration))
    ll = float(np.sum(observed * np.log(rate) - rate * duration))
    rows.append({"model": "exponential", "shape": 1.0, "scale": 1.0 / rate,
                 "log_likelihood": ll, "parameters": 1, "AIC": 2.0 - 2.0 * ll})

    def fit_model(name: str):
        distribution = weibull_min if name == "weibull" else gamma

        def objective(theta: np.ndarray) -> float:
            shape, scale = np.exp(theta)
            logp = distribution.logpdf(duration, shape, loc=0.0, scale=scale)
            logs = distribution.logsf(duration, shape, loc=0.0, scale=scale)
            value = np.where(observed, logp, logs)
            return -float(np.sum(value)) if np.all(np.isfinite(value)) else 1e300

        initial = np.log([1.0, max(float(np.median(duration)), np.finfo(float).tiny)])
        fitted = minimize(objective, initial, method="Nelder-Mead")
        shape, scale = np.exp(fitted.x)
        log_likelihood = -float(fitted.fun)
        rows.append({"model": name, "shape": float(shape), "scale": float(scale),
                     "log_likelihood": log_likelihood, "parameters": 2,
                     "AIC": 4.0 - 2.0 * log_likelihood})

    fit_model("gamma")
    fit_model("weibull")
    result = pd.DataFrame(rows).sort_values("AIC").reset_index(drop=True)
    result["delta_AIC"] = result["AIC"] - result["AIC"].min()
    result["preferred"] = result.index == 0
    result["episode_count"] = len(episodes)
    result["observed_releases"] = failures
    result["censored_count"] = len(episodes) - failures
    return result


def arrhenius_work_summary(work: pd.DataFrame, temperature: float) -> pd.DataFrame:
    """Summarize signed rate multipliers and magnitude-of-change tail fractions."""
    if work.empty:
        return pd.DataFrame()
    kbt = 8.617333262145e-5 * float(temperature)
    rows = []
    for mechanism, column in (("S", "work_shear"), ("C", "work_free_volume")):
        all_values = work[column].to_numpy(float)
        populations = (("all_events", all_values), ("nonzero_work", all_values[np.abs(all_values) > 1e-15]))
        for population, values in populations:
            if not len(values):
                continue
            row = {
                "mechanism": mechanism, "population": population,
                "events": len(values), "all_events": len(all_values), "kBT_eV": kbt,
            }
            for percentile in (50, 90, 95, 99):
                value = float(np.quantile(values, percentile / 100.0))
                row[f"work_p{percentile}_eV"] = value
                row[f"multiplier_p{percentile}"] = float(np.exp(np.clip(value / kbt, -700, 700)))
            magnitude = np.exp(np.clip(np.abs(values) / kbt, 0, 700))
            signed = np.exp(np.clip(values / kbt, -700, 700))
            for threshold in (2, 5, 10):
                row[f"fraction_change_gt_{threshold}x"] = float(np.mean(magnitude > threshold))
                row[f"fraction_accelerate_gt_{threshold}x"] = float(np.mean(signed > threshold))
            rows.append(row)
    return pd.DataFrame(rows)


def reconstruct_gb_occupancy(
    episodes: pd.DataFrame,
    events: pd.DataFrame,
    *,
    has_g: bool,
    has_c: bool,
) -> dict[str, float]:
    """Reconstruct G/C-limited domain-time from blocked tracks and completions.

    Completion events identify which gates were active in a pin episode.  When
    both complete, their times split the episode into multiple- and single-gate
    portions.  Missing completions remain explicitly unresolved.
    """
    totals = {"G_limited": 0.0, "C_limited": 0.0, "multiple": 0.0, "unresolved": 0.0}
    if episodes.empty:
        return {**totals, "total_blocked_domain_time": 0.0}
    completed = events[events["event_type"].isin({"compatibility_release", "climb_quota_completion"})]
    by_entity = {key: group.sort_values("time") for key, group in completed.groupby("entity_id")}
    for episode in episodes.itertuples(index=False):
        start, end = float(episode.start), float(episode.end)
        if not (has_g or has_c):
            totals["unresolved"] += end - start
            continue
        group = by_entity.get(episode.entity_id)
        times = {}
        if group is not None:
            inside = group[(group["time"] >= start) & (group["time"] <= end + 1e-12)]
            for event_type, label in (("compatibility_release", "G"), ("climb_quota_completion", "C")):
                found = inside.loc[inside["event_type"] == event_type, "time"]
                if len(found):
                    times[label] = float(found.iloc[0])
        # A module being enabled does not mean its gate was pending in every
        # episode.  A completion inside the episode is the persisted evidence
        # that identifies an active gate.  No-completion episodes stay unresolved.
        active = set(times)
        if not active:
            totals["unresolved"] += end - start
            continue
        cursor = start
        for label, completion in sorted(times.items(), key=lambda item: item[1]):
            completion = min(end, max(cursor, completion))
            key = "multiple" if len(active) > 1 else f"{next(iter(active))}_limited" if active else "unresolved"
            totals[key] += completion - cursor
            active.discard(label)
            cursor = completion
        key = "multiple" if len(active) > 1 else f"{next(iter(active))}_limited" if active else "unresolved"
        totals[key] += end - cursor
    total = float(sum(totals.values()))
    result = {**totals, "total_blocked_domain_time": total}
    for key in tuple(totals):
        result[f"fraction_{key}"] = totals[key] / total if total else np.nan
    return result
