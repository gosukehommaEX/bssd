"""Collect the reference values of ref_exptime.json, ref_binbeta.json,
ref_canary.json, ref_mapbern.json and ref_multibin.json into
ref_methods.csv (columns id, quantity, reference), which
inst/validation/validate_methods_python.R reads.
"""
import csv
import json


def rows_of(name, prefix, mapping):
    res = json.load(open(name))
    out = []
    for r in res["results"]:
        for key, label in mapping:
            if key not in r:
                continue
            v = r[key]
            if isinstance(v, list):
                for lab, x in zip(("lower", "median", "upper"), v):
                    out.append(["%s_%s_%s_%s" % (prefix, r["id"], key, lab),
                                "%s %s" % (label, lab), x])
            else:
                out.append(["%s_%s_%s" % (prefix, r["id"], key), label, v])
    return out


rows = []
rows += rows_of("ref_exptime.json", "PY_EXPTIME", [
    ("prob_rate", "P(lambda_d > rate0)"), ("prob_ratio", "P(r > rr0)"),
    ("mean_lambda_c", "posterior mean of lambda_c"),
    ("mean_lambda_d", "posterior mean of lambda_d"),
    ("mean_r", "posterior mean of r"),
    ("q_lambda_d", "quantile of lambda_d"), ("q_r", "quantile of r")])
rows += rows_of("ref_binbeta.json", "PY_BINBETA", [
    ("prob_rate", "P(theta_d > p0)"), ("prob_ratio", "P(rr > rr0)"),
    ("mean_theta_c", "posterior mean of theta_c"),
    ("mean_theta_d", "posterior mean of theta_d"),
    ("mean_rr", "posterior mean of rr"),
    ("q_theta_d", "quantile of theta_d"), ("q_rr", "quantile of rr")])
rows += rows_of("ref_canary.json", "PY_CANARY", [
    ("prob", "P(M > crit)"), ("mean_metric", "posterior mean of M"),
    ("mean_p_c", "posterior mean of p_c"),
    ("mean_lambda_c", "posterior mean of lambda_c"),
    ("q", "quantile of M")])
rows += rows_of("ref_mapbern.json", "PY_MAPBERN", [
    ("prob_ratio", "P(pi_d > rr0 pi_c)"),
    ("mean_pi_c", "posterior mean of pi_c"),
    ("mean_pi_d", "posterior mean of pi_d"),
    ("q_pi_d", "quantile of pi_d"), ("q_rr", "quantile of rr")])
rows += rows_of("ref_multibin.json", "PY_MULTIBIN", [
    ("prob", "P(delta > rr0)"),
    ("mean_log_delta", "posterior mean of log(delta)"),
    ("q", "quantile of delta")])
mb = json.load(open("ref_mapbern.json"))
for r in mb["map_results"]:
    rows.append(["PY_MAP_%s_mean" % r["id"], "MAP prior mean", r["mean"]])
    rows.append(["PY_MAP_%s_sd" % r["id"], "MAP prior SD", r["sd"]])
    for x, c in zip(r["cdf_at"], r["cdf"]):
        rows.append(["PY_MAP_%s_cdf_%g" % (r["id"], x),
                     "MAP prior distribution function at %g" % x, c])
rows.append(["PY_MAP_SUBJ_mean", "MAP prior mean, subject-level data",
             mb["map_subject"]["mean"]])
rows.append(["PY_MAP_SUBJ_sd", "MAP prior SD, subject-level data",
             mb["map_subject"]["sd"]])
with open("ref_methods.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "quantity", "reference"])
    for r in rows:
        w.writerow([r[0], r[1], repr(float(r[2]))])
print(len(rows), "rows")

# data sets of the subject-level cases
import ref_common as rc
data = rc.make_data()
with open("data_methods.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["dataset", "y", "t", "study"])
    for name in ("E1", "E2", "E3", "M1"):
        for yy, tt in zip(data[name]["y"], data[name]["t"]):
            w.writerow([name, int(yy), repr(float(tt)), 1])
    d = mb["map_subject_data"]
    for yy, tt, ss in zip(d["y"], d["t"], d["study"]):
        w.writerow(["MAP_SUBJ", int(yy), repr(float(tt)), int(ss)])
print("data_methods.csv written")
