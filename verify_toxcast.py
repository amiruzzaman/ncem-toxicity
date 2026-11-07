"""Independently reconcile a published ToxCast snapshot against native result IDs."""
import argparse
import csv
import json
import platform
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

from download_toxcast import fingerprint


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        yield from csv.DictReader(stream, strict=True)


def validate(snapshot, reviewed_master):
    manifest = json.loads((snapshot / "toxcast_import.json").read_text())
    for name, digest in manifest["output_sha256"].items():
        if fingerprint(snapshot / name) != digest:
            raise ValueError(f"Published output hash mismatch: {name}")
    ids, original_ids, tox21 = set(), set(), set()
    summaries = {}
    current = rows(snapshot / "master_chemicals.csv")
    for original in rows(reviewed_master):
        row = next(current, None)
        if row is None or any(row.get(k) != v for k,v in original.items()):
            raise ValueError("Original master values, IDs, order, or reviewed keys changed")
        cid = row["chemical_id"]
        if cid in ids:
            raise ValueError("Duplicate master ID")
        ids.add(cid); original_ids.add(cid)
        if any(v != "" for k,v in original.items() if k.startswith(("NR-","SR-"))):
            tox21.add(cid)
        counts = tuple(int(row[k]) for k in ("toxcast_mc_result_count","toxcast_sc_result_count","toxcast_endpoint_count"))
        if any(counts):summaries[cid] = counts
    for row in current:
        cid = row["chemical_id"]
        if cid in ids or row["identity_status"] != "toxcast_substance_only" or row["standard_inchikey"]:
            raise ValueError("Invalid newly added substance identity")
        ids.add(cid)
        counts = tuple(int(row[k]) for k in ("toxcast_mc_result_count","toxcast_sc_result_count","toxcast_endpoint_count"))
        if any(counts):summaries[cid] = counts
    counts, endpoints = Counter(), defaultdict(set)
    catalog = {str(int(Decimal(r["aeid"]))):r["annotation_status"] for r in rows(snapshot/"toxcast_endpoints.csv")}
    documented = {str(int(Decimal(r["aeid"]))) for r in rows(snapshot/"toxcast_assays.csv")}
    if any((status=="documented") != (aeid in documented) for aeid,status in catalog.items()) or not documented <= set(catalog):
        raise ValueError("Complete endpoint catalog disagrees with released annotations")
    links=set()
    for row in rows(snapshot/"toxcast_targets.csv"):
        aeid=str(int(Decimal(row["aeid"])))
        if aeid not in catalog:raise ValueError("Dangling target annotation")
        if row["target_id"].strip().upper() not in {"","NA","N/A","#N/A","NULL","NONE","NAN"}:
            links.add((aeid,row["target_type"],str(int(Decimal(row["target_id"])))))
    annotated={t[0] for t in links}
    results_with_target=0
    interpretable, present = set(), set()
    result_ids = {"mc":set(),"sc":set()}
    status_counts, raw_counts = Counter(), Counter()
    unresolved = 0
    for screen, id_field in (("mc","m5id"),("sc","s2id")):
        for row in rows(snapshot / f"toxcast_{screen}_observations.csv"):
            raw_counts[screen] += 1
            rid = str(int(Decimal(row[id_field])))
            if row["result_id"] != rid:
                raise ValueError("Native processing ID differs from exported result ID")
            # Re-evaluate scientific interpretation from native fields, not the exported status.
            raw = row.get("hitc",row.get("hitcall",""))
            if raw.strip().upper() in {"","NA","N/A","#N/A","NULL","NONE","NAN"}:
                expected="missing"
            else:
                try:
                    v=Decimal(raw)
                    if not v.is_finite():expected="invalid"
                    elif screen=="mc" and row.get("modl","").lower() in {"none","nofit"}:expected="unfit"
                    elif v<0:expected="negative_direction" if screen=="mc" else "unfit"
                    elif v>1 or (screen=="sc" and v not in (0,1)):expected="invalid"
                    else:expected="active" if v >= (Decimal("0.9") if screen=="mc" else 1) else "inactive"
                except InvalidOperation:expected="invalid"
            if row["activity_status"] != expected:
                raise ValueError(f"Incorrect activity interpretation at {screen} result {rid}")
            if rid in result_ids[screen]:continue
            result_ids[screen].add(rid)
            aeid=str(int(Decimal(row["aeid"])))
            if aeid not in catalog:raise ValueError("Result missing from complete endpoint catalog")
            results_with_target+=aeid in annotated
            status_counts[screen+"_"+expected+"_results"]+=1
            cid=row["chemical_id"]
            if not cid:
                unresolved+=1;continue
            if cid not in ids:raise ValueError("Dangling chemical reference")
            counts[cid,screen]+=1
            endpoints[cid].add(str(int(Decimal(row["aeid"]))))
            present.add(cid)
            if expected in ("active","inactive"):interpretable.add(cid)
    calculated = {cid:(counts[cid,"mc"],counts[cid,"sc"],len(endpoints[cid])) for cid in present}
    if calculated != summaries:
        raise ValueError("Master counts disagree with independent native-ID aggregation")
    m=manifest["metrics"]
    calculated_metrics={
        "baseline_chemicals":len(original_ids), "expanded_chemicals":len(ids),
        "baseline_with_tox21":len(tox21),"baseline_newly_interpretable":len((interpretable & original_ids)-tox21),
        "fixed_catalog_covered_after":len(tox21 | (interpretable & original_ids)),
        "expanded_catalog_covered_after":len(tox21 | interpretable),
        "processed_results":sum(map(len,result_ids.values())),"unresolved_processed_results":unresolved,
        "mc_results":len(result_ids["mc"]),"sc_results":len(result_ids["sc"]),
        "assay_endpoints_documented":len(documented),
        "assay_endpoints_missing_annotations":len(catalog)-len(documented),
        "assay_target_links":len(links),"results_with_target_annotation":results_with_target,
        "toxcast_interpretable_chemicals":len(interpretable),"toxcast_source_result_chemicals":len(present),
        **status_counts}
    for key,value in calculated_metrics.items():
        if m[key]!=value:raise ValueError(f"Coverage reconciliation failed: {key}")
    for table in manifest["tables"]:
        if table["output"] in ("toxcast_mc_observations.csv","toxcast_sc_observations.csv"):
            screen=table["kind"]
            if raw_counts[screen]!=table["rows"]:raise ValueError("Source row count mismatch")
    return {"all_output_hashes_verified":True,"original_master_values_and_reviewed_keys_preserved":True,
            "native_processing_ids_reconciled":True,"replicates_do_not_inflate_coverage":True,
            "endpoint_catalog_and_target_links_verified":True,
            "independent_coverage_metrics":calculated_metrics}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot",type=Path,required=True)
    parser.add_argument("--reviewed-master",type=Path,required=True)
    parser.add_argument("--compare",type=Path)
    parser.add_argument("--report",type=Path,help="Optional JSON verification report outside the immutable snapshot")
    args=parser.parse_args()
    result=validate(args.snapshot,args.reviewed_master)
    if args.compare:
        first=json.loads((args.snapshot/"toxcast_import.json").read_text())
        second=json.loads((args.compare/"toxcast_import.json").read_text())
        if first!=second:raise ValueError("Rerun manifests differ")
        for name,digest in first["output_sha256"].items():
            if fingerprint(args.compare/name)!=digest:raise ValueError(f"Rerun differs: {name}")
        result["full_rerun_byte_identical"]=True
    result["release"]=json.loads((args.snapshot/"toxcast_import.json").read_text())["release"]
    result["python_version"]=platform.python_version()
    result["implementation_sha256"]={name:fingerprint(Path(__file__).with_name(name)) for name in (
        "download_toxcast.py","import_toxcast.py","verify_toxcast.py","toxcast_sources.json")}
    if args.report:
        report=args.report.resolve()
        if report.is_relative_to(args.snapshot.resolve()) or (args.compare and report.is_relative_to(args.compare.resolve())):
            raise ValueError("Verification report must be outside immutable snapshots")
        args.report.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
