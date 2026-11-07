"""Scientific-grain, provenance, publication, and download recovery checks."""
import csv
import hashlib
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch
from xml.sax.saxutils import escape

from download_toxcast import download, fingerprint
from import_toxcast import Importer, activity_status, cytotox_value_status, init_db, potency_unit, run, workbook_tables
from verify_toxcast import validate


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def xlsx(sheets):
    """Tiny source fixtures with deliberately wrong A1:A1 dimensions."""
    out = io.BytesIO()
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    with zipfile.ZipFile(out, "w") as z:
        z.writestr("[Content_Types].xml", '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="xml" ContentType="application/xml"/><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' + ''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1)) + '</Types>')
        z.writestr("_rels/.rels", f'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr("xl/workbook.xml", f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>' + ''.join(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>' for i,(name,_) in enumerate(sheets,1)) + '</sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' + ''.join(f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheets)+1)) + '</Relationships>')
        for i, (_, rows) in enumerate(sheets,1):
            xml = []
            for n,row in enumerate(rows,1):
                cells = []
                for j,val in enumerate(row):
                    column, k = "", j+1
                    while k:
                        k, rem = divmod(k-1,26); column = chr(65+rem)+column
                    cells.append(f'<c r="{column}{n}" t="inlineStr"><is><t>{escape(str(val))}</t></is></c>')
                xml.append(f'<row r="{n}">' + ''.join(cells) + '</row>')
            z.writestr(f"xl/worksheets/sheet{i}.xml", f'<worksheet xmlns="{ns}"><dimension ref="A1:A1"/><sheetData>' + ''.join(xml) + '</sheetData></worksheet>')
    return out.getvalue()


class InterpretationTests(unittest.TestCase):
    def test_calls(self):
        cases = [("0", "inactive"), ("0.899999", "inactive"), ("0.90", "active"),
                 ("1", "active"), ("-0.98", "negative_direction"), ("NA", "missing"),
                 ("1.1", "invalid"), ("Infinity", "invalid")]
        for raw, expected in cases:
            self.assertEqual(activity_status({"hitc": raw, "modl": "hill"}, "mc"), expected)
        self.assertEqual(activity_status({"hitc": "0", "modl": "none"}, "mc"), "unfit")
        self.assertEqual(activity_status({"hitc": "0.9"}, "sc"), "invalid")
        self.assertEqual(activity_status({"hitc": "0"}, "sc"), "inactive")

    def test_units(self):
        self.assertEqual(potency_unit({"conc_unit": "uM"}, "mc"), "uM")
        self.assertEqual(potency_unit({"conc_unit": "NA"}, "mc"), "")
        self.assertEqual(potency_unit({"conc_unit": "nM"}, "mc"), "nM")
        self.assertEqual(potency_unit({"conc_unit": "CF"}, "mc"), "CF")

    def test_bad_excel_dimensions(self):
        b = xlsx([("data", [("one","two"), ("a","b"), ("c","d")])])
        for _, fields, rows in workbook_tables(io.BytesIO(b)):
            self.assertEqual(fields, ("one","two"))
            self.assertEqual(len(list(rows)),2)

    def test_cytotoxicity_default_is_not_a_measurement(self):
        self.assertEqual(cytotox_value_status({"cytotox_median_um":"1000","ntested":"80","nhit":"1"}),"default_1000_um")
        self.assertEqual(cytotox_value_status({"cytotox_median_um":"8","ntested":"80","nhit":"10"}),"estimated_from_assays")

    def test_conflicting_processing_id_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);db=init_db(p/"x.sqlite");imp=Importer(db,p,{})
            one={"m5id":"1","aeid":"1","spid":"a","dsstox_substance_id":"DTXSID1000001","hitc":"1","modl":"hill"}
            two={**one,"hitc":"0"}
            try:
                with self.assertRaisesRegex(ValueError,"Conflicting values"):
                    imp.process("input.csv","",tuple(one),[(2,one),(3,two)],"mc.csv","mc","a"*64)
            finally:imp.close();db.close()

    def test_identity_ambiguity_new_substances_and_no_name_matching(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp); db = init_db(p/"x.sqlite")
            db.executemany("INSERT INTO chemicals VALUES (?,1,0)", [("a",),("b",)])
            db.executemany("INSERT INTO sids VALUES (?,?)", [("DTXSID1000001","a"),("DTXSID1000001","b")])
            imp = Importer(db,p,{})
            prov = dict(source_file="x",source_sheet="",source_row_number=2,source_sha256="abc")
            self.assertEqual(imp.resolve({"dsstox_substance_id":"DTXSID1000001"},prov),("","ambiguous_dtxsid"))
            first = imp.resolve({"dsstox_substance_id":"DTXSID1000002"},prov)
            self.assertEqual(first[1],"new_substance")
            self.assertEqual(first,imp.resolve({"dsstox_substance_id":"DTXSID1000002"},prov))
            self.assertEqual(imp.resolve({"casn":"123-45-6","chnm":"same name"},prov)[0],"")
            imp.close(); db.close()


class DownloadTests(unittest.TestCase):
    class Response(io.BytesIO):
        def __init__(self, data, status=200, headers=None):
            super().__init__(data); self.status=status; self.headers=headers or {}

    def test_resume_and_ignore_range(self):
        data = b"public EPA fixture"
        item = {"filename":"x.dat","url":"https://example.test/x","size":len(data),"sha256":hashlib.sha256(data).hexdigest()}
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp); (p/"x.dat.part").write_bytes(data[:5])
            with patch("urllib.request.urlopen",return_value=self.Response(data[5:],206,{"Content-Range":f"bytes 5-{len(data)-1}/{len(data)}"})) as call:
                download(item,p)
                self.assertEqual(call.call_args[0][0].get_header("Range"),"bytes=5-")
            self.assertEqual((p/"x.dat").read_bytes(),data)
            (p/"x.dat").unlink(); (p/"x.dat.part").write_bytes(data[:5])
            with patch("urllib.request.urlopen",return_value=self.Response(data,200)):
                download(item,p)
            self.assertEqual((p/"x.dat").read_bytes(),data)

    def test_corrupt_download_not_published(self):
        with tempfile.TemporaryDirectory() as temp:
            item={"filename":"x.dat","url":"https://example.test/x","size":3,"sha256":hashlib.sha256(b"yes").hexdigest()}
            with patch("urllib.request.urlopen",return_value=self.Response(b"bad")):
                with self.assertRaises(ValueError): download(item,Path(temp))
            self.assertFalse((Path(temp)/"x.dat").exists())

    def test_interrupted_response_resumes(self):
        data=b"complete release"
        item={"filename":"x.dat","url":"https://example.test/x","size":len(data),"sha256":hashlib.sha256(data).hexdigest()}
        replies=[self.Response(data[:3]),self.Response(data[3:],206,{"Content-Range":f"bytes 3-{len(data)-1}/{len(data)}"})]
        with tempfile.TemporaryDirectory() as temp:
            with patch("urllib.request.urlopen",side_effect=replies),patch("download_toxcast.time.sleep"):
                download(item,Path(temp),attempts=2)
            self.assertEqual((Path(temp)/"x.dat").read_bytes(),data)


class EndToEndTests(unittest.TestCase):
    def fixture(self,p):
        base,review,source = (p/n for n in ("base","review","source"))
        for d in (base,review,source):d.mkdir()
        fields=("chemical_id","standard_inchikey","representative_smiles","identity_status","NR-AR")
        rows=[dict(zip(fields,("baseline-a","KEY-A","C","matched_by_structure","1"))),dict(zip(fields,("baseline-b","","O","dsstox_unresolved","")))]
        write_csv(base/"master_chemicals.csv",fields,rows)
        write_csv(base/"dsstox_records.csv",("chemical_id","DTXSID","source_sha256"),[
            {"chemical_id":"baseline-a","DTXSID":"DTXSID1000001","source_sha256":"release"},
            {"chemical_id":"baseline-b","DTXSID":"DTXSID1000002","source_sha256":"release"}])
        write_csv(base/"tox21_observations.csv",("chemical_id","endpoint","result"),[{"chemical_id":"baseline-a","endpoint":"NR-AR","result":"1"}])
        (base/"dsstox_import.json").write_text(json.dumps({"source_sha256":"release","linked_rows":2}))
        write_csv(review/"key_updates.csv",("chemical_id",),[{"chemical_id":"baseline-b"}])
        rows[1].update(standard_inchikey="KEY-B",identity_status="key_derived_from_source_smiles")
        write_csv(review/"master_chemicals_reviewed.csv",fields,rows)
        audit={"dsstox_release_sha256":"release","master_rows":2,"derived_keys_applied":1,
               "input_sha256":{f.name:fingerprint(f) for f in base.iterdir()},
               "output_sha256":{f.name:fingerprint(f) for f in review.iterdir()}}
        (review/"key_review_manifest.json").write_text(json.dumps(audit))
        annotation=xlsx([("annotations_combined",[("aeid","assay_component_endpoint_name","assay_source_name","data_usability"),(1,"test_assay","Tox21",1)])])
        targets=xlsx([("Sheet 1",[("aeid","target_type","target_id","official_symbol"),(1,"entrez_gene_id",101,"GENE1"),(1,"entrez_gene_id",102,"GENE2")])])
        cyto=xlsx([("Sheet 1",[("dsstox_substance_id","cytotox_median_um","ntested","nhit"),("DTXSID1000001",1000,1,0)])])
        qc=xlsx([("Sheet 1",[("dsstox_substance_id","pass_or_caution"),("DTXSID1000002","caution")])])
        works={"assay_annotations":annotation,"assay_target_mappings":targets,"cytotox":cyto,"analytical_qc":qc}
        for prefix,data in works.items():(source/f"{prefix}_invitrodb_v4_3_AUG2025.xlsx").write_bytes(data)
        header=("m5id","m4id","aeid","spid","dsstox_substance_id","chnm","hitc","modl","conc_unit","ac50","mc6_flags")
        records=[(1,11,1,"sample1","DTXSID1000002","B",0.90,"hill","uM",2,"caution"),
                 (2,12,1,"sample2","DTXSID1000002","B",0,"cnst","uM","NA",""),
                 (3,13,1,"sample3","DTXSID1000003","C",1,"hill","uM",1,""),
                 (4,14,1,"","NA","",0,"none","NA","NA","")]
        buf=io.StringIO();w=csv.writer(buf);w.writerow(header);w.writerows(records)
        sc=xlsx([("sc1",[("aeid","spid","dsstox_substance_id","s1id","resp"),(1,"sample1","DTXSID1000002",1,4),(1,"sample1","DTXSID1000002",2,5)]),
                 ("sc2",[("aeid","spid","dsstox_substance_id","s2id","hitc"),(1,"sample1","DTXSID1000002",1,1),(1,"sample1","DTXSID1000002",1,1)])])
        with zipfile.ZipFile(source/"INVITRODB_SUMMARY.zip","w") as z:
            for prefix,data in works.items():z.writestr(f"{prefix}_invitrodb_v4_3_AUG2024.xlsx",data)
            z.writestr("mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv",buf.getvalue())
            z.writestr("sc1_sc2_invitrodb_v4_3_AUG2024.xlsx",sc)
            for name in ("methods","assay_quality","mc7"):
                z.writestr(f"{name}_invitrodb_v4_3_AUG2024.xlsx",xlsx([("Sheet 1",[("description",),("fixture",)])]))
        lock={"release":"invitrodb_v4.3","files":[{"filename":f.name,"sha256":fingerprint(f),"size":f.stat().st_size,"url":"https://example.test/"+f.name} for f in sorted(source.iterdir())]}
        lp=p/"lock.json";lp.write_text(json.dumps(lock))
        return base,review,source,lp

    def test_full_import_grain_preservation_and_determinism(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);base,review,source,lock=self.fixture(p)
            first=run(base,review,source,p/"first",lock)
            second=run(base,review,source,p/"second",lock)
            self.assertEqual(first,second)
            self.assertTrue(validate(p/"first",review/"master_chemicals_reviewed.csv")["native_processing_ids_reconciled"])
            m=first["metrics"]
            self.assertEqual(m["baseline_newly_interpretable"],1)
            self.assertEqual(m["new_entries_with_interpretable_results"],1)
            self.assertEqual(m["mc_results"],4)
            self.assertEqual(m["sc_results"],1)
            self.assertEqual(m["sample_disagreement_groups"],1)
            self.assertEqual(m["assay_target_links"],2)
            self.assertEqual(m["results_with_target_annotation"],5)
            with (p/"first/master_chemicals.csv").open() as f:rows=list(csv.DictReader(f))
            b=next(r for r in rows if r["chemical_id"]=="baseline-b")
            self.assertEqual((b["standard_inchikey"],b["toxcast_mc_result_count"],b["toxcast_sc_result_count"],b["toxcast_endpoint_count"]),("KEY-B","2","1","1"))
            with self.assertRaises(FileExistsError):run(base,review,source,p/"first",lock)

    def test_changed_input_and_failure_leave_no_snapshot(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);base,review,source,lock=self.fixture(p)
            (base/"tox21_observations.csv").write_text("changed")
            with self.assertRaises(ValueError):run(base,review,source,p/"output",lock)
            self.assertFalse((p/"output").exists())
            self.assertFalse(list(p.glob(".toxcast-*")))

    def test_missing_released_annotation_keeps_results_and_marks_review(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);base,review,source,lock=self.fixture(p)
            path=source/"INVITRODB_SUMMARY.zip"
            with zipfile.ZipFile(path) as z:members={n:z.read(n) for n in z.namelist()}
            name=next(n for n in members if n.startswith("mc5-6"))
            reader=csv.DictReader(io.StringIO(members[name].decode())); fields=reader.fieldnames; records=list(reader)
            records[0]["aeid"]="2"
            out=io.StringIO();writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader();writer.writerows(records)
            members[name]=out.getvalue().encode()
            with zipfile.ZipFile(path,"w") as z:
                for n,b in members.items():z.writestr(n,b)
            data=json.loads(lock.read_text())
            for item in data["files"]:
                if item["filename"]==path.name:item.update(size=path.stat().st_size,sha256=fingerprint(path))
            lock.write_text(json.dumps(data))
            result=run(base,review,source,p/"output",lock)
            self.assertEqual(result["metrics"]["assay_endpoints_missing_annotations"],1)
            self.assertEqual(result["metrics"]["processed_results"],5)
            self.assertTrue(validate(p/"output",review/"master_chemicals_reviewed.csv")["native_processing_ids_reconciled"])

    def test_input_changed_during_import_not_published(self):
        import import_toxcast
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);base,review,source,lock=self.fixture(p)
            original=import_toxcast.import_sources
            def change(*args):
                inventory=original(*args)
                (base/"tox21_observations.csv").write_text("changed during import")
                return inventory
            with patch("import_toxcast.import_sources",side_effect=change):
                with self.assertRaises(ValueError):run(base,review,source,p/"output",lock)
            self.assertFalse((p/"output").exists())


if __name__=="__main__":
    unittest.main()
