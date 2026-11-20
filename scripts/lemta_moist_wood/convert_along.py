"""Convert the LEMTA along-grain workbooks to CSV, one file per sheet.

The workbooks (L. Terrei, LEMTA, by email, 2026-11-20) hold the along-grain
runs of Diarra et al., Fire Saf. J. 165 (2026) 104951, and the slice profiles
of its Fig. 11. Each sheet is written as stored: every cell's text unchanged,
rows in order, an empty cell where the sheet has none, so that the output is
the same on every run and can be pinned (data_files_along.toml). Standard
library only.

    python3 -I scripts/lemta_moist_wood/convert_along.py [DATA_DIR]

DATA_DIR (default $LEMTA_ALONG_DATA, else ~/Documents/data/lemta-parallel)
holds temperature_mass_parallel/*.xlsx; the CSVs go to DATA_DIR/csv/ as
<workbook>__<sheet>.csv.
"""

import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006/relationships"


def column_index(ref):
    """Zero-based column of a cell reference such as 'AB12'."""
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def shared_strings(book):
    if "xl/sharedStrings.xml" not in book.namelist():
        return []
    root = ET.fromstring(book.read("xl/sharedStrings.xml"))
    return [
        "".join(t.text or "" for t in si.iter(f"{{{MAIN}}}t"))
        for si in root.findall(f"{{{MAIN}}}si")
    ]


def cell_text(c, strings):
    kind = c.get("t")
    v = c.find(f"{{{MAIN}}}v")
    if kind == "s" and v is not None:
        return strings[int(v.text)]
    if kind == "inlineStr":
        return "".join(t.text or "" for t in c.iter(f"{{{MAIN}}}t"))
    return v.text if v is not None and v.text is not None else ""


def sheets(path):
    """(sheet name, rows) for every sheet of the workbook, in workbook order."""
    book = zipfile.ZipFile(path)
    strings = shared_strings(book)
    rels = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels.findall(f"{{{PKG}}}Relationship")}
    wb = ET.fromstring(book.read("xl/workbook.xml"))
    for sh in wb.find(f"{{{MAIN}}}sheets").findall(f"{{{MAIN}}}sheet"):
        t = target[sh.get(f"{{{REL}}}id")]
        part = t.lstrip("/") if t.startswith("/") else "xl/" + t
        rows = []
        for row in ET.fromstring(book.read(part)).iter(f"{{{MAIN}}}row"):
            cells = {
                column_index(c.get("r")): cell_text(c, strings)
                for c in row.findall(f"{{{MAIN}}}c")
            }
            index = int(row.get("r")) - 1
            rows.extend([] for _ in range(index - len(rows)))
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
        yield sh.get("name").strip(), rows


def main():
    root = (
        sys.argv[1]
        if len(sys.argv) > 1
        else os.environ.get(
            "LEMTA_ALONG_DATA",
            os.path.join(os.path.expanduser("~"), "Documents", "data", "lemta-parallel"),
        )
    )
    src = os.path.join(root, "temperature_mass_parallel")
    out = os.path.join(root, "csv")
    os.makedirs(out, exist_ok=True)
    for name in sorted(f for f in os.listdir(src) if f.endswith(".xlsx")):
        stem = os.path.splitext(name)[0]
        for sheet, rows in sheets(os.path.join(src, name)):
            path = os.path.join(out, f"{stem}__{sheet}.csv")
            with open(path, "w", newline="") as f:
                csv.writer(f, lineterminator="\n").writerows(rows)
            print(f"{path}: {len(rows)} rows")


if __name__ == "__main__":
    main()
