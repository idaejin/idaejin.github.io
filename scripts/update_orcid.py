#!/usr/bin/env python3
"""Sync publications from ORCID into Hugo Academic content + data/publications.json.

Usage:
  python3 scripts/update_orcid.py
  ORCID_TOKEN=... python3 scripts/update_orcid.py   # optional read-public token

Env:
  ORCID_ID     default 0000-0002-8995-8535
  ORCID_TOKEN  optional Bearer token for /read-public
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ORCID_ID = os.environ.get("ORCID_ID", "0000-0002-8995-8535")
ORCID_TOKEN = os.environ.get("ORCID_TOKEN", "").strip()
API = f"https://pub.orcid.org/v3.0/{ORCID_ID}"

ROOT = Path(__file__).resolve().parents[1]
PUB_DIR = ROOT / "content" / "publication"
JSON_PATH = ROOT / "data" / "publications.json"
OVERRIDES_PATH = ROOT / "data" / "orcid_overrides.json"

TYPE_MAP = {
    "journal-article": "2",
    "conference-paper": "1",
    "book-chapter": "6",
    "book": "5",
    "report": "4",
    "working-paper": "3",
    "preprint": "3",
    "dissertation-thesis": "7",
    "online-resource": "3",
    "other": "0",
}

# Prefer summaries that look more complete when ORCID groups duplicates.
SOURCE_RANK = {
    "crossref": 0,
    "scopus": 1,
    "pubmed": 2,
    "europe pubmed central": 3,
    "datacite": 4,
}


def headers() -> dict[str, str]:
    h = {"Accept": "application/vnd.orcid+json", "User-Agent": "idaejin.github.io-orcid-sync/1.0"}
    if ORCID_TOKEN:
        h["Authorization"] = f"Bearer {ORCID_TOKEN}"
    return h


def api_get(url: str, retries: int = 4) -> dict:
    last_err: Exception | None = None
    for attempt in range(retries):
        req = urllib.request.Request(url, headers=headers())
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code in {429, 500, 502, 503, 504}:
                time.sleep(1.5 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as e:
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"ORCID request failed for {url}: {last_err}")


def text(node, *keys, default=""):
    cur = node
    for k in keys:
        if cur is None:
            return default
        cur = cur.get(k) if isinstance(cur, dict) else None
    if isinstance(cur, dict):
        return (cur.get("value") or default) or default
    return cur or default


def normalize_doi(doi: str) -> str:
    doi = (doi or "").strip()
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi, flags=re.I)
    return doi.lower().rstrip(".")


def slugify(title: str, year: str, put: int | str) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")[:60].strip("-")
    return f"{year}-{s}-{put}"


def is_self(name: str) -> bool:
    n = re.sub(r"[^a-z]", "", name.lower())
    return ("daejin" in n and "lee" in n) or n in {"daejinlee", "leedaejin"}


def yaml_escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def pick_summary(summaries: list[dict]) -> dict:
    def score(s: dict) -> tuple:
        src = ((s.get("source") or {}).get("source-name") or {}).get("value", "")
        src_rank = SOURCE_RANK.get(src.lower(), 50)
        has_journal = 0 if text(s, "journal-title") else 1
        has_doi = 0 if extract_doi(s) else 1
        return (src_rank, has_journal, has_doi, -(s.get("put-code") or 0))

    return sorted(summaries, key=score)[0]


def extract_doi(work: dict) -> str:
    ids = ((work.get("external-ids") or {}).get("external-id")) or []
    for e in ids:
        if (e.get("external-id-type") or "").lower() == "doi":
            return normalize_doi(e.get("external-id-value") or "")
    return ""


def extract_url(work: dict) -> str:
    url = text(work, "url")
    if url:
        return url
    ids = ((work.get("external-ids") or {}).get("external-id")) or []
    for e in ids:
        if (e.get("external-id-type") or "").lower() == "uri":
            return e.get("external-id-value") or ""
    doi = extract_doi(work)
    return f"https://doi.org/{doi}" if doi else ""


def date_parts(work: dict) -> tuple[str, str, str]:
    pd = work.get("publication-date") or {}
    year = text(pd, "year") or "1900"
    month = text(pd, "month") or "01"
    day = text(pd, "day") or "01"
    return year, month.zfill(2), day.zfill(2)


def title_key(title: str, year: str) -> str:
    return f"{year}:" + re.sub(r"[^a-z0-9]+", "", title.lower())


def fetch_contributors(put_code: int) -> list[str]:
    detail = api_get(f"{API}/work/{put_code}")
    contribs = ((detail.get("contributors") or {}).get("contributor")) or []
    names = []
    for c in contribs:
        name = text(c, "credit-name")
        if name:
            names.append(name.strip())
    return names


def to_authors(names: list[str]) -> list[str]:
    authors: list[str] = []
    self_added = False
    for n in names:
        if is_self(n):
            if not self_added:
                authors.append("admin")
                self_added = True
        else:
            authors.append(n)
    if not authors:
        return ["admin"]
    if not self_added:
        authors.append("admin")
    return authors


def load_overrides() -> dict:
    if not OVERRIDES_PATH.exists():
        return {}
    return json.loads(OVERRIDES_PATH.read_text())


def apply_author_overrides(put_code: int | str, authors: list[str], overrides: dict) -> list[str]:
    cfg = (overrides.get("by_put_code") or {}).get(str(put_code)) or {}
    remove = {re.sub(r"[^a-z0-9]+", "", a.lower()) for a in (cfg.get("remove_authors") or [])}
    if remove:
        authors = [a for a in authors if re.sub(r"[^a-z0-9]+", "", a.lower()) not in remove]
    if cfg.get("authors"):
        authors = list(cfg["authors"])
    if not authors:
        authors = ["admin"]
    return authors


def collect_works() -> list[dict]:
    data = api_get(f"{API}/works")
    groups = data.get("group") or []
    candidates: list[dict] = []

    for g in groups:
        summaries = g.get("work-summary") or []
        if not summaries:
            continue
        s = pick_summary(summaries)
        title = text(s, "title", "title")
        if not title:
            continue
        year, month, day = date_parts(s)
        wtype = s.get("type") or "other"
        put = s.get("put-code")
        doi = extract_doi(s)
        journal = text(s, "journal-title")
        url = extract_url(s)
        src = ((s.get("source") or {}).get("source-name") or {}).get("value", "")
        candidates.append(
            {
                "put_code": put,
                "title": title,
                "year": year,
                "month": month,
                "day": day,
                "type": wtype,
                "doi": doi,
                "journal": journal,
                "url": url,
                "source": src,
            }
        )

    # Deduplicate by DOI, then by normalized title (ignore year mismatches).
    by_doi: dict[str, dict] = {}
    by_title: dict[str, dict] = {}

    def richer(a: dict, b: dict) -> dict:
        score_a = (1 if a["doi"] else 0, 1 if a["journal"] else 0, len(a["title"]), a["year"])
        score_b = (1 if b["doi"] else 0, 1 if b["journal"] else 0, len(b["title"]), b["year"])
        return a if score_a >= score_b else b

    def norm_title(title: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", title.lower())

    for c in candidates:
        if c["doi"]:
            prev = by_doi.get(c["doi"])
            by_doi[c["doi"]] = richer(prev, c) if prev else c
        else:
            key = norm_title(c["title"])
            prev = by_title.get(key)
            by_title[key] = richer(prev, c) if prev else c

    doi_titles = {norm_title(v["title"]) for v in by_doi.values()}
    unique = list(by_doi.values())
    for key, v in by_title.items():
        if key not in doi_titles:
            unique.append(v)

    unique.sort(key=lambda x: (x["year"], x["month"], x["day"], x["title"]), reverse=True)

    print(f"ORCID groups: {len(groups)}; after dedupe: {len(unique)}")
    overrides = load_overrides()
    for i, w in enumerate(unique):
        names = fetch_contributors(w["put_code"])
        w["authors"] = apply_author_overrides(w["put_code"], to_authors(names), overrides)
        if (i + 1) % 10 == 0:
            print(f"  fetched authors {i + 1}/{len(unique)}")
        time.sleep(0.05)
    return unique


def write_json(works: list[dict]) -> None:
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "orcid": ORCID_ID,
        "source": f"{API}/works",
        "count": len(works),
        "publications": works,
    }
    JSON_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {JSON_PATH.relative_to(ROOT)}")


def write_hugo(works: list[dict]) -> None:
    PUB_DIR.mkdir(parents=True, exist_ok=True)
    for p in PUB_DIR.iterdir():
        if p.name == "_index.md":
            continue
        if p.is_dir():
            shutil.rmtree(p)

    for w in works:
        slug = slugify(w["title"], w["year"], w["put_code"])
        folder = PUB_DIR / slug
        folder.mkdir(parents=True, exist_ok=True)
        date = f'{w["year"]}-{w["month"]}-{w["day"]}'
        pub_type = TYPE_MAP.get(w["type"], "0")
        pub_line = f'*{w["journal"]}*' if w["journal"] else ""
        doi = w["doi"]
        auth_yaml = "\n".join(
            '  - admin' if a == "admin" else f'  - "{yaml_escape(a)}"' for a in w["authors"]
        )
        md = f'''---
title: "{yaml_escape(w["title"])}"
authors:
{auth_yaml}
date: "{date}"
publishDate: "{date}"
doi: "{doi}"
publication_types: ["{pub_type}"]
publication: "{yaml_escape(pub_line)}"
publication_short: ""
abstract: ""
featured: false
url_pdf: ""
url_code: ""
url_dataset: ""
url_poster: ""
url_project: ""
url_slides: ""
url_source: ""
url_video: ""
orcid_put_code: {w["put_code"]}
---
'''
        (folder / "index.md").write_text(md)

    index = PUB_DIR / "_index.md"
    index.write_text(
        f'''---
title: Publications
cms_exclude: true

# View.
#   1 = List
#   2 = Compact
#   3 = Card
#   4 = Citation
view: 4

filters:
  publication_type: -1
---

Synced automatically from [ORCID](https://orcid.org/{ORCID_ID}) (weekly). Also on [Google Scholar](https://scholar.google.com/citations?user=hsCRT4EAAAAJ).
'''
    )
    print(f"Wrote {len(works)} Hugo publication pages")


def main() -> int:
    works = collect_works()
    write_json(works)
    write_hugo(works)
    return 0


if __name__ == "__main__":
    sys.exit(main())
