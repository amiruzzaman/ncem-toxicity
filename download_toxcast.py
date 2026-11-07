"""Download the pinned public EPA invitrodb v4.3 inputs, resuming partial files."""
import argparse
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

RELEASE = "invitrodb_v4.3"
DEFAULT_LOCK = Path(__file__).with_name("toxcast_sources.json")


def fingerprint(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_lock(path=DEFAULT_LOCK):
    data = json.loads(Path(path).read_text())
    if data.get("release") != RELEASE or not data.get("files"):
        raise ValueError("Unsupported or empty ToxCast source lock")
    names = set()
    for item in data["files"]:
        name = item["filename"]
        if Path(name).name != name or name in names:
            raise ValueError("Unsafe or duplicate source filename")
        if not re.fullmatch(r"[a-f0-9]{64}", item["sha256"]) or item["size"] <= 0:
            raise ValueError("Invalid source size or SHA-256")
        names.add(name)
    return data


def verify(path, item):
    if path.stat().st_size != item["size"] or fingerprint(path) != item["sha256"]:
        raise ValueError(f"Source size/hash mismatch: {path.name}")


def download(item, directory, attempts=5):
    target = directory / item["filename"]
    partial = target.with_name(target.name + ".part")
    if target.exists():
        verify(target, item)
        print(f"Verified cached {target.name}", flush=True)
        return target
    for attempt in range(attempts):
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > item["size"]:
            raise ValueError(f"Oversized partial download: {partial}")
        if offset == item["size"]:
            break
        request = urllib.request.Request(item["url"], headers={
            "User-Agent": "ResearchChem-ToxCast/1", "Accept-Encoding": "identity",
            **({"Range": f"bytes={offset}-"} if offset else {}),
        })
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if response.status == 206:
                    expected = f"bytes {offset}-"
                    if not response.headers.get("Content-Range", "").startswith(expected):
                        raise ValueError("Server returned the wrong resume range")
                    mode = "ab" if offset else "wb"
                elif response.status == 200:
                    # Some servers ignore Range. Restart rather than append duplicate bytes.
                    mode, offset = "wb", 0
                else:
                    raise ValueError(f"Unexpected HTTP status {response.status}")
                print(f"Downloading {target.name} from byte {offset:,}", flush=True)
                next_progress = offset + 256 * 1024**2
                with partial.open(mode) as stream:
                    while block := response.read(1024**2):
                        stream.write(block)
                        offset += len(block)
                        if offset > item["size"]:
                            raise ValueError("Source grew beyond its pinned size")
                        if offset >= next_progress:
                            print(f"  {offset:,}/{item['size']:,} bytes", flush=True)
                            next_progress = offset + 256 * 1024**2
            if offset != item["size"]:
                raise OSError("Incomplete HTTP response")
            break
        except (OSError, urllib.error.URLError) as error:
            if attempt + 1 == attempts:
                raise RuntimeError(f"Download interrupted; rerun to resume {partial}") from error
            print(f"Retry {attempt + 1}: {error}", flush=True)
            time.sleep(min(2 ** attempt, 8))
    verify(partial, item)
    if target.suffix == ".zip":
        with zipfile.ZipFile(partial) as archive:
            archive.infolist()  # Full pinned hash protects contents; selected members also get CRC checks.
    partial.rename(target)
    return target


def run(source_dir, lock_file=DEFAULT_LOCK):
    lock = load_lock(lock_file)
    source_dir.mkdir(parents=True, exist_ok=True)
    for item in lock["files"]:
        download(item, source_dir)
    manifest = {"release": RELEASE, "source_lock_sha256": fingerprint(lock_file), "files": lock["files"]}
    (source_dir / "download_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--source-lock", type=Path, default=DEFAULT_LOCK)
    args = parser.parse_args()
    run(args.source_dir, args.source_lock)
