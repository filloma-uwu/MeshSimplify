from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil


OFFICIAL_SIMPLIFIED_MODEL_IDS = (2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)
BASE_MODEL_IDS = tuple(range(1, 24))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def object_id(token: str) -> int:
    name = Path(token).name
    if not name.endswith(".obj"):
        raise ValueError(f"invalid model token: {token}")
    return int(name[:-4])


def filter_queries(source: Path, output: Path) -> tuple[int, int]:
    selected = set(OFFICIAL_SIMPLIFIED_MODEL_IDS)
    total = 0
    retained = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="ascii") as input_stream, output.open(
        "w", encoding="ascii", newline="\n"
    ) as output_stream:
        for line_number, line in enumerate(input_stream, 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 26:
                raise ValueError(f"invalid query at {source}:{line_number}")
            first_id = object_id(fields[0])
            second_id = object_id(fields[13])
            total += 1
            if first_id in selected or second_id in selected:
                output_stream.write(line.rstrip("\r\n") + "\n")
                retained += 1
    return total, retained


def prepare_hybrid_pool(original_dir: Path, candidate_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(f"hybrid output directory is not empty: {output_dir}")
    selected = set(OFFICIAL_SIMPLIFIED_MODEL_IDS)
    for model_id in BASE_MODEL_IDS:
        source_dir = candidate_dir if model_id in selected else original_dir
        source = source_dir / f"{model_id}.obj"
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, output_dir / source.name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare the fixed real-scene scope matching the officially simplified models"
    )
    parser.add_argument("--query-source", type=Path, required=True)
    parser.add_argument("--query-output", type=Path, required=True)
    parser.add_argument("--original-model-dir", type=Path, required=True)
    parser.add_argument("--candidate-model-dir", type=Path, required=True)
    parser.add_argument("--hybrid-output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.query_output.exists():
        raise FileExistsError(f"query output already exists: {args.query_output}")
    total, retained = filter_queries(args.query_source, args.query_output)
    prepare_hybrid_pool(args.original_model_dir, args.candidate_model_dir, args.hybrid_output_dir)

    manifest = {
        "scope": "at_least_one_officially_simplified_model",
        "official_simplified_model_ids": list(OFFICIAL_SIMPLIFIED_MODEL_IDS),
        "query_source": str(args.query_source.resolve()),
        "query_source_sha256": file_sha256(args.query_source),
        "query_output": str(args.query_output.resolve()),
        "query_output_sha256": file_sha256(args.query_output),
        "unfiltered_query_count": total,
        "retained_query_count": retained,
        "excluded_query_count": total - retained,
        "hybrid_pool": str(args.hybrid_output_dir.resolve()),
        "candidate_models": list(OFFICIAL_SIMPLIFIED_MODEL_IDS),
        "original_passthrough_models": [
            model_id for model_id in BASE_MODEL_IDS
            if model_id not in OFFICIAL_SIMPLIFIED_MODEL_IDS
        ],
    }
    manifest_path = args.query_output.parent / "scope_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    print(f"unfiltered_queries={total}")
    print(f"retained_queries={retained}")
    print(f"excluded_queries={total - retained}")
    print(f"query_output={args.query_output}")
    print(f"hybrid_pool={args.hybrid_output_dir}")
    print(f"manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
