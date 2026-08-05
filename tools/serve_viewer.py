"""Serve the original/proxy mesh comparison viewer and its model metadata."""

from __future__ import annotations

import csv
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import sys
from urllib.parse import urlparse


PROJECT_DIR = Path(__file__).resolve().parent.parent
GENERATION_REPORT = (
    PROJECT_DIR
    / "outputs"
    / "real_scene"
    / "conservative_outer_adaptive_depth8_strict"
    / "generation_report.json"
)
BVH_REPORT = (
    PROJECT_DIR
    / "reports"
    / "real_scene"
    / "adaptive_depth8_validated"
    / "model_bvh_stats.csv"
)


def _number(value: str) -> int | float | str:
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def model_metadata() -> dict[str, object]:
    generation = json.loads(GENERATION_REPORT.read_text(encoding="ascii"))
    with BVH_REPORT.open("r", encoding="utf-8-sig", newline="") as stream:
        bvh_by_id = {
            int(row["model_id"]): {key: _number(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
            if row["category"] == "base"
        }

    models = []
    for item in generation["models"]:
        model_id = int(Path(item["model"]).stem)
        models.append(
            {
                "id": model_id,
                "name": item["model"],
                "originalUrl": f"/test_data/real_scene/source_pool/{model_id}.obj",
                "proxyUrl": (
                    "/outputs/real_scene/conservative_outer_adaptive_depth8_strict/"
                    f"{model_id}.obj"
                ),
                "generation": item,
                "bvh": bvh_by_id[model_id],
            }
        )
    models.sort(key=lambda item: item["id"])
    return {
        "requestedMaxDepth": generation["options"]["max_pqss_bvh_depth"],
        "models": models,
    }


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, directory=str(PROJECT_DIR), **kwargs)

    def log_message(self, *_: object) -> None:
        pass

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self) -> None:
        if urlparse(self.path).path == "/api/models":
            try:
                body = json.dumps(model_metadata(), separators=(",", ":")).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
                body = json.dumps({"error": str(error)}).encode("utf-8")
                self.send_response(500)
                self.send_header("Content-Type", "application/json; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            return
        super().do_GET()


class Server(ThreadingHTTPServer):
    allow_reuse_address = True


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8090
    print(f"PQSSProxyMesh viewer: http://127.0.0.1:{port}/viewer/", flush=True)
    Server(("127.0.0.1", port), Handler).serve_forever()
