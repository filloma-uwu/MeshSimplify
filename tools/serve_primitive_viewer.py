"""Serve primitive-analysis output and its Three.js viewer."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit


PROJECT_DIRECTORY = Path(__file__).resolve().parent.parent


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, directory=str(PROJECT_DIRECTORY), **kwargs)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def send_json(self, status: int, value: dict[str, object]) -> None:
        payload = json.dumps(value, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:  # noqa: N802 - HTTP handler API
        if self.path == "/":
            manifest = self.server.manifest_path.relative_to(
                PROJECT_DIRECTORY).as_posix()
            location = (
                "/viewer/primitive_analysis.html?manifest=/" +
                quote(manifest, safe="/")
            )
            self.send_response(302)
            self.send_header("Location", location)
            self.end_headers()
            return
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802 - HTTP handler API
        if self.path != "/api/regenerate":
            self.send_json(404, {"error": "unknown API endpoint"})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 64 * 1024:
                raise ValueError("invalid request size")
            request = json.loads(self.rfile.read(length).decode("utf-8"))
            model_id = str(request["model_id"])
            maximum_error = float(request["maximum_error"])
            if not math.isfinite(maximum_error) or maximum_error < 0.0:
                raise ValueError("maximum_error must be finite and non-negative")

            server = self.server
            requested_manifest = request.get("manifest")
            if requested_manifest is None:
                manifest_path = server.manifest_path
            else:
                manifest_url_path = unquote(
                    urlsplit(str(requested_manifest)).path).lstrip("/")
                manifest_path = (PROJECT_DIRECTORY / manifest_url_path).resolve()
                manifest_path.relative_to(PROJECT_DIRECTORY)
                if not manifest_path.is_file():
                    raise ValueError(f"manifest does not exist: {manifest_url_path}")
            manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
            model = next(
                (item for item in manifest_data.get("models", [])
                 if str(item.get("id")) == model_id),
                None,
            )
            if model is None:
                raise ValueError(f"model {model_id} is not present in the served manifest")
            metadata_path = (manifest_path.parent / model["metadata"]).resolve()
            metadata_path.relative_to(PROJECT_DIRECTORY)
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            source_path = (metadata_path.parent / metadata["source"]).resolve()
            source_path.relative_to(PROJECT_DIRECTORY)
            if not source_path.is_file():
                raise ValueError(f"source OBJ does not exist: {source_path}")

            error_key = format(maximum_error, ".12g").replace("-", "m").replace(".", "p")
            executable_version = server.executable.stat().st_mtime_ns
            output = (PROJECT_DIRECTORY / "outputs" / "interactive_max_error" /
                      f"binary_{executable_version}" / f"model_{model_id}" /
                      f"error_{error_key}").resolve()
            output.relative_to(PROJECT_DIRECTORY)
            generated_metadata = output / "model.json"
            cached = generated_metadata.is_file()
            with server.generation_lock:
                if not generated_metadata.is_file():
                    cached = False
                    output.mkdir(parents=True, exist_ok=True)
                    command = [
                        str(server.executable),
                        "--input", str(source_path),
                        "--output-dir", str(output),
                        "--maximum-open-error-distance", format(maximum_error, ".17g"),
                    ]
                    completed = subprocess.run(
                        command, cwd=PROJECT_DIRECTORY, capture_output=True,
                        text=True, timeout=30 * 60, check=False,
                    )
                    (output / "analysis.stdout.txt").write_text(
                        completed.stdout, encoding="utf-8")
                    (output / "analysis.stderr.txt").write_text(
                        completed.stderr, encoding="utf-8")
                    if completed.returncode != 0 or not generated_metadata.is_file():
                        raise RuntimeError(
                            completed.stderr.strip() or
                            f"analyzer exited with code {completed.returncode}")
            relative = generated_metadata.relative_to(PROJECT_DIRECTORY).as_posix()
            self.send_json(200, {"metadata": "/" + relative, "cached": cached})
        except (ValueError, KeyError, json.JSONDecodeError) as error:
            self.send_json(400, {"error": str(error)})
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
            self.send_json(500, {"error": str(error)})


class Server(ThreadingHTTPServer):
    # On Windows SO_REUSEADDR can let a stale viewer keep receiving requests
    # after a new process is launched on the same port. Fail fast instead.
    allow_reuse_address = False

    def __init__(self, address: tuple[str, int], handler: type[Handler],
                 manifest_path: Path, executable: Path) -> None:
        super().__init__(address, handler)
        self.manifest_path = manifest_path
        self.executable = executable
        self.generation_lock = threading.Lock()


def main() -> int:
    parser = argparse.ArgumentParser(description="Serve the primitive-analysis viewer")
    parser.add_argument("--port", type=int, default=8091)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("outputs/direct_coverage_uniform_closed_responsibility_full_v13/viewer_manifest.json"),
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=Path("build/Release/pqss-primitive-mesh-analyze.exe"),
    )
    args = parser.parse_args()
    manifest = args.manifest.resolve()
    try:
        relative_manifest = manifest.relative_to(PROJECT_DIRECTORY)
    except ValueError as error:
        raise SystemExit(f"manifest must be inside {PROJECT_DIRECTORY}: {manifest}") from error
    if not manifest.is_file():
        raise SystemExit(f"manifest does not exist: {manifest}")
    executable = args.executable.resolve()
    if not executable.is_file():
        raise SystemExit(f"analyzer executable does not exist: {executable}")
    manifest_url = "/" + quote(relative_manifest.as_posix(), safe="/")
    url = (
        f"http://127.0.0.1:{args.port}/viewer/primitive_analysis.html"
        f"?manifest={manifest_url}"
    )
    print(f"PQSS primitive-analysis viewer: {url}", flush=True)
    Server(("127.0.0.1", args.port), Handler, manifest, executable).serve_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
