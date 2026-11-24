"""P1-1b targeted checks (spec 17 P1-1b): CSV temporary files are always removed.

Every test points ``tempfile.tempdir`` at its own empty directory and asserts
the directory is empty afterwards.  Successful and Range responses are checked
against the golden C01 body and headers; error cases check that the original
error is what surfaces.
"""
import ast
import asyncio
import hashlib
import json
import logging
import os
import tempfile
from pathlib import Path

import polars as pl
import pytest

from src import responses
from src.responses import TempFileResponse

API = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((API / "tests/golden/manifest.json").read_text())
C01 = MANIFEST["cases"]["C01"]
C01_URL = C01["url"]
C01_HEADERS = C01["headers"]  # content-type, content-disposition, accept-ranges, content-length
SIZE = C01["size"]


@pytest.fixture()
def tmpdir_only(tmp_path, monkeypatch):
    d = tmp_path / "tmp"
    d.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(d))
    yield d
    assert sorted(os.listdir(d)) == [], "temporary file left behind"


def _golden_headers(r):
    return {k: r.headers.get(k) for k in C01_HEADERS}


# ------------------------------------------------------------- via the endpoint


def test_200_full_download(client, tmpdir_only):
    r = client.get(C01_URL)
    assert r.status_code == 200
    assert hashlib.sha256(r.content).hexdigest() == C01["sha256"]
    assert _golden_headers(r) == C01_HEADERS


def test_206_range_0_9(client, tmpdir_only):
    full = client.get(C01_URL).content
    r = client.get(C01_URL, headers={"Range": "bytes=0-9"})
    assert r.status_code == 206
    assert r.content == full[:10]
    assert r.headers["content-range"] == f"bytes 0-9/{SIZE}"
    assert r.headers["content-length"] == "10"
    assert r.headers["content-type"] == C01_HEADERS["content-type"]
    assert r.headers["content-disposition"] == C01_HEADERS["content-disposition"]
    assert r.headers["accept-ranges"] == "bytes"


def test_416_range_not_satisfiable(client, tmpdir_only):
    r = client.get(C01_URL, headers={"Range": "bytes=999999999-"})
    assert r.status_code == 416
    assert r.headers["content-range"] == f"*/{SIZE}"


def test_400_malformed_range(client, tmpdir_only):
    r = client.get(C01_URL, headers={"Range": "items=0-1"})
    assert r.status_code == 400
    assert r.text == "Only support bytes range"


@pytest.mark.parametrize("exc,status,detail", [
    (RuntimeError("disk full"), 500, "Internal server error. Please try it later or inform admin"),
    (ValueError("bad value from writer"), 400, "bad value from writer"),  # the original error, unchanged
])
def test_write_csv_failure_before_response(client, tmpdir_only, monkeypatch, exc, status, detail):
    seen = {}

    def failing_write_csv(self, path, *a, **k):
        seen["path"] = path
        Path(path).write_text("partial")  # the file exists and has content when it fails
        raise exc

    monkeypatch.setattr(pl.DataFrame, "write_csv", failing_write_csv)
    r = client.get(C01_URL)
    assert r.status_code == status
    assert r.json() == {"detail": detail}
    assert Path(seen["path"]).parent == tmpdir_only and not Path(seen["path"]).exists()


def test_no_background_unlink_task():
    src = (API / "mhw_app.py").read_text()
    names = {n.id for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Name)}
    assert "BackgroundTask" not in names and "BackgroundTasks" not in names


# ------------------------------------------------------------ direct ASGI calls


def _scope(range_header=None, method="GET"):
    headers = [(b"range", range_header.encode())] if range_header else []
    return {"type": "http", "method": method, "path": "/", "headers": headers}


async def _receive():
    return {"type": "http.request", "body": b"", "more_body": False}


def _make_file(d: Path, data: bytes = b"0123456789abcdef") -> str:
    f = tempfile.NamedTemporaryFile(delete=False)
    with f:
        f.write(data)
    assert Path(f.name).parent == d
    return f.name


class BodySendError(Exception):
    pass


@pytest.mark.parametrize("range_header", [None, "bytes=0-9"])
def test_send_error_during_body_removes_file_and_keeps_error(tmpdir_only, range_header):
    path = _make_file(tmpdir_only)
    err = BodySendError("client went away")
    sent = []

    async def send(msg):
        sent.append(msg["type"])
        if msg["type"] == "http.response.body":
            raise err

    resp = TempFileResponse(path, media_type="text/csv", filename="x.csv")
    with pytest.raises(BodySendError) as ei:
        asyncio.run(resp(_scope(range_header), _receive, send))
    assert ei.value is err  # the original exception object, not replaced by cleanup
    assert sent == ["http.response.start", "http.response.body"]
    assert resp.background is None


def test_direct_call_removes_file_on_success(tmpdir_only):
    path = _make_file(tmpdir_only)
    body = []

    async def send(msg):
        if msg["type"] == "http.response.body":
            body.append(msg.get("body", b""))

    asyncio.run(TempFileResponse(path)(_scope(), _receive, send))
    assert b"".join(body) == b"0123456789abcdef"


# ---------------------------------------------------------- cleanup error policy


def test_missing_file_is_ignored_and_original_error_kept(tmpdir_only):
    path = _make_file(tmpdir_only)
    os.unlink(path)  # already gone before the response runs

    async def send(msg):
        pass

    with pytest.raises(RuntimeError, match="does not exist"):  # starlette's own error, not FileNotFoundError
        asyncio.run(TempFileResponse(path)(_scope(), _receive, send))


def test_cleanup_error_after_success_is_raised(tmpdir_only, monkeypatch):
    path = _make_file(tmpdir_only)
    real_unlink = os.unlink

    def deny(p):
        raise PermissionError(13, "denied", p)

    monkeypatch.setattr(responses.os, "unlink", deny)

    async def send(msg):
        pass

    try:
        with pytest.raises(PermissionError):
            asyncio.run(TempFileResponse(path)(_scope(), _receive, send))
    finally:
        real_unlink(path)


def test_cleanup_error_during_send_error_is_logged_not_replacing(tmpdir_only, monkeypatch, caplog):
    path = _make_file(tmpdir_only)
    real_unlink = os.unlink
    monkeypatch.setattr(responses.os, "unlink", lambda p: (_ for _ in ()).throw(PermissionError(13, "denied", p)))
    err = BodySendError("client went away")

    async def send(msg):
        if msg["type"] == "http.response.body":
            raise err

    try:
        with caplog.at_level(logging.ERROR, logger="src.responses"):
            with pytest.raises(BodySendError) as ei:
                asyncio.run(TempFileResponse(path)(_scope(), _receive, send))
        assert ei.value is err
        assert any("could not remove temporary file" in r.getMessage() for r in caplog.records)
        assert any(isinstance(r.exc_info[1], PermissionError) for r in caplog.records if r.exc_info)
    finally:
        real_unlink(path)
