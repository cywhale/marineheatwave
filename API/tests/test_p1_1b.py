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
import starlette

from src import responses
from src.responses import TempFileResponse

API = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((API / "tests/golden/manifest.json").read_text())
C01 = MANIFEST["cases"]["C01"]
C01_URL = C01["url"]
C01_HEADERS = C01["headers"]  # content-type, content-disposition, accept-ranges, content-length
SIZE = C01["size"]
# Starlette 1.x writes the unit in a 416 Content-Range (RFC 9110: "bytes */<size>");
# owner accepted this header change for P1-2 on 2026-11-25 (doc 29).  The body and
# every other response stay as recorded in golden.
STARLETTE_1 = int(starlette.__version__.split(".")[0]) >= 1


def unsatisfiable_range(size):
    return f"bytes */{size}" if STARLETTE_1 else f"*/{size}"


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
    # nothing else besides the validators FileResponse always sends
    assert set(r.headers) == {"content-range", "content-length", "content-type", "content-disposition",
                              "accept-ranges", "etag", "last-modified"}


def test_416_range_not_satisfiable(client, tmpdir_only):
    """The whole 416 response.  The only accepted change from Starlette 0.x to
    1.x is the unit in Content-Range; status, body and every other header must
    stay exactly as before (checked on both Starlette versions)."""
    r = client.get(C01_URL, headers={"Range": "bytes=999999999-"})
    assert r.status_code == 416
    assert r.content == b""
    assert dict(r.headers) == {
        "content-length": "0",
        "content-range": unsatisfiable_range(SIZE),
        "content-type": "text/plain; charset=utf-8",
    }


def test_400_malformed_range(client, tmpdir_only):
    r = client.get(C01_URL, headers={"Range": "items=0-1"})
    assert r.status_code == 400
    assert r.text == "Only support bytes range"
    assert dict(r.headers) == {"content-length": "24", "content-type": "text/plain; charset=utf-8"}


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


# The ASGI HTTP scope spec_version real servers send.  uvicorn (0.34.1 in the
# P1-1 venv, 0.54.0 in the P1-2 candidate) sends "2.3" from every HTTP protocol
# implementation (h11, httptools, zttp); test_uvicorn_http_spec_version pins it.
# Below 2.4 Starlette 1.x FileResponse listens for http.disconnect while it
# streams, which is the branch production takes.
SPEC_VERSION = "2.3"


def _scope(range_header=None, method="GET"):
    headers = [(b"range", range_header.encode())] if range_header else []
    return {"type": "http", "asgi": {"version": "3.0", "spec_version": SPEC_VERSION},
            "http_version": "1.1", "method": method, "path": "/", "headers": headers}


class _Client:
    """ASGI connection double that follows uvicorn's HTTP cycle.

    * The first ``receive()`` returns the (empty) request body at once, as
      uvicorn does for a GET whose body has been read.
    * Later ``receive()`` calls wait until the client disconnects or the
      response is complete (last body chunk sent, or sending raised), then
      return ``http.disconnect``.
    * ``disconnect_after_body_chunks=n`` simulates the client going away after
      the n-th body chunk: like uvicorn, later ``send()`` calls are dropped
      silently (no exception).

    A ``receive`` that always answered ``http.request`` at once made Starlette
    1.x FileResponse busy-loop forever; older Starlette never calls receive
    here.
    """

    def __init__(self, on_send=None, disconnect_after_body_chunks=None):
        self.on_send = on_send
        self.sent: list[dict] = []      # delivered to the client
        self.dropped: list[dict] = []   # sent after the client went away
        self.receives: list[str] = []
        self.done = asyncio.Event()
        self.disconnect_after = disconnect_after_body_chunks
        self.disconnected = False
        self._request_delivered = False

    async def receive(self):
        if not self._request_delivered:
            self._request_delivered = True
            self.receives.append("http.request")
            return {"type": "http.request", "body": b"", "more_body": False}
        await self.done.wait()
        self.receives.append("http.disconnect")
        return {"type": "http.disconnect"}

    async def send(self, msg):
        if self.disconnected:
            self.dropped.append(msg)
            return
        self.sent.append(msg)
        try:
            if self.on_send is not None:
                await self.on_send(msg)
        except BaseException:
            self.done.set()
            raise
        if msg["type"] == "http.response.body":
            n = sum(m["type"] == "http.response.body" for m in self.sent)
            if self.disconnect_after is not None and n >= self.disconnect_after:
                self.disconnected = True
                self.done.set()
            elif not msg.get("more_body", False):
                self.done.set()

    def body(self):
        return b"".join(m.get("body", b"") for m in self.sent if m["type"] == "http.response.body")


def _run(resp, scope, on_send=None, **client_kw):
    async def go():
        c = _Client(on_send, **client_kw)
        await resp(scope, c.receive, c.send)
        return c
    return asyncio.run(go())


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
        _run(resp, _scope(range_header), send)
    assert ei.value is err  # the original exception object, not replaced by cleanup
    assert sent == ["http.response.start", "http.response.body"]
    assert resp.background is None


def test_direct_call_removes_file_on_success(tmpdir_only):
    path = _make_file(tmpdir_only)
    body = []

    async def send(msg):
        if msg["type"] == "http.response.body":
            body.append(msg.get("body", b""))

    _run(TempFileResponse(path), _scope(), send)
    assert b"".join(body) == b"0123456789abcdef"


def test_direct_call_206_removes_file(tmpdir_only):
    path = _make_file(tmpdir_only)
    c = _run(TempFileResponse(path, media_type="text/csv"), _scope("bytes=0-9"))
    start = c.sent[0]
    assert start["status"] == 206
    assert b"".join(m.get("body", b"") for m in c.sent if m["type"] == "http.response.body") == b"0123456789"


def test_direct_call_416_removes_file_with_starlette_header(tmpdir_only):
    path = _make_file(tmpdir_only)  # 16 bytes
    c = _run(TempFileResponse(path, media_type="text/csv"), _scope("bytes=999-"))
    start = c.sent[0]
    assert start["status"] == 416
    headers = {k.decode().lower(): v.decode() for k, v in start["headers"]}
    assert headers == {"content-length": "0", "content-range": unsatisfiable_range(16),
                       "content-type": "text/plain; charset=utf-8"}
    assert c.body() == b""
    # tmpdir_only asserts that the file is gone


def test_uvicorn_http_spec_version():
    """The scope used above matches what the installed uvicorn sends (all HTTP
    protocol implementations); if uvicorn moves to 2.4, Starlette 1.x stops
    listening for disconnects and the tests below must be revisited."""
    import re

    import uvicorn.protocols.http as h

    found = {}
    for f in Path(h.__file__).parent.glob("*.py"):
        for v in re.findall(r'"spec_version":\s*"([0-9.]+)"', f.read_text()):
            found.setdefault(v, []).append(f.name)
    assert set(found) == {SPEC_VERSION}, found


BIG = 300_000  # more than four 64 KiB FileResponse chunks


@pytest.mark.parametrize("range_header", [None, "bytes=0-199999"])
def test_client_disconnect_mid_body_removes_file(tmpdir_only, range_header):
    """The client goes away after the first body chunk.  No exception reaches
    the server (uvicorn drops later sends), only the first chunk was delivered,
    and the temporary file is removed.  Starlette 1.x notices the disconnect
    (it was waiting in receive) and stops streaming; older Starlette never
    listens and keeps sending into the void."""
    data = bytes(range(256)) * (BIG // 256) + b"x" * (BIG % 256)
    path = _make_file(tmpdir_only, data)
    c = _run(TempFileResponse(path, media_type="text/csv"), _scope(range_header), disconnect_after_body_chunks=1)
    assert c.sent[0]["status"] == (206 if range_header else 200)
    body_msgs = [m for m in c.sent if m["type"] == "http.response.body"]
    assert len(body_msgs) == 1 and body_msgs[0]["more_body"] is True
    assert data.startswith(c.body()) and 0 < len(c.body()) < len(data)
    if STARLETTE_1:
        assert c.receives == ["http.request", "http.disconnect"]
        assert len(c.dropped) <= 1  # cancelled at its next await (0 in every local run)
    else:
        assert c.receives == [] and c.dropped  # kept sending after the disconnect
    # tmpdir_only asserts that the file is gone


# ---------------------------------------------------------- cleanup error policy


def test_missing_file_is_ignored_and_original_error_kept(tmpdir_only):
    path = _make_file(tmpdir_only)
    os.unlink(path)  # already gone before the response runs

    with pytest.raises(RuntimeError, match="does not exist"):  # starlette's own error, not FileNotFoundError
        _run(TempFileResponse(path), _scope())


def test_cleanup_error_after_success_is_raised(tmpdir_only, monkeypatch):
    path = _make_file(tmpdir_only)
    real_unlink = os.unlink

    def deny(p):
        raise PermissionError(13, "denied", p)

    monkeypatch.setattr(responses.os, "unlink", deny)

    try:
        with pytest.raises(PermissionError):
            _run(TempFileResponse(path), _scope())
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
                _run(TempFileResponse(path), _scope(), send)
        assert ei.value is err
        assert any("could not remove temporary file" in r.getMessage() for r in caplog.records)
        assert any(isinstance(r.exc_info[1], PermissionError) for r in caplog.records if r.exc_info)
    finally:
        real_unlink(path)
