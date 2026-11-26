"""Minimal ASGI app that returns what the server put in the scope (P1-2 V-b).

Run under the same gunicorn + worker class as the API to see the real
``asgi.spec_version`` / ``http_version`` (tests/test_p1_1b.py pins the value its
ASGI double uses).  Never part of the API.
"""
import json


async def app(scope, receive, send):
    if scope["type"] == "lifespan":
        while True:
            msg = await receive()
            if msg["type"] == "lifespan.startup":
                await send({"type": "lifespan.startup.complete"})
            elif msg["type"] == "lifespan.shutdown":
                await send({"type": "lifespan.shutdown.complete"})
                return
    first = await receive()
    body = json.dumps({"asgi": scope.get("asgi"), "http_version": scope.get("http_version"),
                       "type": scope["type"], "first_receive": first["type"]}).encode()
    await send({"type": "http.response.start", "status": 200,
                "headers": [(b"content-type", b"application/json"), (b"content-length", str(len(body)).encode())]})
    await send({"type": "http.response.body", "body": body})
