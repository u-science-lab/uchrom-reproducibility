"""Generate illustration assets with gpt-image-2.5, through OpenRouter (default) or the OpenAI Images API.

    python gen_image.py OUT_PREFIX "prompt" [--model openai/gpt-image-2.5-sunburst] [--provider openrouter|openai]
                        [--size 1536x1024] [--n 1] [--background transparent|opaque] [--image ref.png ...]

Reads OPENROUTER_API_KEY / OPENAI_API_KEY from the environment (source ~/.env in a subshell; the key is never
written anywhere).  OpenRouter: its images endpoint (OpenAI-compatible); OpenAI: images/generations (or images/edits with
--image).
With --image, the reference image(s) are edited / used as style references (images/edits).  The prompt and the
model are written next to the image (OUT_PREFIX.json) so that every asset can be regenerated.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import sys
import urllib.request
import urllib.error
import uuid


def post_json(url, body, key):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=900))


def post_multipart(url, fields, files, key):
    b = uuid.uuid4().hex
    out = []
    for k, v in fields.items():
        out += [f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode()]
    for k, path in files:
        data = open(path, "rb").read()
        out += [f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"; filename=\"{os.path.basename(path)}\"\r\n"
                f"Content-Type: image/png\r\n\r\n".encode(), data, b"\r\n"]
    out.append(f"--{b}--\r\n".encode())
    req = urllib.request.Request(url, data=b"".join(out),
                                 headers={"Authorization": f"Bearer {key}", "Content-Type": f"multipart/form-data; boundary={b}"})
    return json.load(urllib.request.urlopen(req, timeout=900))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("prompt")
    ap.add_argument("--model", default="openai/gpt-image-2.5-sunburst")
    ap.add_argument("--provider", default="openrouter", choices=["openrouter", "openai"])
    ap.add_argument("--size", default="1536x1024")
    ap.add_argument("--quality", default="high")
    ap.add_argument("--background", default="transparent")
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--image", nargs="*", default=[])
    a = ap.parse_args()
    if a.provider == "openrouter":
        openrouter(a)
        return
    key = os.environ["OPENAI_API_KEY"]
    params = {"model": a.model, "prompt": a.prompt, "size": a.size, "quality": a.quality, "n": a.n,
              "background": a.background}
    try:
        if a.image:
            r = post_multipart("https://api.openai.com/v1/images/edits", {k: str(v) for k, v in params.items()},
                               [("image[]", p) for p in a.image], key)
        else:
            r = post_json("https://api.openai.com/v1/images/generations", params, key)
    except urllib.error.HTTPError as e:
        print("HTTP", e.code, e.read().decode()[:800], file=sys.stderr)
        sys.exit(1)
    data = r.get("data") or []
    for i, d in enumerate(data):
        path = f"{a.out}_{i}.png" if len(data) > 1 else f"{a.out}.png"
        open(path, "wb").write(base64.b64decode(d["b64_json"]))
        print(path)
    json.dump({**{k: v for k, v in params.items()}, "image": a.image, "usage": r.get("usage")},
              open(f"{a.out}.json", "w"), indent=1)
    if not data:
        print("no image:", json.dumps(r)[:500], file=sys.stderr)
        sys.exit(1)


def openrouter(a) -> None:
    """OpenRouter's images endpoint (OpenAI-compatible request and response)."""
    key = os.environ["OPENROUTER_API_KEY"]
    params = {"model": a.model, "prompt": a.prompt, "size": a.size, "quality": a.quality, "n": a.n,
              "background": a.background}
    r = None
    for url in ("https://openrouter.ai/api/v1/images/generations", "https://openrouter.ai/api/v1/images"):
        try:
            r = post_json(url, params, key)
            break
        except urllib.error.HTTPError as e:
            err = f"HTTP {e.code} {url}: {e.read().decode()[:600]}"
            print(err, file=sys.stderr)
    if r is None:
        sys.exit(1)
    data = r.get("data") or []
    for i, d in enumerate(data):
        path = f"{a.out}_{i}.png" if len(data) > 1 else f"{a.out}.png"
        if d.get("b64_json"):
            open(path, "wb").write(base64.b64decode(d["b64_json"]))
        else:
            open(path, "wb").write(urllib.request.urlopen(d["url"], timeout=300).read())
        print(path)
    json.dump({"provider": "openrouter", **params, "usage": r.get("usage")}, open(f"{a.out}.json", "w"), indent=1)
    print("usage:", r.get("usage"))
    if not data:
        print("no image:", json.dumps(r)[:500], file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
