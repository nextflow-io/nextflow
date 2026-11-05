# Copyright 2013-2026, Seqera Labs
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Scripted stand-in for the Anthropic Messages API, used by smoke-test.sh.

It drives the harness through the runner-native tools: `read` of a small PNG, then `bash`. Then it
checks what came back and answers with `SMOKE-OK`, or with `SMOKE-FAIL: <reason>`, as the agent's
final text.
"""

import json
import struct
import sys
import zlib
from http.server import BaseHTTPRequestHandler, HTTPServer

IMAGE = "smoke.png"
BASH_MARKER = "smoke-bash-ok"


def write_png(path, size=64):
    raw = b"".join(b"\x00" + bytes(v for x in range(size) for v in (x * 4, y * 4, 128)) for y in range(size))

    def chunk(tag, data):
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    header = struct.pack(">IIBBBBB", size, size, 8, 2, 0, 0, 0)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", header) + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


def message(blocks, stop_reason):
    events = [{"type": "message_start", "message": {
        "id": "msg", "type": "message", "role": "assistant", "model": "mock", "content": [],
        "stop_reason": None, "usage": {"input_tokens": 1, "output_tokens": 1}}}]
    for i, (start, delta) in enumerate(blocks):
        events += [
            {"type": "content_block_start", "index": i, "content_block": start},
            {"type": "content_block_delta", "index": i, "delta": delta},
            {"type": "content_block_stop", "index": i}]
    events += [
        {"type": "message_delta", "delta": {"stop_reason": stop_reason}, "usage": {"output_tokens": 1}},
        {"type": "message_stop"}]
    return "".join(f"event: {e['type']}\ndata: {json.dumps(e)}\n\n" for e in events).encode()


def tool_use(call_id, name, args):
    return message([(
        {"type": "tool_use", "id": call_id, "name": name, "input": {}},
        {"type": "input_json_delta", "partial_json": json.dumps(args)})], "tool_use")


def text(value):
    return message([({"type": "text", "text": ""}, {"type": "text_delta", "text": value})], "end_turn")


def tool_results(body):
    return {
        block["tool_use_id"]: block
        for msg in body["messages"] if msg["role"] == "user" and isinstance(msg["content"], list)
        for block in msg["content"] if block.get("type") == "tool_result"}


def parts(result):
    content = result.get("content")
    return content if isinstance(content, list) else [{"type": "text", "text": content or ""}]


def verdict(results):
    read = results["read-1"]
    if read.get("is_error") or not any(p["type"] == "image" for p in parts(read)):
        return f"SMOKE-FAIL: read returned no image: {json.dumps(parts(read))[:300]}"
    bash = results["bash-1"]
    if bash.get("is_error") or BASH_MARKER not in json.dumps(parts(bash)):
        return f"SMOKE-FAIL: bash did not run: {json.dumps(parts(bash))[:300]}"
    return "SMOKE-OK"


class Handler(BaseHTTPRequestHandler):

    def do_POST(self):
        body = json.loads(self.rfile.read(int(self.headers["content-length"])))
        results = tool_results(body)
        if "read-1" not in results:
            out = tool_use("read-1", "read", {"path": IMAGE})
        elif "bash-1" not in results:
            out = tool_use("bash-1", "bash", {"command": f"echo {BASH_MARKER}"})
        else:
            out = text(verdict(results))
        self.send_response(200)
        self.send_header("content-type", "text/event-stream")
        self.send_header("content-length", str(len(out)))
        self.end_headers()
        self.wfile.write(out)

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    write_png(f"{sys.argv[1]}/{IMAGE}")
    HTTPServer(("0.0.0.0", 8080), Handler).serve_forever()
