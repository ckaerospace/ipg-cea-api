"""CORS allowlist and optional API-key gate for the public Grok Build API."""

from __future__ import annotations

import os
import re

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

# https://*.grok.me, https://grok.me, https://grok.com, http://localhost:<any port>
CORS_ORIGIN_REGEX = (
    r"https://([a-zA-Z0-9-]+\.)*grok\.me"
    r"|https://([a-zA-Z0-9-]+\.)*grok\.com"
    r"|http://localhost(:\d+)?"
    r"|http://127\.0\.0\.1(:\d+)?"
    r"|https://([a-zA-Z0-9-]+\.)*onrender\.com"
)

_ORIGIN_RE = re.compile(rf"^(?:{CORS_ORIGIN_REGEX})$")

CORS_ALLOW_HEADERS = "Authorization, Content-Type, X-API-Key"
CORS_ALLOW_METHODS = "GET, POST, OPTIONS"

GRID_MIN = 17
GRID_MAX = 97


def origin_allowed(origin: str | None) -> bool:
    if not origin:
        return False
    return _ORIGIN_RE.fullmatch(origin) is not None


def cors_headers(origin: str | None) -> dict[str, str]:
    headers = {
        "Access-Control-Allow-Methods": CORS_ALLOW_METHODS,
        "Access-Control-Allow-Headers": CORS_ALLOW_HEADERS,
        "Access-Control-Max-Age": "86400",
        "Vary": "Origin",
    }
    if origin_allowed(origin):
        headers["Access-Control-Allow-Origin"] = origin  # type: ignore[assignment]
    return headers


def api_key() -> str:
    return os.environ.get("API_KEY", "").strip()


class ApiKeyMiddleware(BaseHTTPMiddleware):
    """If API_KEY is set, require matching X-API-Key on /api/* (except health + OPTIONS)."""

    async def dispatch(self, request: Request, call_next) -> Response:
        if request.method == "OPTIONS":
            return await call_next(request)
        key = api_key()
        if not key:
            return await call_next(request)
        path = request.url.path
        if path == "/api/health" or path == "/api/health/":
            return await call_next(request)
        if path.startswith("/api/"):
            provided = request.headers.get("x-api-key", "")
            if provided != key:
                body = JSONResponse({"detail": "Invalid or missing X-API-Key"}, status_code=401)
                origin = request.headers.get("origin")
                for k, v in cors_headers(origin).items():
                    body.headers[k] = v
                return body
        return await call_next(request)
