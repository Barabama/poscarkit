"""Network utilities — lightweight HTTP helpers using stdlib only."""

import json
import logging
import urllib.request
import urllib.error

log = logging.getLogger(__name__)

_GITHUB_TIMEOUT = 5  # seconds


def get_latest_version(repo: str) -> dict | None:
    """Query GitHub API for the latest release.

    Args:
        repo: 'Owner/Repo' string, e.g. 'Barabama/poscarkit'.

    Returns:
        dict with keys 'version', 'url', 'published_at' or None on failure.
    """
    url = f"https://api.github.com/repos/{repo}/releases/latest"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "poscarkit"})
    try:
        with urllib.request.urlopen(req, timeout=_GITHUB_TIMEOUT) as resp:
            data = json.loads(resp.read().decode())
            tag = data.get("tag_name", "")
            return {
                "version": tag.lstrip("v"),
                "url": data.get("html_url", url),
                "published_at": data.get("published_at", ""),
            }
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError, OSError) as exc:
        log.debug(f"Update check failed: {exc}")
        return None
