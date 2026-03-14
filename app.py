"""ViNatX Experiment Monitor — real-time dashboard for Team 2 / Monomer Bio Hackathon"""

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import streamlit as st

DATA_DIR = Path(__file__).parent / "data"

# ─── Constants ────────────────────────────────────────────────────────────────

AUTOPLAT_URL     = "http://192.168.68.60:8080/mcp"
CLOUD_URL        = "https://backend-staging.monomerbio.com/mcp"
WORKFLOW_UUID    = "Jwx5nWNaOhTIHJVO3rhZ"
PLATE_BARCODE    = "ViNatX Tutorial Experiment Plate"
PLATE_ID         = "PLT1YUUDGALEIZHKHKSQPMF6PZLQ2Q"
REFRESH_INTERVAL = 30  # seconds

CONDITION_LABELS = {
    "A": "50% cells",
    "B": "30% cells",
    "C": "20% cells",
    "D": "15% cells",
    "E": "10% cells",
    "F": "5%  cells",
}

STATUS_ICON = {
    "completed":   "✅",
    "in_progress": "🔵",
    "initialized": "⚪",
    "failed":      "🔴",
    "canceled":    "🟡",
    "skipped":     "⚫",
}

STATUS_COLOR = {
    "completed":   "green",
    "in_progress": "blue",
    "initialized": "gray",
    "failed":      "red",
    "canceled":    "orange",
}


# ─── MCP HTTP client (sync, no asyncio) ───────────────────────────────────────

class MCPClient:
    """Minimal synchronous MCP StreamableHTTP client using requests."""

    def __init__(self, url: str, timeout: int = 20):
        self.url     = url
        self.timeout = timeout
        self._sid: str | None = None          # Mcp-Session-Id

    def _headers(self) -> dict:
        h = {"Content-Type": "application/json", "Accept": "application/json, text/event-stream"}
        if self._sid:
            h["Mcp-Session-Id"] = self._sid
        return h

    def _post(self, payload: dict) -> dict:
        resp = requests.post(self.url, json=payload, headers=self._headers(), timeout=self.timeout)
        resp.raise_for_status()
        if "Mcp-Session-Id" in resp.headers:
            self._sid = resp.headers["Mcp-Session-Id"]
        return self._parse(resp.text)

    @staticmethod
    def _parse(text: str) -> dict:
        text = text.strip()
        # SSE stream: pick first data: line
        if text.startswith("data:") or "\ndata:" in text:
            for line in text.splitlines():
                if line.startswith("data: "):
                    payload = line[6:].strip()
                    if payload and payload != "[DONE]":
                        try:
                            return json.loads(payload)
                        except json.JSONDecodeError:
                            pass
        # Plain JSON fallback
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return {}

    def call_tool(self, name: str, arguments: dict | None = None) -> dict | None:
        # initialize
        self._post({
            "jsonrpc": "2.0", "id": 0, "method": "initialize",
            "params": {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": {"name": "vinatx-dashboard", "version": "1.0"},
            },
        })
        # initialized notification (fire-and-forget, ignore errors)
        try:
            self._post({"jsonrpc": "2.0", "method": "notifications/initialized"})
        except Exception:
            pass
        # tool call
        result = self._post({
            "jsonrpc": "2.0", "id": 1, "method": "tools/call",
            "params": {"name": name, "arguments": arguments or {}},
        })
        content = result.get("result", {}).get("content", [])
        for item in content:
            if item.get("type") == "text":
                try:
                    return json.loads(item["text"])
                except json.JSONDecodeError:
                    return item["text"]
        return None

    def read_resource(self, uri: str) -> str | None:
        result = self._post({
            "jsonrpc": "2.0", "id": 2, "method": "resources/read",
            "params": {"uri": uri},
        })
        contents = result.get("result", {}).get("contents", [])
        for item in contents:
            if "text" in item:
                return item["text"]
        return None


def call_autoplat(name: str, arguments: dict | None = None):
    return MCPClient(AUTOPLAT_URL).call_tool(name, arguments)


def call_cloud(name: str, arguments: dict | None = None):
    return MCPClient(CLOUD_URL).call_tool(name, arguments)


# ─── Cached data fetchers ─────────────────────────────────────────────────────

@st.cache_data(ttl=REFRESH_INTERVAL)
def fetch_workflow():
    return call_autoplat(
        "get_workflow_instance_details",
        {"instance_uuid": WORKFLOW_UUID},
    )


@st.cache_data(ttl=10)
def load_obs_cache() -> tuple[pd.DataFrame | None, dict]:
    """Read OD data from local cache written by Claude Code."""
    csv_path  = DATA_DIR / "obs_cache.csv"
    meta_path = DATA_DIR / "obs_meta.json"

    meta = {}
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)

    if not csv_path.exists():
        return None, meta

    try:
        df = pd.read_csv(csv_path)
        df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True, errors="coerce")

        # Melt wide → long (one column per well → one row per well per timepoint)
        df_long = df.melt(id_vars=["timestamp"], var_name="well", value_name="OD600")
        df_long = df_long.dropna(subset=["OD600"])
        df_long["OD600"] = pd.to_numeric(df_long["OD600"], errors="coerce")
        df_long = df_long.dropna(subset=["OD600"])

        # Filter to our experiment wells (rows A-F, cols 2-4)
        exp_wells = {f"{r}{c}" for r in "ABCDEF" for c in (2, 3, 4)}
        df_long = df_long[df_long["well"].isin(exp_wells)].copy()
        df_long["row"]       = df_long["well"].str[0]
        df_long["condition"] = df_long["row"].map(CONDITION_LABELS)
        df_long = df_long.sort_values("timestamp")

        return df_long, meta
    except Exception as e:
        return None, {**meta, "parse_error": str(e)}


# ─── Page layout ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ViNatX Monitor",
    page_icon="🧫",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# Header
st.title("🧫 ViNatX Experiment Monitor")
st.caption(
    f"Team 2 · Monomer Bio Hackathon · "
    f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"
)

# Refresh controls
col_btn, col_tog, col_status = st.columns([1, 2, 4])
with col_btn:
    if st.button("🔄 Refresh now"):
        st.cache_data.clear()
        st.rerun()
with col_tog:
    auto_refresh = st.toggle("Auto-refresh (30s)", value=True)

# Auto-refresh countdown logic
if "last_refresh" not in st.session_state:
    st.session_state.last_refresh = time.time()

if auto_refresh:
    elapsed   = time.time() - st.session_state.last_refresh
    remaining = max(0, REFRESH_INTERVAL - int(elapsed))
    with col_status:
        st.caption(f"Next refresh in {remaining}s")
    if remaining == 0:
        st.session_state.last_refresh = time.time()
        st.cache_data.clear()
        st.rerun()
    else:
        time.sleep(1)
        st.rerun()

st.divider()

# ─── Section 1: Workflow progress ─────────────────────────────────────────────

st.subheader("Workflow Progress")

workflow = fetch_workflow()

if workflow is None:
    st.error("Could not connect to monomer-autoplat MCP server.")
else:
    routines  = workflow.get("workflow_routines", [])
    status    = workflow.get("status", "unknown")
    total     = len(routines)
    done      = sum(1 for r in routines if r["status"] == "completed")
    in_prog   = sum(1 for r in routines if r["status"] == "in_progress")

    # Top metrics
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Workflow status", status.upper())
    m2.metric("Steps completed", f"{done} / {total}")
    m3.metric("In progress",     str(in_prog))

    start_at = workflow.get("start_at")
    if start_at:
        start_dt = datetime.fromisoformat(start_at.replace("Z", "+00:00"))
        elapsed_min = int((datetime.now(timezone.utc) - start_dt).total_seconds() // 60)
        m4.metric("Elapsed", f"{elapsed_min} min")

        # Estimate finish: last OD read scheduled at ~80 min from start
        est_finish = start_dt.timestamp() + 90 * 60
        remaining_min = max(0, int((est_finish - time.time()) // 60))
        m5.metric("Est. remaining", f"~{remaining_min} min")

    # Step-by-step timeline
    st.markdown("**Steps:**")
    for i, r in enumerate(routines):
        icon  = STATUS_ICON.get(r["status"], "⚪")
        t_des = (r.get("desired_execution_time") or "")[:16].replace("T", " ")
        t_fin = (r.get("finished_execution_at") or "")[:16].replace("T", " ")
        label = r["routine_name"]

        if r["status"] == "completed":
            detail = f"finished {t_fin} UTC"
        elif r["status"] == "in_progress":
            detail = f"started ~{t_des} UTC — **running now**"
        else:
            detail = f"scheduled {t_des} UTC"

        st.markdown(f"{icon} **Step {i+1}** — {label} · {detail}")

    # Progress bar
    progress = done / total if total else 0
    st.progress(progress, text=f"{done}/{total} steps complete")

st.divider()

# ─── Section 2: OD600 growth curves ──────────────────────────────────────────

st.subheader("OD600 Growth Curves")

df_obs, obs_meta = load_obs_cache()

# Cache freshness info
cache_time = obs_meta.get("cache_written_at") or obs_meta.get("latest_observation_at")
n_datasets = obs_meta.get("total_datasets", 0)
if cache_time:
    cache_dt   = datetime.fromisoformat(cache_time.replace("Z", "+00:00"))
    age_min    = int((datetime.now(timezone.utc) - cache_dt).total_seconds() // 60)
    st.caption(
        f"📡 {n_datasets} timepoints · last fetched {age_min} min ago "
        f"({cache_dt.strftime('%H:%M UTC')}) · "
        f"*say \"refresh OD data\" to Claude Code to pull latest*"
    )

if "parse_error" in obs_meta:
    st.error(f"Cache parse error: {obs_meta['parse_error']}")
elif df_obs is None or df_obs.empty:
    routines_done = workflow.get("workflow_routines", []) if workflow else []
    transfer_done = any(
        r["routine_name"] == "Hackathon Transfer Samples" and r["status"] == "completed"
        for r in routines_done
    )
    if transfer_done:
        st.info("Transfers complete — waiting for first OD600 read to be cached.")
    else:
        st.info("Waiting for liquid handling to complete before first OD600 read.")
else:
    # Growth curves
    fig = px.line(
        df_obs,
        x="timestamp",
        y="OD600",
        color="condition",
        line_group="well",
        markers=True,
        title=f"OD600 over time — {PLATE_BARCODE}",
        labels={"timestamp": "Time (UTC)", "OD600": "Absorbance (OD600)"},
        height=420,
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(legend_title_text="Condition", hovermode="x unified")
    st.plotly_chart(fig, use_container_width=True)

    # Latest readings table
    latest = (
        df_obs.sort_values("timestamp")
        .groupby("well")
        .last()
        .reset_index()[["well", "condition", "OD600", "timestamp"]]
        .sort_values("well")
    )
    latest.columns = ["Well", "Condition", "OD600", "Last read (UTC)"]
    st.dataframe(latest, use_container_width=True, hide_index=True)

st.divider()

# ─── Section 3: Platemap ──────────────────────────────────────────────────────

st.subheader("Platemap — Tutorial Seeding Density Titration")

ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))

# Build z-values and annotations for the heatmap
z_vals = []
annots = []
hover  = []

for row in ROWS:
    z_row, a_row, h_row = [], [], []
    for col in COLS:
        if row in CONDITION_LABELS and col in (2, 3, 4):
            label   = CONDITION_LABELS[row]
            pct_str = label.split("%")[0].strip()
            pct     = float(pct_str)
            z_row.append(pct)
            a_row.append(label.replace(" cells", ""))
            h_row.append(f"<b>{row}{col}</b><br>{label}")
        else:
            z_row.append(None)
            a_row.append("")
            h_row.append(f"<b>{row}{col}</b><br>empty")
    z_vals.append(z_row)
    annots.append(a_row)
    hover.append(h_row)

fig2 = go.Figure(data=go.Heatmap(
    z=z_vals,
    x=[str(c) for c in COLS],
    y=ROWS,
    text=annots,
    customdata=hover,
    texttemplate="%{text}",
    hovertemplate="%{customdata}<extra></extra>",
    colorscale="Blues",
    showscale=True,
    colorbar=dict(title="Cell %"),
    zmin=0,
    zmax=55,
))

fig2.update_layout(
    title="96-well plate — experiment wells highlighted (cols 2–4, rows A–F)",
    xaxis_title="Column",
    yaxis_title="Row",
    yaxis_autorange="reversed",
    height=360,
)
st.plotly_chart(fig2, use_container_width=True)

st.caption(
    "**Conditions (triplicates):** "
    + " · ".join(f"{k}: {v}" for k, v in CONDITION_LABELS.items())
)
