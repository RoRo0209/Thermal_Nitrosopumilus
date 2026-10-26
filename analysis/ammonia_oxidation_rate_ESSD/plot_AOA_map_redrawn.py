#!/usr/bin/env python3
"""Redraw the AOA sampling map from source data without using the old PDF."""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'ammonia_oxidation_rate_ESSD')


import pandas as pd
import plotly.graph_objects as go


DATA = "Database of nitrification and nitrifiers in the global ocean.xlsx"
RATE_COLOR = "#0072B2"
AOA_COLOR = "#D55E00"


def find_col(columns, token):
    for col in columns:
        if token.lower() in col.lower():
            return col
    raise KeyError(f"Could not find column containing {token!r}")


def load_filtered_points():
    df_aoa = pd.read_excel(DATA, sheet_name="ammonia oxidizer", header=0)
    df_rate = pd.read_excel(DATA, sheet_name="volumetric ammonia oxidation", header=0)
    df_aoa.columns = [c.strip() for c in df_aoa.columns]
    df_rate.columns = [c.strip() for c in df_rate.columns]

    temp_aoa = find_col(df_aoa.columns, "Temperature")
    temp_rate = find_col(df_rate.columns, "Temperature")

    df_aoa = df_aoa.rename(
        columns={
            "Latitude": "lat",
            "Longitude": "lon",
            "Depth (m)": "depth",
            "total AOA (copies L-1)": "total_aoa",
            "shallow clade AOA (copies L-1)": "shallow_aoa",
            temp_aoa: "temp",
        }
    )
    df_rate = df_rate.rename(
        columns={
            "Latitude": "lat",
            "Longitude": "lon",
            "Depth (m)": "depth",
            "Rate (nmol N L-1 d-1)": "rate",
            temp_rate: "temp",
        }
    )

    df_aoa["aoa"] = df_aoa["total_aoa"].fillna(df_aoa["shallow_aoa"])

    for col in ["lat", "lon", "depth", "temp", "aoa"]:
        df_aoa[col] = pd.to_numeric(df_aoa[col], errors="coerce")
    for col in ["lat", "lon", "depth", "temp", "rate"]:
        df_rate[col] = pd.to_numeric(df_rate[col], errors="coerce")

    df_aoa_s = df_aoa[
        (df_aoa["depth"] <= 200)
        & df_aoa["temp"].notna()
        & df_aoa["aoa"].notna()
        & (df_aoa["aoa"] > 0)
        & df_aoa["lon"].notna()
        & df_aoa["lat"].notna()
    ].copy()
    df_rate_s = df_rate[
        (df_rate["depth"] <= 200)
        & df_rate["temp"].notna()
        & df_rate["rate"].notna()
        & (df_rate["rate"] >= 0)
        & df_rate["lon"].notna()
        & df_rate["lat"].notna()
    ].copy()
    return df_rate_s, df_aoa_s


def main():
    df_rate_s, df_aoa_s = load_filtered_points()

    fig = go.Figure()
    fig.add_trace(
        go.Scattergeo(
            lon=df_rate_s["lon"],
            lat=df_rate_s["lat"],
            mode="markers",
            name=f"AO rate (n={len(df_rate_s)})",
            marker=dict(
                size=7.5,
                color=RATE_COLOR,
                opacity=0.82,
                symbol="circle",
                line=dict(width=0.7, color="white"),
            ),
        )
    )
    fig.add_trace(
        go.Scattergeo(
            lon=df_aoa_s["lon"],
            lat=df_aoa_s["lat"],
            mode="markers",
            name=f"AOA abundance (n={len(df_aoa_s)})",
            marker=dict(
                size=8.0,
                color=AOA_COLOR,
                opacity=0.82,
                symbol="square",
                line=dict(width=0.7, color="white"),
            ),
        )
    )
    fig.update_geos(
        domain=dict(x=[0.0, 1.0], y=[0.0, 0.94]),
        projection_type="equirectangular",
        showland=True,
        landcolor="#F1F1F1",
        showocean=True,
        oceancolor="white",
        showcountries=False,
        showcoastlines=True,
        coastlinecolor="#303030",
        coastlinewidth=0.8,
        lonaxis=dict(showgrid=True, gridcolor="#E6E6E6", gridwidth=0.6, dtick=60),
        lataxis=dict(showgrid=True, gridcolor="#E6E6E6", gridwidth=0.6, dtick=30),
    )
    fig.update_layout(
        width=1800,
        height=820,
        title=dict(text="Sampling locations (depth <= 200 m)", x=0.5, y=0.985, font=dict(size=24)),
        margin=dict(l=18, r=18, t=42, b=14),
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=dict(family="Arial, Helvetica, sans-serif", size=18, color="black"),
        legend=dict(
            x=0.015,
            y=0.035,
            xanchor="left",
            yanchor="bottom",
            bgcolor="rgba(255,255,255,0.92)",
            bordercolor="#D0D0D0",
            borderwidth=1,
            font=dict(size=15),
            itemsizing="constant",
        ),
    )

    fig.write_image("Fig_AOA_map_redrawn.pdf")
    fig.write_image("Fig_AOA_map_redrawn.png", scale=2)
    print("Saved: Fig_AOA_map_redrawn.pdf/png")


if __name__ == "__main__":
    main()
