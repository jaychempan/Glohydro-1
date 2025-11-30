# -*- coding: utf-8 -*-
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ============ Paths ============
HYDRO_PLANTS_CSV = ""
RIVER_SHP = "HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp"
WORLD_SHP = "ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"

OUT_DIR = "./map_outputs"
os.makedirs(OUT_DIR, exist_ok=True)

# ============ Parameters (Adjustable) ============
BUFFER_M = 1000  # Hydropower plant-river matching distance (meters)
BORDER_BUFFER_DEG = 0.01  # Country border buffer (degrees) to avoid missing border crossings
PLANT_COLOR = "#5b3660"

# River filtering: Only plot "main stems" for readability (choose one or both)
USE_ORD_FLOW_FILTER = True
MAX_ORD_FLOW = 4  # Smaller = more "main stem", typical range 3~5
USE_FLOW_FILTER = False
MIN_DIS_AV_CMS = 200.0  # Only plot rivers with average discharge ≥ threshold (m3/s)

# Line simplification after dissolve (meters), recommended global range: 500~1500
SIMPLIFY_M = 800

# Visual styles
RIVER_UNDER = "#cfe2ee"  # Underlying river line color
RIVER_COLOR = "#15658a"  # Main river line color
ALPHA_UNDER = 0.55  # Transparency for underlying line
ALPHA_MAIN = 0.70  # Transparency for main line
W_UNDER = 2.4  # Line width for underlying river
W_MAIN = 1.2  # Line width for main river
PLANT_SIZE = 38  # Hydropower plant marker size

# Regional extents (lon_min, lon_max, lat_min, lat_max)
REGIONS = {
    "NA": [-170, -50, 5, 80],  # North America
    "SA": [-90, -30, -60, 15],  # South America
    "EU": [-15, 40, 32, 72],  # Europe
    "AF": [-20, 55, -35, 38],  # Africa
    "WA-SA": [35, 90, 5, 45],  # West Asia - South Asia
    "EA-SEA": [90, 140, 0, 55],  # East Asia - Southeast Asia
}


# ============ Utility Functions ============
def pick_river_id_col(rivers):
    """Select the most appropriate river ID column from common options"""
    for c in ["HYRIV_ID", "OBJECTID", "ID", "FID"]:
        if c in rivers.columns:
            return c
    rivers.reset_index(inplace=True)
    rivers.rename(columns={"index": "RIVER_ID"}, inplace=True)
    return "RIVER_ID"


def pick_country_col(world):
    """Select the country name column from common options"""
    if "ADMIN" in world.columns:
        return "ADMIN"
    if "NAME" in world.columns:
        return "NAME"
    raise ValueError("World boundary shapefile missing country name column (ADMIN/NAME).")


def dissolve_main_riv(rivers):
    """Dissolve multiple segments of the same MAIN_RIV into a single geometry (MultiLineString/LineString)"""
    keep_cols = [c for c in ["MAIN_RIV", "DIS_AV_CMS", "ORD_FLOW"] if c in rivers.columns]
    agg_dict = {}

    if "DIS_AV_CMS" in keep_cols:
        agg_dict["DIS_AV_CMS"] = "max"  # Preserve representative discharge (maximum value)
    if "ORD_FLOW" in keep_cols:
        agg_dict["ORD_FLOW"] = "min"  # Preserve main stem order (smaller = more main; take minimum for the MAIN_RIV)

    dissolved = rivers.dissolve(by="MAIN_RIV", aggfunc=agg_dict, as_index=False)
    return dissolved


def simplify_meter(gdf, tol_m):
    """Simplify geometry in metric CRS (3857) then convert back to WGS84 (4326)"""
    if tol_m <= 0:
        return gdf
    g3857 = gdf.to_crs(3857)  # Web Mercator (metric)
    g3857["geometry"] = g3857.geometry.simplify(tol_m, preserve_topology=True)
    return g3857.to_crs(4326)  # Convert back to WGS84


# ============ Main Workflow ============
def build_crossborder_fullrivers_with_plants():
    print("[1/5] Loading data ...")
    plants_df = pd.read_csv(HYDRO_PLANTS_CSV)
    plants = gpd.GeoDataFrame(
        plants_df,
        geometry=gpd.points_from_xy(plants_df["Longitude"], plants_df["Latitude"]),
        crs="EPSG:4326",
    )
    rivers = gpd.read_file(RIVER_SHP).to_crs(4326)
    world = gpd.read_file(WORLD_SHP).to_crs(4326)

    # Get key column names
    river_id_col = pick_river_id_col(rivers)
    country_col = pick_country_col(world)
    has_main_riv = "MAIN_RIV" in rivers.columns

    if not has_main_riv:
        raise ValueError("This script requires the MAIN_RIV field from HydroRIVERS dataset.")

    print("[2/5] Identifying transboundary rivers (with border buffer to avoid missing crossings) ...")
    # Buffer country boundaries to prevent missing river crossings at borders
    world_buf = world.copy()
    world_buf["geometry"] = world_buf.buffer(BORDER_BUFFER_DEG)

    # Spatial join between rivers and buffered countries
    river_country_join = gpd.sjoin(
        rivers[[river_id_col, "MAIN_RIV", "ORD_FLOW", "DIS_AV_CMS", "geometry"]],
        world_buf[[country_col, "geometry"]],
        how="left", predicate="intersects"
    ).dropna(subset=[country_col])

    # Count number of countries per river segment
    river_country_count = (river_country_join.drop_duplicates([river_id_col, country_col])
                           .groupby(river_id_col, as_index=False)
                           .agg(n_countries=(country_col, "nunique")))

    # Identify transboundary rivers (crossing ≥2 countries)
    transboundary_river_ids = set(river_country_count.loc[river_country_count["n_countries"] > 1, river_id_col])
    transboundary_main_rivs = set(rivers.loc[rivers[river_id_col].isin(transboundary_river_ids), "MAIN_RIV"])
    transboundary_full_rivers = rivers[rivers["MAIN_RIV"].isin(transboundary_main_rivs)].copy()

    print(
        f"  · Number of transboundary full rivers (by MAIN_RIV): {len(transboundary_main_rivs)}; Total segments: {len(transboundary_full_rivers)}")

    if transboundary_full_rivers.empty:
        return transboundary_full_rivers.iloc[0:0], plants.iloc[0:0]

    print(
        "[3/5] Matching hydropower plants to transboundary rivers (1000m buffer) → retaining full rivers and associated plants ...")
    # Convert to metric CRS for accurate distance calculations
    rivers_3857 = transboundary_full_rivers.to_crs(3857)[["MAIN_RIV", "geometry"]]
    plants_3857 = plants.to_crs(3857)[["geometry"]].copy()
    plants_3857["plant_id"] = np.arange(len(plants_3857))

    # Find nearest plants within buffer distance
    nearest_matches = gpd.sjoin_nearest(
        rivers_3857, plants_3857, how="left",
        max_distance=BUFFER_M, distance_col="distance_m"
    )
    matched_records = nearest_matches.dropna(subset=["plant_id"])

    if matched_records.empty:
        print("  · No hydropower plants found within 1000m of transboundary rivers.")
        return transboundary_full_rivers.iloc[0:0], plants.iloc[0:0]

    # Retain only rivers with matched plants and their associated plants
    retained_main_rivs = set(matched_records["MAIN_RIV"].unique())
    rivers_retained = transboundary_full_rivers[transboundary_full_rivers["MAIN_RIV"].isin(retained_main_rivs)].copy()
    plants_retained = plants.iloc[matched_records["plant_id"].astype(int).unique()].copy()

    print(
        f"  · Retained full rivers: {len(retained_main_rivs)}; Segments: {len(rivers_retained)}; Hydropower plants: {len(plants_retained)}")

    print("[4/5] River simplification (main stem filtering + dissolve + geometry simplification) ...")
    # Optional: Filter main stems based on flow order
    if USE_ORD_FLOW_FILTER and "ORD_FLOW" in rivers_retained.columns:
        rivers_retained = rivers_retained[rivers_retained["ORD_FLOW"] <= MAX_ORD_FLOW].copy()
        print(f"  · Segments after ORD_FLOW≤{MAX_ORD_FLOW} filtering: {len(rivers_retained)}")

    # Optional: Filter based on minimum discharge
    if USE_FLOW_FILTER and "DIS_AV_CMS" in rivers_retained.columns:
        rivers_retained = rivers_retained[rivers_retained["DIS_AV_CMS"] >= MIN_DIS_AV_CMS].copy()
        print(f"  · Segments after discharge≥{MIN_DIS_AV_CMS} filtering: {len(rivers_retained)}")

    # Dissolve segments into full rivers
    rivers_full = dissolve_main_riv(rivers_retained)
    # Simplify geometry (metric simplification)
    rivers_full = simplify_meter(rivers_full, SIMPLIFY_M)

    print(f"  · Number of full rivers after dissolve: {len(rivers_full)}")

    return rivers_full, plants_retained


# ============ Plotting ============
def draw_global_and_regions(rivers_gdf, plants_gdf):
    # ---- Global Overview Map ----
    fig = plt.figure(figsize=(22, 11))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())

    # Add base map features
    ax.add_feature(cfeature.LAND, facecolor="#f2f2f2")
    ax.add_feature(cfeature.OCEAN, facecolor="white")
    ax.coastlines(resolution="50m", linewidth=0.3, color="gray")
    ax.axis("off")

    # Plot rivers (double line effect for visibility)
    if not rivers_gdf.empty:
        river_geoms = list(rivers_gdf.geometry)
        ax.add_geometries(river_geoms, crs=ccrs.PlateCarree(),
                          edgecolor=RIVER_UNDER, facecolor="none",
                          linewidth=W_UNDER, alpha=ALPHA_UNDER, zorder=2, rasterized=True)
        ax.add_geometries(river_geoms, crs=ccrs.PlateCarree(),
                          edgecolor=RIVER_COLOR, facecolor="none",
                          linewidth=W_MAIN, alpha=ALPHA_MAIN, zorder=3, rasterized=True)

    # Plot hydropower plants
    if not plants_gdf.empty:
        ax.scatter(plants_gdf.geometry.x, plants_gdf.geometry.y,
                   s=PLANT_SIZE, color=PLANT_COLOR, alpha=0.85,
                   transform=ccrs.PlateCarree(), zorder=4,
                   edgecolors="white", linewidth=0.35)

    # Save global map
    global_out_path = os.path.join(OUT_DIR, "crossborder_fullrivers_overview.png")
    plt.savefig(global_out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {global_out_path}")

    # ---- Regional Detailed Maps ----
    for region_name, (lon_min, lon_max, lat_min, lat_max) in REGIONS.items():
        fig = plt.figure(figsize=(12, 9))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

        # Add base map features
        ax.add_feature(cfeature.LAND, facecolor="#f3f3f3")
        ax.add_feature(cfeature.OCEAN, facecolor="white")
        ax.coastlines(resolution="50m", linewidth=0.4, color="gray")
        ax.axis("off")

        # Subset data to current region
        regional_rivers = rivers_gdf.cx[lon_min:lon_max, lat_min:lat_max]
        regional_plants = plants_gdf.cx[lon_min:lon_max, lat_min:lat_max]

        # Plot regional rivers
        if not regional_rivers.empty:
            ax.add_geometries(list(regional_rivers.geometry), crs=ccrs.PlateCarree(),
                              edgecolor=RIVER_UNDER, facecolor="none",
                              linewidth=W_UNDER, alpha=ALPHA_UNDER, zorder=2, rasterized=True)
            ax.add_geometries(list(regional_rivers.geometry), crs=ccrs.PlateCarree(),
                              edgecolor=RIVER_COLOR, facecolor="none",
                              linewidth=W_MAIN, alpha=ALPHA_MAIN, zorder=3, rasterized=True)

        # Plot regional plants
        if not regional_plants.empty:
            ax.scatter(regional_plants.geometry.x, regional_plants.geometry.y,
                       s=PLANT_SIZE, color=PLANT_COLOR, alpha=0.9,
                       transform=ccrs.PlateCarree(), zorder=4,
                       edgecolors="white", linewidth=0.4)

        # Save regional map
        regional_out_path = os.path.join(OUT_DIR, f"region_{region_name}.png")
        plt.savefig(regional_out_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"Saved: {regional_out_path}")


# ============ Entry Point ============
def main():
    rivers_full, plants_keep = build_crossborder_fullrivers_with_plants()
    if rivers_full.empty:
        print("❗ No data available for plotting.")
        return
    draw_global_and_regions(rivers_full, plants_keep)
    print("✅ Completed.")


if __name__ == "__main__":
    main()