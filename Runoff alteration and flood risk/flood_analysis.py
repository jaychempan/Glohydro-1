#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flood Risk Analysis - World Map Visualization
Supports elliptical projection, displaying high/medium/low risk rivers and hydropower plants
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.validation import make_valid
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import json
import warnings

warnings.filterwarnings('ignore')

# ==== Configuration ====
FLOOD_GEOJSON_DIR = r""
RIVER_SHP_PATH = r"/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp"
HYDRO_CSV_PATH = r"v"
OUT_DIR = "./world_risk_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

# Risk level mapping (standardize various input formats)
RISK_MAPPING = {
    'high_risk': 'High Risk',
    'high risk': 'High Risk',
    'high': 'High Risk',
    'medium_risk': 'Medium Risk',
    'medium risk': 'Medium Risk',
    'medium': 'Medium Risk',
    'low_risk': 'Low Risk',
    'low risk': 'Low Risk',
    'low': 'Low Risk'
}

# Plot configuration
PLOT_CONFIG = {
    'risk_colors': {
        'High Risk': 'red',
        'Medium Risk': 'orange',
        'Low Risk': 'lightgreen'
    },
    'dpi': 300,
    'figsize': (16, 10),
    'projection': ccrs.Robinson()  # Use Robinson projection (ideal for world maps)
}


def make_geometry_valid(geom):
    """Fix invalid geometries"""
    try:
        if geom is None or geom.is_empty:
            return None
        if not geom.is_valid:
            # Attempt to fix geometry with make_valid
            fixed_geom = make_valid(geom)
            if fixed_geom.is_valid:
                return fixed_geom
            else:
                # If make_valid fails, try buffer(0)
                buffered = geom.buffer(0)
                if buffered.is_valid:
                    return buffered
                return None
        return geom
    except:
        return None


def extract_risk_level(properties):
    """Extract risk level from feature properties"""
    if properties is None or not hasattr(properties, 'get'):
        return 'Low Risk'

    # Check common risk level fields
    risk_fields = ['name', 'risk', 'risk_level', 'level', 'type', 'category']

    for field in risk_fields:
        try:
            if field in properties and properties[field] is not None:
                field_value = str(properties[field]).strip().lower()

                # Match risk level using mapping
                for key, risk_level in RISK_MAPPING.items():
                    if key in field_value:
                        return risk_level
        except:
            continue

    return 'Low Risk'


def load_flood_data():
    """Load flood risk data from GeoJSON files"""
    files = glob.glob(os.path.join(FLOOD_GEOJSON_DIR, "*.geojson"))

    all_features = []

    print("📊 Loading flood risk data...")

    for file_path in files:
        filename = os.path.basename(file_path)
        print(f"  Processing file: {filename}")

        try:
            # Read file with geopandas
            gdf = gpd.read_file(file_path)

            # Ensure correct CRS (WGS84)
            if gdf.crs is None or gdf.crs.to_epsg() != 4326:
                gdf = gdf.to_crs(4326)

            # Filter valid geometries
            gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
            gdf = gdf[gdf.geometry.geom_type.isin(['Polygon', 'MultiPolygon'])]

            if len(gdf) == 0:
                continue

            # Extract risk level for each feature
            risk_levels = []
            valid_geoms = []

            for idx, row in gdf.iterrows():
                # Fix geometry
                fixed_geom = make_geometry_valid(row.geometry)
                if fixed_geom is not None:
                    valid_geoms.append(fixed_geom)
                    # Extract risk level
                    risk_level = extract_risk_level(row)
                    risk_levels.append(risk_level)

            if valid_geoms:
                # Create valid GeoDataFrame
                valid_gdf = gpd.GeoDataFrame({
                    'geometry': valid_geoms,
                    'risk_level': risk_levels,
                    'source_file': filename
                }, crs=4326)

                all_features.append(valid_gdf)

                # Count risk distribution for this file
                risk_counts = pd.Series(risk_levels).value_counts()
                print(f"    Valid features: {len(valid_geoms)}")
                for risk, count in risk_counts.items():
                    print(f"      {risk}: {count} features")
            else:
                print(f"    ⚠️ No valid geometries")

        except Exception as e:
            print(f"⚠️ Skipping file {filename}: {e}")
            continue

    if all_features:
        # Merge all data
        combined_gdf = gpd.GeoDataFrame(
            pd.concat(all_features, ignore_index=True),
            crs=4326
        )

        print(f"✅ Successfully loaded {len(combined_gdf)} flood risk features")

        # Overall risk statistics
        risk_counts = combined_gdf['risk_level'].value_counts()
        print("📈 Flood risk feature statistics:")
        for risk, count in risk_counts.items():
            print(f"  {risk}: {count} features")

        return combined_gdf
    else:
        print("❌ No valid flood risk data found")
        return gpd.GeoDataFrame(columns=['geometry', 'risk_level', 'source_file'], crs=4326)


def assign_river_risk(rivers, flood_data):
    """Assign risk levels to rivers based on flood risk data"""
    print("\n🌊 Assigning risk levels to rivers...")

    # Initialize all rivers as Low Risk
    rivers = rivers.copy()
    rivers['risk_level'] = 'Low Risk'

    if len(flood_data) == 0:
        print("⚠️ No flood risk data available - all rivers marked as Low Risk")
        return rivers

    # Separate flood features by risk level
    high_risk_elements = flood_data[flood_data['risk_level'] == 'High Risk']
    medium_risk_elements = flood_data[flood_data['risk_level'] == 'Medium Risk']

    print(f"  High risk features: {len(high_risk_elements)}")
    print(f"  Medium risk features: {len(medium_risk_elements)}")

    # 1. Process high risk features first
    if len(high_risk_elements) > 0:
        print("  Processing high risk features...")
        try:
            # Spatial join to find intersecting rivers
            high_joined = gpd.sjoin(rivers, high_risk_elements, how='inner', predicate='intersects')
            rivers.loc[high_joined.index, 'risk_level'] = 'High Risk'
            print(f"    Marked {len(high_joined)} rivers as High Risk")
        except Exception as e:
            print(f"    Failed to process high risk features: {e}")

    # 2. Process medium risk features (only for rivers not yet marked as High Risk)
    if len(medium_risk_elements) > 0:
        print("  Processing medium risk features...")
        try:
            # Only consider rivers still marked as Low Risk
            low_risk_rivers = rivers[rivers['risk_level'] == 'Low Risk']
            if len(low_risk_rivers) > 0:
                medium_joined = gpd.sjoin(low_risk_rivers, medium_risk_elements, how='inner', predicate='intersects')
                rivers.loc[medium_joined.index, 'risk_level'] = 'Medium Risk'
                print(f"    Marked {len(medium_joined)} rivers as Medium Risk")
        except Exception as e:
            print(f"    Failed to process medium risk features: {e}")

    # Final river risk statistics
    risk_counts = rivers['risk_level'].value_counts()
    print(f"\n📊 River risk statistics:")
    for risk, count in risk_counts.items():
        print(f"  {risk}: {count} rivers")

    return rivers


def assign_hydro_risk(hydro_csv, rivers_risk):
    """Assign risk levels to hydropower plants based on nearby risk rivers"""
    print("\n⚡ Assigning risk levels to hydropower plants...")

    # Read hydropower plant data
    try:
        hydro_df = pd.read_csv(hydro_csv)
    except Exception as e:
        print(f"❌ Failed to read hydropower data: {e}")
        return gpd.GeoDataFrame(columns=['geometry', 'risk_level'], crs=4326)

    # Create valid geometry points
    valid_points = []
    valid_indices = []

    for i, row in hydro_df.iterrows():
        try:
            if (not pd.isna(row.Longitude) and not pd.isna(row.Latitude) and
                    isinstance(row.Longitude, (int, float)) and isinstance(row.Latitude, (int, float))):

                point = Point(float(row.Longitude), float(row.Latitude))
                if not point.is_empty:
                    valid_points.append(point)
                    valid_indices.append(i)
        except:
            continue

    if not valid_points:
        print("❌ No valid hydropower plant coordinates found")
        return gpd.GeoDataFrame(columns=['geometry', 'risk_level'], crs=4326)

    # Create GeoDataFrame for hydropower plants
    hydro_gdf = gpd.GeoDataFrame(
        hydro_df.iloc[valid_indices].copy(),
        geometry=valid_points,
        crs=4326
    )

    # Initialize all hydropower plants as Low Risk
    hydro_gdf['risk_level'] = 'Low Risk'

    # Only consider rivers with significant risk
    risk_rivers = rivers_risk[rivers_risk['risk_level'].isin(['High Risk', 'Medium Risk'])]

    if len(risk_rivers) == 0:
        print("⚠️ No risk rivers available - all hydropower plants marked as Low Risk")
        return hydro_gdf

    print("  Using grid-based fuzzy matching for hydropower plants...")

    # Grid-based fuzzy matching parameters
    grid_size = 0.1  # 0.1 degree grid (~11 km)

    # Create grid index for risk rivers
    river_grid = {}
    for idx, river in risk_rivers.iterrows():
        # Get river centroid
        center = river.geometry.centroid
        grid_x = int(center.x / grid_size)
        grid_y = int(center.y / grid_size)
        grid_key = (grid_x, grid_y)

        if grid_key not in river_grid:
            river_grid[grid_key] = []
        river_grid[grid_key].append((river.geometry, river.risk_level))

    # Assign risk levels to hydropower plants
    for i, hydro_row in hydro_gdf.iterrows():
        hydro_point = hydro_row.geometry

        # Determine grid cell for hydropower plant
        grid_x = int(hydro_point.x / grid_size)
        grid_y = int(hydro_point.y / grid_size)

        # Search for risk rivers in surrounding grid cells
        found_risk = 'Low Risk'
        min_distance = float('inf')

        # Search 3x3 grid area
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                search_key = (grid_x + dx, grid_y + dy)
                if search_key in river_grid:
                    for river_geom, risk_level in river_grid[search_key]:
                        try:
                            distance = hydro_point.distance(river_geom)
                            if distance < min_distance and distance <= 0.045:  # Within ~5 km
                                min_distance = distance
                                found_risk = risk_level
                        except:
                            continue

        if found_risk != 'Low Risk':
            hydro_gdf.at[i, 'risk_level'] = found_risk

    # Final hydropower plant risk statistics
    hydro_counts = hydro_gdf['risk_level'].value_counts()
    print(f"\n📊 Hydropower plant risk statistics:")
    for risk, count in hydro_counts.items():
        print(f"  {risk}: {count} hydropower plants")

    return hydro_gdf


def create_world_map_visualization(rivers_risk, hydro_risk, flood_data):
    """Create world map visualization of flood risks"""
    print("\n🌍 Creating world map visualization...")

    # 1. Create figure
    fig = plt.figure(figsize=PLOT_CONFIG['figsize'])

    # Create axes with projection
    ax = plt.axes(projection=PLOT_CONFIG['projection'])
    ax.set_global()  # Set global view

    # Add base map features
    ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
    ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, alpha=0.5)

    # Set title
    ax.set_title('Global Flood Risk Distribution Map', fontsize=16, fontweight='bold', pad=20)

    # 2. Transform data to projection CRS
    rivers_proj = rivers_risk.to_crs(PLOT_CONFIG['projection'].proj4_init)
    hydro_proj = hydro_risk.to_crs(PLOT_CONFIG['projection'].proj4_init)
    flood_proj = flood_data.to_crs(PLOT_CONFIG['projection'].proj4_init)

    # 3. Plot flood risk areas (sampled to avoid overcrowding)
    for risk_level, color in PLOT_CONFIG['risk_colors'].items():
        risk_elements = flood_proj[flood_proj['risk_level'] == risk_level]
        if len(risk_elements) > 0:
            # Sample features to improve rendering performance
            sample_size = min(100, len(risk_elements))
            sample = risk_elements.sample(sample_size, random_state=42)
            sample.plot(ax=ax, color=color, alpha=0.3,
                        label=f'Flood-{risk_level}', transform=PLOT_CONFIG['projection'])

    # 4. Plot rivers (sampled to avoid overcrowding)
    for risk_level, color in PLOT_CONFIG['risk_colors'].items():
        risk_rivers = rivers_proj[rivers_proj['risk_level'] == risk_level]
        if len(risk_rivers) > 0:
            # Sample rivers for better visualization
            sample_size = min(5000, len(risk_rivers))
            sample = risk_rivers.sample(sample_size, random_state=42)

            # Set line width based on risk level
            if risk_level == 'High Risk':
                linewidth = 1.2
            elif risk_level == 'Medium Risk':
                linewidth = 0.8
            else:  # Low Risk
                linewidth = 0.4

            sample.plot(ax=ax, color=color, linewidth=linewidth,
                        label=f'River-{risk_level}', transform=PLOT_CONFIG['projection'])

    # 5. Plot hydropower plants
    for risk_level, color in PLOT_CONFIG['risk_colors'].items():
        risk_hydro = hydro_proj[hydro_proj['risk_level'] == risk_level]
        if len(risk_hydro) > 0:
            # Set marker style based on risk level
            if risk_level == 'High Risk':
                markersize = 8
                marker = '^'  # Triangle
            elif risk_level == 'Medium Risk':
                markersize = 6
                marker = 's'  # Square
            else:  # Low Risk
                markersize = 4
                marker = 'o'  # Circle

            risk_hydro.plot(ax=ax, color=color, markersize=markersize, marker=marker,
                            edgecolor='white', linewidth=0.5,
                            label=f'Hydropower-{risk_level}', transform=PLOT_CONFIG['projection'])

    # 6. Add legend
    ax.legend(loc='lower left', bbox_to_anchor=(0, 0), frameon=True,
              fancybox=True, framealpha=0.8, fontsize=10)

    # 7. Add grid lines
    gl = ax.gridlines(draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "world_risk_distribution.png"),
                dpi=PLOT_CONFIG['dpi'], bbox_inches='tight')
    plt.close()

    # 8. Create risk statistics charts
    create_risk_statistics(rivers_risk, hydro_risk, flood_data)

    print("✅ World map visualization completed")


def create_risk_statistics(rivers_risk, hydro_risk, flood_data):
    """Create statistical charts of risk distributions"""
    print("📊 Creating risk statistics charts...")

    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Flood risk feature statistics
    if len(flood_data) > 0:
        flood_counts = flood_data['risk_level'].value_counts()
        bars1 = ax1.bar(flood_counts.index, flood_counts.values,
                        color=[PLOT_CONFIG['risk_colors'][r] for r in flood_counts.index])
        ax1.set_title('Flood Risk Area Statistics', fontweight='bold')
        ax1.set_ylabel('Number of Areas')
        ax1.tick_params(axis='x', rotation=45)

        # Add value labels
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height,
                     f'{int(height)}', ha='center', va='bottom', fontweight='bold')

    # 2. River risk statistics
    river_counts = rivers_risk['risk_level'].value_counts()
    bars2 = ax2.bar(river_counts.index, river_counts.values,
                    color=[PLOT_CONFIG['risk_colors'][r] for r in river_counts.index])
    ax2.set_title('River Risk Statistics', fontweight='bold')
    ax2.set_ylabel('Number of Rivers')
    ax2.tick_params(axis='x', rotation=45)

    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{int(height)}', ha='center', va='bottom', fontweight='bold')

    # 3. Hydropower plant risk statistics
    hydro_counts = hydro_risk['risk_level'].value_counts()
    bars3 = ax3.bar(hydro_counts.index, hydro_counts.values,
                    color=[PLOT_CONFIG['risk_colors'][r] for r in hydro_counts.index])
    ax3.set_title('Hydropower Plant Risk Statistics', fontweight='bold')
    ax3.set_ylabel('Number of Hydropower Plants')
    ax3.tick_params(axis='x', rotation=45)

    for bar in bars3:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{int(height)}', ha='center', va='bottom', fontweight='bold')

    # 4. Risk level distribution comparison
    risk_categories = ['High Risk', 'Medium Risk', 'Low Risk']
    river_values = [river_counts.get(r, 0) for r in risk_categories]
    hydro_values = [hydro_counts.get(r, 0) for r in risk_categories]

    # Calculate totals for validation
    river_total = sum(river_values)
    hydro_total = sum(hydro_values)

    if river_total > 0 and hydro_total > 0:
        # Create grouped bar chart
        x = np.arange(len(risk_categories))
        width = 0.35

        bars_river = ax4.bar(x - width / 2, river_values, width,
                             label='Rivers', alpha=0.7)
        bars_hydro = ax4.bar(x + width / 2, hydro_values, width,
                             label='Hydropower Plants', alpha=0.7)

        ax4.set_xlabel('Risk Level')
        ax4.set_ylabel('Count')
        ax4.set_title('Risk Level Distribution Comparison', fontweight='bold')
        ax4.set_xticks(x)
        ax4.set_xticklabels(risk_categories)
        ax4.legend()

        # Add value labels
        for bar in bars_river:
            height = bar.get_height()
            if height > 0:
                ax4.text(bar.get_x() + bar.get_width() / 2., height,
                         f'{int(height)}', ha='center', va='bottom')

        for bar in bars_hydro:
            height = bar.get_height()
            if height > 0:
                ax4.text(bar.get_x() + bar.get_width() / 2., height,
                         f'{int(height)}', ha='center', va='bottom')
    else:
        ax4.text(0.5, 0.5, 'No Data Available', ha='center', va='center',
                 transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Risk Level Distribution Comparison', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "risk_statistics.png"),
                dpi=PLOT_CONFIG['dpi'], bbox_inches='tight')
    plt.close()


def save_results(rivers_risk, hydro_risk, flood_data):
    """Save analysis results to files"""
    print("\n💾 Saving results...")

    try:
        # 1. Save flood risk data
        if len(flood_data) > 0:
            # Fix geometries before saving
            flood_fixed = flood_data.copy()
            valid_geoms = []
            valid_indices = []

            for idx, geom in enumerate(flood_data.geometry):
                fixed_geom = make_geometry_valid(geom)
                if fixed_geom is not None:
                    valid_geoms.append(fixed_geom)
                    valid_indices.append(idx)

            if valid_geoms:
                flood_fixed = flood_data.iloc[valid_indices].copy()
                flood_fixed['geometry'] = valid_geoms

                flood_fixed.to_file(os.path.join(OUT_DIR, "flood_risk.geojson"), driver='GeoJSON')
                print(f"✅ Saved {len(flood_fixed)} flood risk features")
    except Exception as e:
        print(f"❌ Failed to save flood risk data: {e}")

    try:
        # 2. Save river risk data (attributes only)
        if len(rivers_risk) > 0:
            rivers_attributes = rivers_risk.drop(columns=['geometry'])
            rivers_attributes.to_csv(os.path.join(OUT_DIR, "rivers_risk.csv"), index=False)
            print(f"✅ Saved risk data for {len(rivers_risk)} rivers")
    except Exception as e:
        print(f"❌ Failed to save river risk data: {e}")

    try:
        # 3. Save hydropower plant risk data
        if len(hydro_risk) > 0:
            # Fix geometries before saving
            hydro_fixed = hydro_risk.copy()
            valid_points = []
            valid_indices = []

            for idx, point in enumerate(hydro_risk.geometry):
                if point is not None and not point.is_empty and point.is_valid:
                    valid_points.append(point)
                    valid_indices.append(idx)

            if valid_points:
                hydro_fixed = hydro_risk.iloc[valid_indices].copy()
                hydro_fixed['geometry'] = valid_points

                hydro_fixed.to_file(os.path.join(OUT_DIR, "hydro_risk.geojson"), driver='GeoJSON')
                hydro_attributes = hydro_fixed.drop(columns=['geometry'])
                hydro_attributes.to_csv(os.path.join(OUT_DIR, "hydro_risk.csv"), index=False)
                print(f"✅ Saved risk data for {len(hydro_fixed)} hydropower plants")
    except Exception as e:
        print(f"❌ Failed to save hydropower plant data: {e}")

    # 4. Save summary statistics
    try:
        summary = {
            'flood_elements': {
                'High Risk': len(flood_data[flood_data['risk_level'] == 'High Risk']),
                'Medium Risk': len(flood_data[flood_data['risk_level'] == 'Medium Risk']),
                'Low Risk': len(flood_data[flood_data['risk_level'] == 'Low Risk']),
                'Total': len(flood_data)
            } if len(flood_data) > 0 else {},
            'rivers': {
                'High Risk': len(rivers_risk[rivers_risk['risk_level'] == 'High Risk']),
                'Medium Risk': len(rivers_risk[rivers_risk['risk_level'] == 'Medium Risk']),
                'Low Risk': len(rivers_risk[rivers_risk['risk_level'] == 'Low Risk']),
                'Total': len(rivers_risk)
            } if len(rivers_risk) > 0 else {},
            'hydro': {
                'High Risk': len(hydro_risk[hydro_risk['risk_level'] == 'High Risk']),
                'Medium Risk': len(hydro_risk[hydro_risk['risk_level'] == 'Medium Risk']),
                'Low Risk': len(hydro_risk[hydro_risk['risk_level'] == 'Low Risk']),
                'Total': len(hydro_risk)
            } if len(hydro_risk) > 0 else {}
        }

        with open(os.path.join(OUT_DIR, "summary.json"), 'w', encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        print("✅ Summary statistics saved")
    except Exception as e:
        print(f"❌ Failed to save summary statistics: {e}")


def main():
    """Main function"""
    print("=== Global Flood Risk Analysis ===")
    print("Strategy: Create world map visualization showing high/medium/low risk rivers and hydropower plants")

    # 1. Load flood risk data
    flood_data = load_flood_data()

    if len(flood_data) == 0:
        print("❌ No valid flood risk data available - program exiting")
        return

    # 2. Load river data
    print("\n🌊 Loading river data...")
    try:
        rivers = gpd.read_file(RIVER_SHP_PATH)
        if rivers.crs.to_epsg() != 4326:
            rivers = rivers.to_crs(4326)

        # Simplify geometries for better performance
        rivers['geometry'] = rivers.geometry.simplify(0.001)
        rivers = rivers[rivers.geometry.notna()]
        print(f"✅ Loaded {len(rivers)} rivers")
    except Exception as e:
        print(f"❌ Failed to load river data: {e}")
        return

    # 3. Assign risk levels to rivers
    rivers_risk = assign_river_risk(rivers, flood_data)

    # 4. Assign risk levels to hydropower plants
    hydro_risk = assign_hydro_risk(HYDRO_CSV_PATH, rivers_risk)

    # 5. Print detailed statistics
    print("\n" + "=" * 50)
    print("📊 Detailed Risk Analysis Results")
    print("=" * 50)

    if len(flood_data) > 0:
        flood_counts = flood_data['risk_level'].value_counts()
        print("\nFlood Risk Feature Statistics:")
        for risk, count in flood_counts.items():
            print(f"  {risk}: {count} features")

    if len(rivers_risk) > 0:
        river_counts = rivers_risk['risk_level'].value_counts()
        print("\nRiver Risk Statistics:")
        for risk, count in river_counts.items():
            print(f"  {risk}: {count} rivers")

    if len(hydro_risk) > 0:
        hydro_counts = hydro_risk['risk_level'].value_counts()
        print("\nHydropower Plant Risk Statistics:")
        for risk, count in hydro_counts.items():
            print(f"  {risk}: {count} hydropower plants")

    # 6. Create world map visualization
    create_world_map_visualization(rivers_risk, hydro_risk, flood_data)

    # 7. Save results
    save_results(rivers_risk, hydro_risk, flood_data)

    print("\n" + "=" * 50)
    print("🎯 Analysis Completed!")
    print(f"📁 Results saved to: {OUT_DIR}")
    print("=" * 50)


if __name__ == "__main__":
    main()