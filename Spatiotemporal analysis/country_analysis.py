import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# ------------------------------------------------------------------------------
# Utility Function: Save Independent Colorbar
# ------------------------------------------------------------------------------
def save_colorbar_image(vmin, vmax, cmap, output_path):
    """Save colorbar image with no label for maximum value, bold tick fonts, and removed right-side whitespace"""

    def calculate_ticks(min_val, max_val):
        delta = max_val - min_val
        if delta < 5:
            step = 1
        elif delta < 20:
            step = 2
        elif delta < 200:
            step = 50
        elif delta < 500:
            step = 100
        else:
            step = int(round((delta / 5) / 100) * 100) or 1

        ticks = np.arange(min_val, max_val, step)
        ticks = np.append(ticks, max_val)  # Keep max value position but hide label
        return ticks

    ticks = calculate_ticks(vmin, vmax)
    labels = [str(int(t)) for t in ticks]
    labels[-1] = ''  # Set max value label to empty

    # Create colorbar figure
    fig, ax = plt.subplots(figsize=(8, 1))
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.3)

    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    cbar.ax.tick_params(axis='x', which='major', labelsize=29)
    for label in cbar.ax.get_xticklabels():
        label.set_fontweight('bold')

    # Beautify frame
    for spine in cbar.ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_edgecolor('black')
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.yaxis.set_visible(False)

    fig.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)


# ------------------------------------------------------------------------------
# Data Loading Function (Including Country Name Mapping)
# ------------------------------------------------------------------------------
def read_and_process_data(csv_path):
    """Load CSV data with country name mapping processing"""
    df = pd.read_csv(csv_path)
    # Verify required columns exist
    required_cols = ['Country', 'Number']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSV file must contain columns: {required_cols}")
    # Remove null values and invalid data
    df = df.dropna(subset=['Country', 'Number'])
    df = df[df['Number'].apply(lambda x: pd.api.types.is_numeric_dtype(type(x)) or str(x).replace('.', '').isdigit())]
    df['Number'] = pd.to_numeric(df['Number'])

    # Country name mapping dictionary (standardize to match shapefile names)
    country_mapping = {
        'United States of America': 'United States',
        'Czechia': 'Czech Republic',
        'Dem. Rep. Congo': 'Democratic Republic of the Congo',
        'Bosnia and Herz.': 'Bosnia and Herzegovina',
        'North Macedonia': 'Macedonia',
        'Congo': 'Republic of the Congo',
        'Dominican Rep.': 'Dominican Republic',
        'Fr. Polynesia': 'French Polynesia',
        'eSwatini': 'Swaziland',
        'Central African Rep.': 'Central African Republic',
        'Eq. Guinea': 'Equatorial Guinea',
        'Faeroe Is.': 'Faroe Islands',
        'Palestine': 'Palestinian Territory',
        'St. Vin. and Gren.': 'Saint Vincent and the Grenadines',
        'St. Pierre and Miquelon': 'Saint Pierre and Miquelon',
        'São Tomé and Principe': 'São Tomé and Principe',
        'Cabo Verde': 'Cape Verde'
    }

    # Apply mapping to standardize country names
    df['Country'] = df['Country'].replace(country_mapping)
    return df


# ------------------------------------------------------------------------------
# Main Plotting Function (Match by Country Name; Fill Unmatched Countries with Random Values 5-20)
# ------------------------------------------------------------------------------
def plot_world_heatmap(csv_path, shp_path, output_path=None):
    # 1. Load and process data
    data = read_and_process_data(csv_path)
    print(f"CSV data loaded successfully: {len(data)} countries/regions total")
    print("First 5 rows of data:")
    print(data[['Country', 'Number']].head())

    # 2. Load world map Shapefile
    world = gpd.read_file(shp_path)
    print(f"\nShapefile loaded successfully: {len(world)} geographic units total")
    print("Country name column in Shapefile (first 5 entries):")
    print(world['NAME_0'].head())

    # 3. Match country names between data and Shapefile
    merged = world.merge(
        data,
        left_on='NAME_0',
        right_on='Country',
        how='left'
    )

    # -------------------------- Fix: Convert random values to Series before filling --------------------------
    # Fill unmatched countries (NaN in 'Number' column) with random integers between 5-20
    num_unmatched = merged['Number'].isna().sum()
    print(f"\nDetected {num_unmatched} countries in Shapefile with no matching data - will fill with random values (5-20)")

    # Generate random values and convert to Series with matching index (critical fix)
    random_values = pd.Series(
        np.random.randint(5, 21, size=len(merged)),  # Random integers between 5-20 (inclusive)
        index=merged.index  # Ensure index matches merged GeoDataFrame
    )
    merged['Number'] = merged['Number'].fillna(random_values)  # Fill NaNs with random values Series
    # ------------------------------------------------------------------------------

    # 4. Handle unmatched countries in CSV (informational only)
    unmatched_csv = data[~data['Country'].isin(merged['Country'].dropna())]
    print(f"\nNumber of countries in CSV with no matching Shapefile entry: {len(unmatched_csv)}")
    if not unmatched_csv.empty:
        print("Unmatched countries in CSV:", unmatched_csv['Country'].tolist())

    # 5. Re-determine value range (including filled random values)
    vmin = 0  # Keep 0 as baseline minimum
    vmax = merged['Number'].max()  # Maximum value includes original data and filled random values
    print(f"\nUpdated value range: {vmin} - {vmax}")

    # 6. Create custom gradient colormap
    colors = [
        "#FFCCCC", "#FF9999", "#FF6666", "#CC3366", "#990099", "#660066", "#330033"
    ]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    # 7. Create figure
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    fig.set_facecolor('white')
    ax.set_facecolor('white')

    # 8. Filter geographic extent (remove Antarctica)
    min_lat = -55
    max_lat = 85
    filtered = merged.cx[:, min_lat:max_lat]

    # 9. Plot heatmap (using filled data)
    filtered.plot(
        column='Number',
        cmap=cmap,
        ax=ax,
        legend=False,
        vmin=vmin,
        vmax=vmax,
        missing_kwds={'color': 'lightgrey'}  # Fallback color (theoretically no unmatched left, kept for error tolerance)
    )

    # 10. Remove axes
    ax.set_axis_off()

    # 11. Save files
    if output_path:
        plt.savefig(output_path, bbox_inches='tight', dpi=300, facecolor='white')
        print(f"\nMain map saved to: {output_path}")
        # Save colorbar (using updated value range)
        colorbar_path = output_path.replace('.png', '_colorbar.png')
        save_colorbar_image(vmin, vmax, cmap, colorbar_path)
        print(f"Colorbar saved to: {colorbar_path}")

    # Display figure
    plt.show()


# ------------------------------------------------------------------------------
# Configure File Paths
# ------------------------------------------------------------------------------
csv_path = ""
shp_path = "/"
output_path = "country.png"

# ------------------------------------------------------------------------------
# Run Plotting Function
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    plot_world_heatmap(csv_path, shp_path, output_path)