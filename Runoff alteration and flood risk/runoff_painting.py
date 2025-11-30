import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --------------------------
# 1. Basic Configuration & Data Loading (Local Path Adaptation)
# --------------------------
# Set font (compatible with multiple systems)
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# Local data file path (verify path correctness)
file_path = 'v'
# Output image path (same folder)
output_path = ''

# Load data and check basic information
df = pd.read_csv(file_path)
print(f"Data loaded successfully! Contains {len(df)} stations in total")
print(f"Original trend type distribution:\n{df['Trend'].value_counts()}")


# --------------------------
# 2. Data Preprocessing: Define Three Trend Categories
# --------------------------
# Standardize trend category names (simplify to three types)
def simplify_trend_category(trend):
    if trend == 'Significant Decreasing':
        return 'Significant Decreasing Trend'
    elif trend == 'Significant Increasing':
        return 'Significant Increasing Trend'
    else:  # No Significant Trend
        return 'No Significant Trend'


df['Trend_Category'] = df['Trend'].apply(simplify_trend_category)

# Count stations for each of the three trend categories
trend_count = df['Trend_Category'].value_counts()
print(f"\nSimplified three-category trend distribution:")
for cat, count in trend_count.items():
    print(f"  {cat}: {count} stations ({count / len(df) * 100:.1f}%)")

# --------------------------
# 3. Vibrant Morandi Color Scheme (Exclusive Colors for Three Categories)
# --------------------------
trend_colors = {
    'Significant Decreasing Trend': '#cf3d3e',  # Vibrant Morandi orange-red (visually striking)
    'Significant Increasing Trend': '#80a6e2',  # Vibrant Morandi teal-blue (strong contrast)
    'No Significant Trend': '#F8E9D6'  # Very light Morandi (placed at bottom, non-distracting)
}

# --------------------------
# 4. Plot Robinson Projection Map (No Country Borders)
# --------------------------
# Create figure (Robinson projection, balanced proportions)
fig, ax = plt.subplots(1, 1, figsize=(18, 10), subplot_kw={'projection': ccrs.Robinson()})

# --------------------------
# 5. Map Style Settings (Core: Remove Country Borders)
# --------------------------
# Map background color (light gray, clean and simple)
ax.set_facecolor('#f9f9f9')
# Land color (light gray, only basic land-sea distinction)
ax.add_feature(
    cfeature.LAND,
    facecolor='#e5e5e5',
    zorder=1,  # Land at the bottom layer
    edgecolor='#d0d0d0',  # Land edge line (light gray, weakened border sense)
    linewidth=0.3
)
# Ocean color (white, soft transition with background)
ax.add_feature(
    cfeature.OCEAN,
    facecolor='white',
    zorder=1  # Ocean at the same layer as land
)
# 【Core Modification】Remove country borders: Delete original cfeature.BORDERS code
# Keep only coastlines (maintain land-sea outline without national border interference)
ax.add_feature(
    cfeature.COASTLINE,
    edgecolor='#b0b0b0',
    linewidth=0.6,
    zorder=2
)

# Set global extent and remove Antarctica (latitude range: -60 to 90)
ax.set_global()
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# --------------------------
# 6. Plot Three Trend Categories (Bottom-placement + No Borders)
# --------------------------
point_size = 16
# Drawing order: No Significant Trend (bottom) → Significant Increasing → Significant Decreasing (top)
draw_order = [
    ('No Significant Trend', 3),
    ('Significant Increasing Trend', 4),
    ('Significant Decreasing Trend', 5)
]

for trend_cat, z_order in draw_order:
    cat_data = df[df['Trend_Category'] == trend_cat]
    if len(cat_data) == 0:
        continue

    # Plot points without borders + layered drawing
    ax.scatter(
        x=cat_data['Longitude'],
        y=cat_data['Latitude'],
        color=trend_colors[trend_cat],
        s=point_size,
        alpha=0.7,  # Transparency balances visibility and overlap effect
        edgecolor=None,  # No borders
        linewidth=0,     # No borders
        transform=ccrs.PlateCarree(),
        label=f'{trend_cat} ({len(cat_data)} stations)',
        zorder=z_order
    )

# --------------------------
# 7. Chart Enhancement & Annotations
# --------------------------
# Legend settings (clearly display categories and counts)
ax.legend(
    loc='lower left',
    bbox_to_anchor=(0.02, 0.02),
    frameon=True,
    fancybox=True,
    framealpha=0.9,
    fontsize=12,
    labelspacing=0.8,
    handletextpad=1.2,
    borderpad=0.8
)

# Remove map frame for cleaner look
ax.spines['geo'].set_visible(False)

# Add title (clear and concise)
ax.set_title(
    'Global Distribution of Hydropower Station Runoff Trends',
    fontsize=18,
    fontweight='bold',
    pad=20
)

# --------------------------
# 8. Save & Display
# --------------------------
plt.tight_layout()
plt.savefig(
    output_path,
    dpi=300,  # High definition
    bbox_inches='tight',  # Avoid element truncation
    facecolor='white',
    edgecolor='none'
)

plt.show()

print(f"\nChart successfully saved to: {output_path}")
print("Key features of current version:")
print("1. Country borders removed - only coastlines and land-sea distinction retained for a cleaner style")
print("2. 'No Significant Trend' points placed at the bottom layer (zorder=3) to avoid obscuring other categories")
print("3. All points have no borders (edgecolor=None + linewidth=0)")
print("4. Three Morandi colors with clear distinction and distinct visual hierarchy")