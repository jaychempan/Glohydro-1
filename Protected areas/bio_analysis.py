import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Reset default matplotlib settings to ensure proper English display
plt.rcParams.update(plt.rcParamsDefault)

# 1. Data Loading and Preprocessing (Please modify the path according to your actual file location; this is an example path)
# If the file is in the Downloads folder, the original path can be used directly: '/Users/lijiahao/Downloads/biomass_data_with_plant_info.csv'
df = pd.read_csv('')
# Calculate relative year (construction year = 0)
df['relative_year'] = df['year'] - df['plant_year']

# 2. Group by Capacity (Using the median 40.5 MW as the classification standard)
capacity_median = 40.5
df['capacity_category'] = df['capacity_mw'].apply(
    lambda x: 'Large Hydropower Plants' if x > capacity_median else 'Small Hydropower Plants'
)

# 3. Data Filtering and Statistics (Starting from -3 years, sample size ≥3 to ensure reliability)
capacity_stats = df.groupby(['capacity_category', 'relative_year']).agg({
    'filled_biomass_Mg': ['mean', 'std', 'count']
}).reset_index()
capacity_stats.columns = ['capacity_category', 'relative_year', 'mean_biomass', 'std_biomass', 'sample_count']
capacity_stats_filtered = capacity_stats[(capacity_stats['relative_year'] >= -3) &
                                         (capacity_stats['sample_count'] >= 3)].copy()

# Overall average line statistics
overall_stats = df[df['relative_year'] >= -3].groupby('relative_year')['filled_biomass_Mg'].mean().reset_index()
overall_stats.columns = ['relative_year', 'overall_mean']

# 4. Plot Biomass Trend Chart (Fixed: Correctly set x and y tick fonts to 16pt bold)
fig, ax = plt.subplots(1, 1, figsize=(14, 8))

# Style Configuration (Maintain professional consistency)
styles = {
    'Large Hydropower Plants': {'color': '#2E86AB', 'marker': 'o', 'linestyle': '-', 'linewidth': 3},
    'Small Hydropower Plants': {'color': '#A23B72', 'marker': 's', 'linestyle': '--', 'linewidth': 3},
    'Overall': {'color': '#333333', 'marker': 'D', 'linestyle': '-.', 'linewidth': 2.5}
}

# Plot Large Hydropower Plants trend
large_data = capacity_stats_filtered[capacity_stats_filtered['capacity_category'] == 'Large Hydropower Plants']
ax.plot(large_data['relative_year'], large_data['mean_biomass'],
        color=styles['Large Hydropower Plants']['color'],
        marker=styles['Large Hydropower Plants']['marker'],
        linestyle=styles['Large Hydropower Plants']['linestyle'],
        linewidth=styles['Large Hydropower Plants']['linewidth'],
        markersize=8, markerfacecolor='white', markeredgewidth=2,
        markeredgecolor=styles['Large Hydropower Plants']['color'])
ax.fill_between(large_data['relative_year'],
                large_data['mean_biomass'] - large_data['std_biomass'],
                large_data['mean_biomass'] + large_data['std_biomass'],
                alpha=0.15, color=styles['Large Hydropower Plants']['color'])

# Plot Small Hydropower Plants trend
small_data = capacity_stats_filtered[capacity_stats_filtered['capacity_category'] == 'Small Hydropower Plants']
ax.plot(small_data['relative_year'], small_data['mean_biomass'],
        color=styles['Small Hydropower Plants']['color'],
        marker=styles['Small Hydropower Plants']['marker'],
        linestyle=styles['Small Hydropower Plants']['linestyle'],
        linewidth=styles['Small Hydropower Plants']['linewidth'],
        markersize=8, markerfacecolor='white', markeredgewidth=2,
        markeredgecolor=styles['Small Hydropower Plants']['color'])
ax.fill_between(small_data['relative_year'],
                small_data['mean_biomass'] - small_data['std_biomass'],
                small_data['mean_biomass'] + small_data['std_biomass'],
                alpha=0.15, color=styles['Small Hydropower Plants']['color'])

# Plot Overall Average line
ax.plot(overall_stats['relative_year'], overall_stats['overall_mean'],
        color=styles['Overall']['color'],
        marker=styles['Overall']['marker'],
        linestyle=styles['Overall']['linestyle'],
        linewidth=styles['Overall']['linewidth'],
        markersize=6, markerfacecolor='white', markeredgewidth=1.5,
        markeredgecolor=styles['Overall']['color'])

# Plot Construction Year marker line
ax.axvline(x=0, color='#C73E1D', linestyle=':', linewidth=2.5, zorder=1)

# Axis and Title Settings
ax.set_xlabel('Relative Year (Newly-built Year = 0)', fontsize=14, fontweight='bold')
ax.set_ylabel('Average Biomass (Mg)', fontsize=14, fontweight='bold')
#ax.set_title('Biomass Variation Trends Before and After Hydropower Construction\nby Hydropower Plant Capacity',
#            fontsize=16, fontweight='bold', pad=20)

# Core Fix: Correctly set x and y axis tick fonts (16pt bold)
# X-axis: Set tick positions + tick label fonts
x_ticks = range(-3, 6)
ax.set_xticks(x_ticks)
ax.set_xticklabels([str(tick) for tick in x_ticks],
                   fontsize=16, fontweight='bold')  # X-ticks: 16pt bold
ax.set_xlim([-3.5, 5.5])

# Y-axis: Get default ticks first, then set fonts (avoid bias from manual tick values)
y_ticks = ax.get_yticks()
ax.set_yticklabels([f'{tick:.0f}' for tick in y_ticks],
                   fontsize=16, fontweight='bold')  # Y-ticks: 16pt bold

# Grid Settings
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8, zorder=0)
ax.set_axisbelow(True)

# Save Trend Chart (Saved to Documents folder; modify path as needed)
save_path_trend = 'final_biomass_trend_chart.png'
plt.tight_layout()
plt.savefig(save_path_trend, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

# 5. Plot Independent Legend Chart (Synchronously fix tick fonts; hide axes but maintain configuration consistency)
fig_legend, ax_legend = plt.subplots(1, 1, figsize=(10, 6))
ax_legend.axis('off')  # Hide axes without affecting display

# Create Legend Elements (Exact match with trend chart styles)
legend_elements = [
    plt.Line2D([0], [0],
               color=styles['Large Hydropower Plants']['color'],
               marker=styles['Large Hydropower Plants']['marker'],
               linestyle=styles['Large Hydropower Plants']['linestyle'],
               linewidth=styles['Large Hydropower Plants']['linewidth'],
               markersize=8, markerfacecolor='white', markeredgewidth=2,
               markeredgecolor=styles['Large Hydropower Plants']['color'],
               label=f'Large Hydropower Plants\n(≥ {capacity_median:.1f} MW)'),
    plt.Line2D([0], [0],
               color=styles['Small Hydropower Plants']['color'],
               marker=styles['Small Hydropower Plants']['marker'],
               linestyle=styles['Small Hydropower Plants']['linestyle'],
               linewidth=styles['Small Hydropower Plants']['linewidth'],
               markersize=8, markerfacecolor='white', markeredgewidth=2,
               markeredgecolor=styles['Small Hydropower Plants']['color'],
               label=f'Small Hydropower Plants\n(< {capacity_median:.1f} MW)'),
    plt.Line2D([0], [0],
               color=styles['Overall']['color'],
               marker=styles['Overall']['marker'],
               linestyle=styles['Overall']['linestyle'],
               linewidth=styles['Overall']['linewidth'],
               markersize=6, markerfacecolor='white', markeredgewidth=1.5,
               markeredgecolor=styles['Overall']['color'],
               label='Overall Average (All Plants)'),
    plt.Line2D([0], [0],
               color='#C73E1D', linestyle=':', linewidth=2.5,
               label='Hydropower Construction Year (Year 0)')
]

# Add Legend and Center It
ax_legend.legend(handles=legend_elements, loc='center', fontsize=12,
                 frameon=True, fancybox=True, shadow=True, framealpha=1.0)

# Save Legend Chart
save_path_legend = 'final_biomass_legend.png'
plt.tight_layout()
plt.savefig(save_path_legend, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Chart generation completed! Fixed tick font setting issue")
print(f"1. Trend chart saved to: {save_path_trend}")
print(f"2. Legend chart saved to: {save_path_legend}")
print(f"\n📊 Configuration Confirmation:")
print(f"   - X-axis ticks: 16pt font, bold")
print(f"   - Y-axis ticks: 16pt font, bold")
print(f"   - Data range: Starting from relative year -3, no n labels, including overall average line")
print(f"   - Save format: 300dpi high-resolution PNG, suitable for academic presentations")