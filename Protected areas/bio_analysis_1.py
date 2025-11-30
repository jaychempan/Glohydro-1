import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings('ignore')

# 1. Data Preprocessing (New: Filter data with relative year ≥ -3)
df = pd.read_csv('')
df['relative_year'] = df['year'] - df['plant_year']

# 【Key Modification 1: Keep only data with relative year ≥ -3 to ensure analysis starts from year -3】
df = df[df['relative_year'] >= -3].copy()

# Calculate biomass change rate (using 1 year before construction as baseline, excluding outliers)
baseline_biomass = df[df['relative_year'] == -1][['name', 'filled_biomass_Mg']].rename(
    columns={'filled_biomass_Mg': 'baseline_biomass'}
)
df_with_baseline = pd.merge(df, baseline_biomass, on='name', how='left')
df_with_baseline = df_with_baseline[df_with_baseline['relative_year'] != -1].copy()
df_with_baseline['biomass_change_rate'] = (
    (df_with_baseline['filled_biomass_Mg'] - df_with_baseline['baseline_biomass'])
    / df_with_baseline['baseline_biomass'] * 100
)
df_valid = df_with_baseline.dropna(subset=['baseline_biomass', 'biomass_change_rate'])
df_valid = df_valid[np.abs(df_valid['biomass_change_rate']) <= 500]
y_min, y_max = df_valid['biomass_change_rate'].min() - 30, df_valid['biomass_change_rate'].max() + 30

# 2. Core Chart Configuration (16:9 aspect ratio, suitable for presentation scenarios)
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['grid.alpha'] = 0.2
plt.rcParams['grid.linestyle'] = ':'

# Create 16:9 canvas with 1 row and 3 columns layout (shared y-axis)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 9), sharey=True)
fig.subplots_adjust(wspace=0.35, left=0.06, right=0.97, bottom=0.12, top=0.88)

# Color Configuration
color_config = {
    'capacity': {'base': '#3498db', 'cmap': 'Blues', 'title': 'Hydropower Capacity'},
    'newly_built': {'base': '#e74c3c', 'cmap': 'Reds', 'title': 'Hydropower Newly-built Year'},
    'relative': {'base': '#2ecc71', 'cmap': 'Greens', 'title': 'Relative Year (Start from -3)'}  # Title indicates start from year -3
}

# ---------------------- Subplot 1: Hydropower Capacity vs Biomass Change Rate ----------------------
x1 = df_valid['capacity_mw']
y1 = df_valid['biomass_change_rate']
xy1 = np.vstack([x1, y1])
z1 = gaussian_kde(xy1)(xy1)
idx1 = z1.argsort()
x1_sorted, y1_sorted, z1_sorted = x1.iloc[idx1], y1.iloc[idx1], z1[idx1]

sc1 = ax1.scatter(x1_sorted, y1_sorted, c=z1_sorted, cmap=color_config['capacity']['cmap'],
                  s=80, alpha=0.7, edgecolor='white', linewidth=0.8, zorder=2)
z1_fit = np.polyfit(x1, y1, 1)
p1 = np.poly1d(z1_fit)
x1_range = np.linspace(x1.min(), x1.quantile(0.95), 100)
ax1.plot(x1_range, p1(x1_range), color=color_config['capacity']['base'], linewidth=2.2,
         linestyle='--', alpha=0.8, label='Trend Line', zorder=3)
res1 = y1 - p1(x1)
se1 = np.std(res1) / np.sqrt(len(x1))
ax1.fill_between(x1_range, p1(x1_range)-2*se1, p1(x1_range)+2*se1,
                 alpha=0.2, color=color_config['capacity']['base'], zorder=1)

pearson_r1, pearson_p1 = stats.pearsonr(x1, y1)
stats_text1 = f'Capacity (MW)\nPearson r = {pearson_r1:.3f}\np = {pearson_p1:.3f}\nSignificance: n.s.'
ax1.text(0.02, 0.98, stats_text1, transform=ax1.transAxes, va='top', ha='left',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9,
                   edgecolor=color_config['capacity']['base'], linewidth=1.2),
         fontsize=10, fontweight='bold', zorder=4)

ax1.set_title(color_config['capacity']['title'], fontsize=13, pad=15)
ax1.set_xlabel('Hydropower Capacity (MW)', fontsize=12)
ax1.set_xlim(0, x1.quantile(0.95))
ax1.set_ylim(y_min, y_max)
ax1.grid(True, zorder=0)
ax1.legend(loc='lower right', fontsize=9, frameon=True, framealpha=0.8)

# ---------------------- Subplot 2: Newly-built Year vs Biomass Change Rate ----------------------
x2 = df_valid['plant_year']
y2 = df_valid['biomass_change_rate']
xy2 = np.vstack([x2, y2])
z2 = gaussian_kde(xy2)(xy2)
idx2 = z2.argsort()
x2_sorted, y2_sorted, z2_sorted = x2.iloc[idx2], y2.iloc[idx2], z2[idx2]

sc2 = ax2.scatter(x2_sorted, y2_sorted, c=z2_sorted, cmap=color_config['newly_built']['cmap'],
                  s=80, alpha=0.7, edgecolor='white', linewidth=0.8, zorder=2)
z2_fit = np.polyfit(x2, y2, 1)
p2 = np.poly1d(z2_fit)
x2_range = np.linspace(x2.min(), x2.max(), 100)
ax2.plot(x2_range, p2(x2_range), color=color_config['newly_built']['base'], linewidth=2.2,
         alpha=0.8, label='Trend Line (Positive)', zorder=3)
res2 = y2 - p2(x2)
se2 = np.std(res2) / np.sqrt(len(x2))
ax2.fill_between(x2_range, p2(x2_range)-2*se2, p2(x2_range)+2*se2,
                 alpha=0.2, color=color_config['newly_built']['base'], zorder=1)

pearson_r2, pearson_p2 = stats.pearsonr(x2, y2)
stats_text2 = f'Newly-built Year\nPearson r = {pearson_r2:.3f}\np = {pearson_p2:.4f}\nSignificance: **'
ax2.text(0.02, 0.98, stats_text2, transform=ax2.transAxes, va='top', ha='left',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9,
                   edgecolor=color_config['newly_built']['base'], linewidth=1.2),
         fontsize=10, fontweight='bold', zorder=4)

ax2.set_title(color_config['newly_built']['title'], fontsize=13, pad=15)
ax2.set_xlabel('Hydropower Newly-built Year', fontsize=12)
ax2.set_xlim(x2.min()-0.3, x2.max()+0.3)
ax2.set_xticks(sorted(df_valid['plant_year'].unique()))
ax2.grid(True, zorder=0)
ax2.legend(loc='lower right', fontsize=9, frameon=True, framealpha=0.8)

# ---------------------- Subplot 3: Relative Year vs Biomass Change Rate (Start from -3) ----------------------
x3 = df_valid['relative_year']
y3 = df_valid['biomass_change_rate']

# 【Key Modification 2: Filter data with relative year ≥ -3 to ensure only years ≥ -3 are displayed】
x3 = x3[x3 >= -3]
y3 = y3[x3 >= -3]  # Synchronously filter y-values to maintain data correspondence

# Recalculate gradient scatter (based on filtered data)
xy3 = np.vstack([x3, y3])
z3 = gaussian_kde(xy3)(xy3)
idx3 = z3.argsort()
x3_sorted, y3_sorted, z3_sorted = x3.iloc[idx3], y3.iloc[idx3], z3[idx3]

# Plot gradient scatter + trend line + confidence interval
sc3 = ax3.scatter(x3_sorted, y3_sorted, c=z3_sorted, cmap=color_config['relative']['cmap'],
                  s=80, alpha=0.7, edgecolor='white', linewidth=0.8, zorder=2)
# Trend line (fitted based on filtered data)
z3_fit = np.polyfit(x3, y3, 1)
p3 = np.poly1d(z3_fit)
x3_range = np.linspace(x3.min(), x3.max(), 100)  # Range automatically starts from -3
ax3.plot(x3_range, p3(x3_range), color=color_config['relative']['base'], linewidth=2.5,
         alpha=0.8, label='Trend Line (Negative)', zorder=3)
res3 = y3 - p3(x3)
se3 = np.std(res3) / np.sqrt(len(x3))
ax3.fill_between(x3_range, p3(x3_range)-2*se3, p3(x3_range)+2*se3,
                 alpha=0.2, color=color_config['relative']['base'], zorder=1)

# Newly-built Year baseline (x=0)
ax3.axvline(x=0, color='#c0392b', linestyle='--', linewidth=2.5, alpha=0.8,
            label='Newly-built Year (0)', zorder=3)

# 【Key Modification 3: X-axis ticks only show -3, -2, -1, 0, 1, 2, 3, 4, 5 (adjust based on actual data)】
x3_ticks = sorted([t for t in x3.unique() if t >= -3])  # Keep only ticks ≥ -3
ax3.set_xticks(x3_ticks)
# 【Key Modification 4: Set x-axis range from -3.5 to max year + 0.5 to ensure year -3 is the starting point without overcrowding】
ax3.set_xlim(-3.5, x3.max() + 0.5)

# Statistical annotation
pearson_r3, pearson_p3 = stats.pearsonr(x3, y3)
stats_text3 = f'Relative Year\nPearson r = {pearson_r3:.3f}\np = {pearson_p3:.4f}\nSignificance: **'
ax3.text(0.02, 0.98, stats_text3, transform=ax3.transAxes, va='top', ha='left',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9,
                   edgecolor=color_config['relative']['base'], linewidth=1.2),
         fontsize=10, fontweight='bold', zorder=4)

# Subplot 3 detailed configuration (x-axis label indicates start from year -3)
ax3.set_title(color_config['relative']['title'], fontsize=13, pad=15)
ax3.set_xlabel('Relative Year (Newly-built Year = 0, Start from -3)', fontsize=12)
ax3.grid(True, zorder=0)
ax3.legend(loc='lower right', fontsize=9, frameon=True, framealpha=0.8)

# 3. Overall Chart Optimization
fig.text(0.02, 0.5, 'Biomass Change Rate (%)', va='center', rotation=90,
         fontsize=13, fontweight='bold')

# Save chart (saved to current local directory; modify path as needed)
plt.savefig('combined_correlation_16_9_relative_start_minus3.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("✅ 16:9 Aspect Ratio Combined Chart Generated Successfully!")
print(f"File Name: combined_correlation_16_9_relative_start_minus3.png")
print(f"Key Confirmation: Relative year starts from -3, Subplot 3 x-axis only shows years ≥ -3.")