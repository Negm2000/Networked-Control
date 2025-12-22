import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
os.chdir("Project")
# --- 1. Parsing Logic ---
def get_file_content(structure):
    fname = f'log_{structure}.txt'
    if structure == 'Distributed':
        fname = 'log_Distributed_2to1.txt'
    with open(fname, 'r') as f:
        return f.read()

STRUCTURES = ['Centralized','Distributed', 'Decentralized']
KS = [0.2, 2.0, 200.0]
METHODS_ORDER = ['Stability LMI', 'Speed LMI', 'Damping LMI', 'Multi-Objective', 'H2']
rows = []
for s in STRUCTURES:
    for k in KS:
        for m in METHODS_ORDER:
            rows.append({'Structure': s, 'Coupling_K': k, 'Method': m})
df = pd.DataFrame(rows)

for idx, row in df.iterrows():
    structure = row['Structure']
    k_target = row['Coupling_K']
    method = row['Method']
    
    content = get_file_content(structure)
    k_str = f"{k_target:.3f}"
    pattern_k = r'RUNNING FOR K = ' + re.escape(k_str) + r'[\s\S]*?(?=RUNNING FOR K = |$)'
    k_match = re.search(pattern_k, content)
    
    if not k_match:
        continue
    k_block = k_match.group(0)
    
    if method == 'H-Infinity':
        method_regex_name = r'H-inf'
    elif method == 'Multi-Objective':
        method_regex_name = r'Multi-Objective'
    else:
        method_regex_name = re.escape(method)

    success_pattern = (
        r'Check .*?\(' + method_regex_name + r'\)[\s\S]*?'
        r'The spectral radius is ([\d\.]+) < 1 \(Stable\)[\s\S]*?'
        r'Theoretical Max Settling Time \(2%\) = ([\d\.]+)[\s\S]*?'
        r'Controller Gain Norm \|K\|[\s\S]*?= ([\d\.]+)'
    )
    
    m_match = re.search(success_pattern, k_block)
    
    if m_match:
        df.at[idx, 'Status'] = 'Stable'
        df.at[idx, 'Spectral_Radius'] = float(m_match.group(1))
        df.at[idx, 'Theo_Settling_Time'] = float(m_match.group(2))
        df.at[idx, 'Gain_Norm'] = float(m_match.group(3))
    else:
        df.at[idx, 'Status'] = 'Infeasible'
        df.at[idx, 'Spectral_Radius'] = np.nan
        df.at[idx, 'Theo_Settling_Time'] = np.nan
        df.at[idx, 'Gain_Norm'] = np.nan

# --- 2. Plotting Logic (Separate Files per K) ---
sns.set_theme(style="whitegrid", context="talk")
palette_dict = {
    'Centralized': '#2ecc71', # Green
    'Distributed': '#3498db', # Blue
    'Decentralized': '#e74c3c' # Red
}

generated_files = []

for k_val in KS:
    # Create a separate figure for each K
    # Increased figsize height to allow more space for titles
    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    
    # Filter data
    df_k = df[(df['Coupling_K'] == k_val) & (df['Status'] == 'Stable')].copy()
    
    if df_k.empty:
        if k_val == 200.0:
            pass 
        else:
            plt.close(fig)
            continue

    # Plot 1: Gain Norm
    ax1 = axes[0]
    sns.barplot(data=df_k, x='Method', y='Gain_Norm', hue='Structure', palette=palette_dict, ax=ax1,
                edgecolor='white', linewidth=1, alpha=0.9, order=METHODS_ORDER)
    
    ax1.set_ylabel('Gain Norm |K|')
    ax1.set_title('Controller Effort (|K|)', fontweight='bold', pad=20) # Added padding
    
    # Plot 2: Spectral Radius
    ax2 = axes[1]
    sns.barplot(data=df_k, x='Method', y='Spectral_Radius', hue='Structure', palette=palette_dict, ax=ax2,
                edgecolor='white', linewidth=1, alpha=0.9, order=METHODS_ORDER)
    ax2.axhline(1.0, color='#e74c3c', linestyle='--', linewidth=2, label='Unstable')
    ax2.set_ylabel('Spectral Radius')
    ax2.set_title('Stability Margin', fontweight='bold', pad=20) # Added padding
    ax2.set_ylim(0, 1.1)
    
    # Plot 3: Settling Time (New)
    ax3 = axes[2]
    sns.barplot(data=df_k, x='Method', y='Theo_Settling_Time', hue='Structure', palette=palette_dict, ax=ax3,
                edgecolor='white', linewidth=1, alpha=0.9, order=METHODS_ORDER)
    ax3.set_ylabel('Time (s)')
    ax3.set_title('Theoretical Max Settling Time', fontweight='bold', pad=20) # Added padding
    
    # Common Formatting
    for ax in axes:
        ax.set_xlabel('Control Method', labelpad=15)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, axis='y', alpha=0.5)
        if ax.get_legend():
            ax.get_legend().remove()
            
    # Add a single legend for the whole figure
    # Adjusted bbox_to_anchor to push legend higher up so it doesn't overlap title
    handles, labels = ax3.get_legend_handles_labels()
    fig.legend(handles, labels, title='Control Structure', loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, frameon=False, fontsize=16, title_fontsize=18)
    
    # Adjust layout with specific rect to leave space at top for legend
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Filename
    filename = f'metrics_K{k_val}.png'
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    generated_files.append(filename)
    plt.close(fig) 

print(generated_files)


rows = []
for s in STRUCTURES:
    for k in KS:
        for m in METHODS_ORDER:
            rows.append({'Structure': s, 'Coupling_K': k, 'Method': m})
df = pd.DataFrame(rows)

for idx, row in df.iterrows():
    structure = row['Structure']
    k_target = row['Coupling_K']
    method = row['Method']
    
    content = get_file_content(structure)
    k_str = f"{k_target:.3f}"
    pattern_k = r'RUNNING FOR K = ' + re.escape(k_str) + r'[\s\S]*?(?=RUNNING FOR K = |$)'
    k_match = re.search(pattern_k, content)
    
    if not k_match:
        continue
    k_block = k_match.group(0)
    
    if method == 'H-Infinity':
        method_regex_name = r'H-inf'
    elif method == 'Multi-Objective':
        method_regex_name = r'Multi-Objective'
    else:
        method_regex_name = re.escape(method)

    success_pattern = (
        r'Check .*?\(' + method_regex_name + r'\)[\s\S]*?'
        r'The spectral radius is ([\d\.]+) < 1 \(Stable\)[\s\S]*?'
    )
    
    m_match = re.search(success_pattern, k_block)
    
    if m_match:
        df.at[idx, 'Status'] = 'Stable'
    else:
        df.at[idx, 'Status'] = 'Infeasible'

# --- 2. Heatmap Generation Only ---
sns.set_theme(style="white", context="talk")

# Create a numeric column for heatmap coloring (1=Pass, 0=Fail)
df['Status_Numeric'] = df['Status'].apply(lambda x: 1 if x == 'Stable' else 0)

# Create the label column for the Y-axis
df['Config'] = df['Structure'] + " (k=" + df['Coupling_K'].astype(str) + ")"

# Pivot the table
# Index: Config, Columns: Method, Values: Status_Numeric
heatmap_data = df.pivot_table(index='Config', columns='Method', values='Status_Numeric', aggfunc='first') # Use 'first' or 'max', doesn't matter as values are unique

# Re-order the Index (Rows) to match Structure order then K order
# We want Decentralized -> Distributed -> Centralized (or alphabetical, but let's stick to logical groupings)
# Let's create a sort key
# Actually, the user's logs have: Decentralized, Distributed, Centralized.
# Let's sort by Structure then K.
df['Structure_Cat'] = pd.Categorical(df['Structure'], categories=['Centralized', 'Distributed', 'Decentralized'], ordered=True)
df = df.sort_values(['Structure_Cat', 'Coupling_K'])
sorted_index = df['Config'].unique()
heatmap_data = heatmap_data.reindex(sorted_index)

# Re-order the Columns (Methods) to match User's Order
heatmap_data = heatmap_data.reindex(columns=METHODS_ORDER)

# Plotting
fig, ax = plt.subplots(figsize=(14, 8))
cmap = sns.color_palette(['#e74c3c', '#2ecc71']) # Red, Green

sns.heatmap(heatmap_data, cmap=cmap, cbar=False, linewidths=2, linecolor='white', square=False, annot=True, fmt='g', 
            annot_kws={"size": 14, "weight": "bold"}, ax=ax)

# Replace 0/1 with Fail/Pass text
for text in ax.texts:
    try:
        val = float(text.get_text())
        text.set_text('Pass' if val == 1 else 'Fail')
        text.set_color('white')
    except:
        pass

ax.set_title('Feasibility Map: Controller Survival', pad=20, fontsize=20, fontweight='bold', color='#333333')
ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='x', rotation=30)
ax.tick_params(axis='y', rotation=0)

plt.tight_layout()
plt.savefig('feasibility_map_only.png', dpi=300)
print("Feasibility map generated")