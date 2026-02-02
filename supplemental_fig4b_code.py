import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']

df = pd.read_csv('merged_proteomics_data.csv')

df = df[~df['PreyGene'].str.startswith('Cont', na=False)].copy()

df = df[df['PreyGene'] != 'TXP0CG47'].copy()

se_wt_ub = df['SE_UB_WT_BR_SaintScore']
se_wt_ub_mg = df['SE_UB_WT_BR_MG_SaintScore']

df_filtered = df[(se_wt_ub > 0.65) | (se_wt_ub_mg > 0.65)].copy()

conditions = [
    ('WT/WT-Ub', 'WT_UB_WT_BR_SaintScore'),
    ('WT/WT-Ub+MG', 'WT_UB_WT_BR_MG_SaintScore'),
    ('WT/K0-Ub', 'WT_KO_BR_SaintScore'),
    ('WT/K0-Ub+MG', 'WT_KO_BR_MG_SaintScore'),
    ('SE/WT-Ub', 'SE_UB_WT_BR_SaintScore'),
    ('SE/WT-Ub+MG', 'SE_UB_WT_BR_MG_SaintScore'),
    ('SE/K0-Ub', 'SE_KO_BR_SaintScore'),
    ('SE/K0-Ub+MG', 'SE_KO_BR_MG_SaintScore'),
]

saint_matrix = pd.DataFrame()
for name, col in conditions:
    saint_matrix[name] = df_filtered[col].values
saint_matrix.index = df_filtered['PreyGene'].values

translation_rqc = ['RPL5', 'EIF4B', 'PABPC1', 'PABPC3', 'PABPC4', 'UPF1', 
                   'RPL13', 'RPL18', 'RPS8', 'RPS27A', 'EIF6', 'HNRNPA2B1', 
                   'HNRNPA1', 'PEG10', 'SRP72']

ups_proteostasis = ['VCP', 'PSMC1', 'PSMC2', 'PSMC5', 'PSMC6', 'PSMA1', 'PSMA6',
                    'PSMB5', 'PSMB6', 'SQSTM1', 'HSP90B1', 'HSP90AA1', 'HSPB1',
                    'HSPA1A', 'CCT8', 'FKBP1A', 'CSNK2A1', 'VBL', 'THOC1',
                    'CYLD', 'SPATA2']

cytoskeletal_membrane = ['MYO6', 'FLNA', 'FLNC', 'MAP4', 'MAP1B', 'MAP2', 
                         'MYH9', 'TLN1', 'CTTN', 'ANXA1', 'ANXA2', 'ACTN4',
                         'TUBB4B', 'SPTA2', 'DPYSL2', 'DPYSL3', 'DPYSL4',
                         'MYOF', 'EZR', 'MYL9', 'MYL12A', 'IGF2R', 'DES', 'TAGLN2']

other_regulators = ['AHNAK', 'IFI16', 'DESI1', 'MVP', 'DSG1', 'ALDH18A1']

ordered_proteins = []
for cat_list in [translation_rqc, ups_proteostasis, cytoskeletal_membrane, other_regulators]:
    for p in cat_list:
        if p in saint_matrix.index and p not in ordered_proteins:
            ordered_proteins.append(p)

remaining = [p for p in saint_matrix.index if p not in ordered_proteins]
if remaining:
    se_scores = saint_matrix.loc[remaining, 'SE/WT-Ub'].fillna(0) + \
                saint_matrix.loc[remaining, 'SE/WT-Ub+MG'].fillna(0)
    remaining_sorted = se_scores.sort_values(ascending=False).index.tolist()
    ordered_proteins.extend(remaining_sorted)

saint_subset = saint_matrix.loc[ordered_proteins]

# Calculate figure size for square cells (compact layout)
n_rows, n_cols = saint_subset.shape
cell_size = 0.3  # size of each cell in inches
fig_width = n_cols * cell_size + 3  # extra space for labels and colorbar
fig_height = n_rows * cell_size + 2  # extra space for title and x-labels

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

colors = ['#deebf7', '#c6dbef', '#9ecae1', '#ffffff', '#fee0d2', '#fc9272', '#de2d26', '#a50f15', '#67000d']
cmap = LinearSegmentedColormap.from_list('blue_white_red', colors, N=100)

im = ax.imshow(saint_subset.values, cmap=cmap, aspect='equal',
               vmin=0, vmax=1, interpolation='nearest')

for i in range(len(saint_subset) + 1):
    ax.axhline(i - 0.5, color='white', linewidth=0.5)
for j in range(len(saint_subset.columns) + 1):
    ax.axvline(j - 0.5, color='white', linewidth=0.5)

ax.set_xticks(range(len(saint_subset.columns)))
ax.set_xticklabels(saint_subset.columns, fontsize=8, rotation=45, ha='right')
ax.xaxis.tick_bottom()

ax.set_yticks(range(len(saint_subset)))
ax.set_yticklabels(saint_subset.index, fontsize=7)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)
    spine.set_edgecolor('black')

cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('SAINT Score', fontsize=11)
cbar.outline.set_linewidth(1.5)

ax.set_title('Split-TurboID: WT-Ub vs K0-Ub Control\nDeSI1 WT and SE Interactions', 
             fontsize=13, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig('supplemental_fig4b_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Heatmap generated: {len(saint_subset)} proteins")
print(f"Filtered for SE DeSI1 + WT-Ub SAINT > 0.65")
print(f"Translation/RQC: rows 0-6 (7 proteins)")
print(f"UPS/Proteostasis: rows 7-16 (10 proteins)")
print(f"Cytoskeletal/Membrane: rows 17-32 (16 proteins)")
print(f"Regulators: rows 33-34 (2 proteins)")
