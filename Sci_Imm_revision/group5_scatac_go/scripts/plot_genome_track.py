import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import matplotlib.font_manager as fm

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'
GENE = 'Irf1'

d = np.load(f'{OUT}/genome_track_{GENE}_coverage.npz')
naive, infected = d['naive'], d['infected']
win_start, win_end, binsize = int(d['win_start']), int(d['win_end']), int(d['bin'])
x = np.arange(win_start, win_end, binsize)[:len(naive)]

exons = pd.read_csv(f'{OUT}/genome_track_{GENE}_exons.csv')
with open(f'{OUT}/genome_track_{GENE}_meta.txt') as f:
    chrom, ws, we, strand, gstart, gend = f.read().strip().split('\t')
gstart, gend = int(gstart), int(gend)

# light smoothing for readability (rolling mean, matches typical browser-track rendering)
def smooth(a, w=5):
    k = np.ones(w)/w
    return np.convolve(a, k, mode='same')

naive_s = smooth(naive, 7)
infected_s = smooth(infected, 7)
ymax = max(naive_s.max(), infected_s.max()) * 1.15

COL_NAIVE = '#5b7fa6'
COL_INFECTED = '#c0562f'
COL_GENE = '#2b2b2b'
COL_PEAK = '#d9d9d9'

fig = plt.figure(figsize=(9, 3.9))
gs = fig.add_gridspec(4, 1, height_ratios=[2.0, 2.0, 0.85, 0.3], hspace=0.15)
fig.subplots_adjust(top=0.80, bottom=0.10, left=0.11, right=0.97)

ax_naive = fig.add_subplot(gs[0])
ax_infected = fig.add_subplot(gs[1], sharex=ax_naive)
ax_gene = fig.add_subplot(gs[2], sharex=ax_naive)
ax_scale = fig.add_subplot(gs[3], sharex=ax_naive)

# --- peak-highlight shading: TSS/promoter region of Irf1 (+ strand -> TSS at gstart) ---
promoter_start, promoter_end = gstart - 800, gstart + 300
for ax in (ax_naive, ax_infected, ax_gene):
    ax.axvspan(promoter_start, promoter_end, color=COL_PEAK, alpha=0.6, zorder=0, lw=0)

# --- coverage tracks ---
ax_naive.fill_between(x, 0, naive_s, color=COL_NAIVE, linewidth=0, alpha=0.9, zorder=2)
ax_naive.set_ylim(0, ymax)
ax_naive.set_ylabel('Naive', rotation=0, ha='right', va='center', fontsize=11, color=COL_NAIVE, fontweight='bold')
ax_naive.text(0.99, 0.85, f'0–{ymax:.0f}', transform=ax_naive.transAxes, ha='right', va='top',
              fontsize=8, color='#666666')

ax_infected.fill_between(x, 0, infected_s, color=COL_INFECTED, linewidth=0, alpha=0.9, zorder=2)
ax_infected.set_ylim(0, ymax)
ax_infected.set_ylabel('Infected', rotation=0, ha='right', va='center', fontsize=11, color=COL_INFECTED, fontweight='bold')
ax_infected.text(0.99, 0.85, f'0–{ymax:.0f}', transform=ax_infected.transAxes, ha='right', va='top',
                  fontsize=8, color='#666666')

for ax in (ax_naive, ax_infected):
    ax.spines[['top','right','left']].set_visible(False)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelbottom=False, bottom=False)
    ax.spines['bottom'].set_color('#cccccc')

# --- gene model track ---
ax_gene.set_xlim(win_start, win_end)
ax_gene.set_ylim(-1, 1)
ax_gene.axis('off')

# intron line
ax_gene.plot([gstart, gend], [0, 0], color=COL_GENE, lw=1.2, zorder=2)

# strand arrows along the intron line
arrow_xs = np.linspace(gstart + (gend-gstart)*0.08, gend - (gend-gstart)*0.08, 7)
for ax_x in arrow_xs:
    dx = (win_end - win_start) * 0.006
    if strand == '+':
        ax_gene.annotate('', xy=(ax_x+dx, 0), xytext=(ax_x-dx, 0),
                          arrowprops=dict(arrowstyle='-|>', color=COL_GENE, lw=1.0), zorder=3)
    else:
        ax_gene.annotate('', xy=(ax_x-dx, 0), xytext=(ax_x+dx, 0),
                          arrowprops=dict(arrowstyle='-|>', color=COL_GENE, lw=1.0), zorder=3)

# exon boxes
exon_h = 0.32
for _, row in exons.iterrows():
    w = max(row['end'] - row['start'], (win_end-win_start)*0.0015)
    rect = Rectangle((row['start'], -exon_h/2), w, exon_h, facecolor=COL_GENE, edgecolor='none', zorder=4)
    ax_gene.add_patch(rect)

strand_lbl = '+ strand' if strand == '+' else '− strand'
ax_gene.annotate(GENE, xy=(win_start, 0.85), xycoords='data', xytext=(0, 0), textcoords='offset points',
                  fontsize=11, style='italic', fontweight='bold', va='bottom', ha='left', color=COL_GENE)
ax_gene.annotate(f'  ({strand_lbl})', xy=(win_start, 0.85), xycoords='data', xytext=(48, 0), textcoords='offset points',
                  fontsize=8, va='bottom', ha='left', color='#666666')

# --- scale bar / axis ---
ax_scale.axis('off')
ax_scale.set_xlim(win_start, win_end)
span = win_end - win_start
bar_len = 2000 if span < 30000 else 5000
bx0 = win_end - span*0.02 - bar_len
bx1 = bx0 + bar_len
ax_scale.plot([bx0, bx1], [0.5, 0.5], color='#333333', lw=1.5)
ax_scale.text((bx0+bx1)/2, 0.05, f'{bar_len//1000} kb', ha='center', va='top', fontsize=8, color='#333333')
ax_scale.text(win_start, 0.05, f'{chrom}:{win_start:,}-{win_end:,} (mm10)', ha='left', va='top',
              fontsize=8, color='#666666')

fig.suptitle(f'Chromatin accessibility at {GENE} — naive vs. infected IA-HSC (scATAC-seq)',
             fontsize=12.5, fontweight='bold', y=0.995)
fig.text(0.5, 0.895, 'Grey band: promoter/TSS region; tracks are fragment-end pileups normalized per 1,000 cells',
          ha='center', fontsize=8.5, color='#666666', style='italic')

plt.savefig(f'{OUT}/assets/genome_track_{GENE}.png', dpi=200, bbox_inches=None, facecolor='white')
print('saved figure')
