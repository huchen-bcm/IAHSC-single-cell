import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'

d = np.load(f'{OUT}/all_gene_tracks_coverage.npz')
meta = json.load(open(f'{OUT}/all_gene_tracks_meta.json'))

COL_NAIVE = '#5b7fa6'
COL_INFECTED = '#c0562f'
COL_GENE = '#2b2b2b'
COL_PEAK = '#f4c542'


def smooth(a, w=7):
    k = np.ones(w) / w
    return np.convolve(a, k, mode='same')


def plot_gene(gu):
    m = meta[gu]
    naive = d[f'{gu}__naive']
    infected = d[f'{gu}__infected']
    ws, we, binsize = m['win_start'], m['win_end'], m['bin']
    n = min(len(naive), len(infected), (we - ws) // binsize)
    naive, infected = naive[:n], infected[:n]
    x = (ws + np.arange(n) * binsize)[:n]
    naive_s = smooth(naive)
    infected_s = smooth(infected)
    ymax = max(naive_s.max(), infected_s.max()) * 1.15

    gname = m['display_name']
    strand = m['strand']
    gstart, gend = m['gene_start'], m['gene_end']
    exons = m['exons']
    peaks = m['peaks']
    chrom = m['chrom']

    fig = plt.figure(figsize=(9, 3.9))
    gs = fig.add_gridspec(4, 1, height_ratios=[2.0, 2.0, 0.85, 0.3], hspace=0.15)
    fig.subplots_adjust(top=0.80, bottom=0.10, left=0.11, right=0.97)

    ax_naive = fig.add_subplot(gs[0])
    ax_infected = fig.add_subplot(gs[1], sharex=ax_naive)
    ax_gene = fig.add_subplot(gs[2], sharex=ax_naive)
    ax_scale = fig.add_subplot(gs[3], sharex=ax_naive)

    for ax in (ax_naive, ax_infected, ax_gene):
        for p in peaks:
            ax.axvspan(p['start'], p['end'], color=COL_PEAK, alpha=0.45, zorder=0, lw=0)

    ax_naive.fill_between(x, 0, naive_s, color=COL_NAIVE, linewidth=0, alpha=0.9, zorder=2)
    ax_naive.set_ylim(0, ymax)
    ax_naive.set_ylabel('Naive', rotation=0, ha='right', va='center', fontsize=11, color=COL_NAIVE, fontweight='bold')
    ax_naive.text(0.99, 0.85, f'0–{ymax:.0f}', transform=ax_naive.transAxes, ha='right', va='top',
                  fontsize=8, color='#666666')

    ax_infected.fill_between(x, 0, infected_s, color=COL_INFECTED, linewidth=0, alpha=0.9, zorder=2)
    ax_infected.set_ylim(0, ymax)
    ax_infected.set_ylabel('Infected', rotation=0, ha='right', va='center', fontsize=11, color=COL_INFECTED,
                            fontweight='bold')
    ax_infected.text(0.99, 0.85, f'0–{ymax:.0f}', transform=ax_infected.transAxes, ha='right', va='top',
                      fontsize=8, color='#666666')

    for ax in (ax_naive, ax_infected):
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.set_yticks([])
        ax.tick_params(axis='x', labelbottom=False, bottom=False)
        ax.spines['bottom'].set_color('#cccccc')

    # peak stat annotations, placed above the naive track
    for p in peaks:
        mid = (p['start'] + p['end']) / 2
        ax_naive.annotate(f"log2FC={p['log2fc']:.2f}\npadj={p['padj']:.3f}",
                           xy=(mid, ymax), xytext=(0, 2), textcoords='offset points',
                           ha='center', va='bottom', fontsize=6.5, color='#8a6d1a')

    # --- gene model ---
    ax_gene.set_xlim(ws, we)
    ax_gene.set_ylim(-1, 1)
    ax_gene.axis('off')
    ax_gene.plot([gstart, gend], [0, 0], color=COL_GENE, lw=1.2, zorder=2)

    arrow_xs = np.linspace(gstart + (gend - gstart) * 0.08, gend - (gend - gstart) * 0.08, max(3, int((gend-gstart)/8000)))
    dx = (we - ws) * 0.006
    for ax_x in arrow_xs:
        if strand == '+':
            ax_gene.annotate('', xy=(ax_x + dx, 0), xytext=(ax_x - dx, 0),
                              arrowprops=dict(arrowstyle='-|>', color=COL_GENE, lw=1.0), zorder=3)
        else:
            ax_gene.annotate('', xy=(ax_x - dx, 0), xytext=(ax_x + dx, 0),
                              arrowprops=dict(arrowstyle='-|>', color=COL_GENE, lw=1.0), zorder=3)

    exon_h = 0.32
    for start, end in exons:
        w = max(end - start, (we - ws) * 0.0015)
        ax_gene.add_patch(Rectangle((start, -exon_h / 2), w, exon_h, facecolor=COL_GENE, edgecolor='none', zorder=4))

    strand_lbl = '+ strand' if strand == '+' else '− strand'
    ax_gene.annotate(gname, xy=(ws, 0.85), xycoords='data', xytext=(0, 0), textcoords='offset points',
                      fontsize=11, style='italic', fontweight='bold', va='bottom', ha='left', color=COL_GENE)
    ax_gene.annotate(f'  ({strand_lbl})', xy=(ws, 0.85), xycoords='data', xytext=(len(gname)*6+6, 0),
                      textcoords='offset points', fontsize=8, va='bottom', ha='left', color='#666666')

    # --- scale bar ---
    ax_scale.axis('off')
    ax_scale.set_xlim(ws, we)
    span = we - ws
    if span < 15000:
        bar_len = 2000
    elif span < 60000:
        bar_len = 10000
    else:
        bar_len = 20000
    bx1 = we - span * 0.02
    bx0 = bx1 - bar_len
    ax_scale.plot([bx0, bx1], [0.5, 0.5], color='#333333', lw=1.5)
    ax_scale.text((bx0 + bx1) / 2, 0.05, f'{bar_len // 1000} kb', ha='center', va='top', fontsize=8, color='#333333')
    ax_scale.text(ws, 0.05, f'{chrom}:{ws:,}-{we:,} (mm10)', ha='left', va='top', fontsize=8, color='#666666')

    fig.suptitle(f'Chromatin accessibility at {gname} — naive vs. infected IA-HSC (scATAC-seq)',
                 fontsize=12.5, fontweight='bold', y=0.995)
    fig.text(0.5, 0.895,
              'Yellow band(s): SnapATAC2 diff_test peak(s) linked to this gene (±100kb); tracks = fragment-end pileups per 1,000 cells',
              ha='center', fontsize=8, color='#666666', style='italic')

    fname = gname.replace('-', '')
    plt.savefig(f'{OUT}/assets/genome_track_{fname}.png', dpi=200, facecolor='white')
    plt.close(fig)


for gu in meta:
    plot_gene(gu)
    print('plotted', meta[gu]['display_name'])
