import matplotlib.pyplot as plt
import numpy as np

def digest_dna(sequence, enzymes):
    """
    Simuleaza digestia enzimatica completa a unei secvente de ADN.
    Returneaza dimensiunile fragmentelor rezultate (in bp).
    """
    cut_positions = set()
    for name, info in enzymes.items():
        pattern = info["site"]
        cut_offset = info["cut_offset"]
        start = 0
        while True:
            idx = sequence.find(pattern, start)
            if idx == -1:
                break
            cut_positions.add(idx + cut_offset)
            start = idx + 1

    sorted_cuts = sorted(list(cut_positions))
    fragments = []
    prev = 0
    for cut in sorted_cuts:
        fragments.append(len(sequence[prev:cut]))
        prev = cut
    fragments.append(len(sequence[prev:]))
    return sorted(fragments, reverse=True)


def draw_gel(title, sample_name, sample_fragments, ladder_bands, band_color="#ff77ff"):
    """
    Deseneaza un gel individual de electroforeza in agaroza.
    """
    fig, ax = plt.subplots(figsize=(5, 7.5), facecolor="#1e1e1e")
    ax.set_facecolor("#111111")

    def size_to_y(bp):
        return np.log10(bp)

    lanes = {
        1: {"name": "DNA Ladder", "fragments": ladder_bands, "color": "#00ffcc"},
        2: {"name": sample_name, "fragments": sample_fragments, "color": band_color}
    }

    for lane_idx, data in lanes.items():
        ax.add_patch(plt.Rectangle((lane_idx - 0.25, 3.25), 0.5, 0.05, color="#444444"))
        
        for bp in data["fragments"]:
            y_pos = size_to_y(bp)
            ax.plot([lane_idx - 0.23, lane_idx + 0.23], [y_pos, y_pos], 
                    color=data["color"], lw=3.8, solid_capstyle='round')
            
            if lane_idx == 1:
                ax.text(lane_idx - 0.28, y_pos, f"{bp} bp", color="#cccccc", 
                        fontsize=8, va="center", ha="right")
            else:
                ax.text(lane_idx + 0.28, y_pos, f"{bp} bp", color="#ffffff", 
                        fontsize=8, va="center", ha="left")

    ax.set_xlim(0.3, 2.7)
    ax.set_ylim(np.log10(70), 3.4)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([lanes[1]["name"], lanes[2]["name"]], color="white", fontsize=9, fontweight="bold")
    ax.set_yticks([])

    ax.text(0.4, 3.27, "(-) Catod", color="#ff5555", fontsize=9, fontweight="bold")
    ax.text(0.4, np.log10(80), "(+) Anod", color="#55ff55", fontsize=9, fontweight="bold")
    ax.set_title(title, color="white", fontsize=12, pad=15, fontweight="bold")

    plt.tight_layout()
    return fig


genome_1 = (
    "ATGGCCTACGTACCGTTAGGCTACGTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGC"
    "TAGCTAGCTAGCTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "ATCGATCGATCGATCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGGATCCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGAAGCTTGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCTAA"
)

genome_2 = (
    "ATGGCCAAGCTTGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "ATCGATCGATCGATCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "GATCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCTAA"
)

enzymes = {
    "EcoRI":   {"site": "GAATTC", "cut_offset": 1},
    "BamHI":   {"site": "GGATCC", "cut_offset": 1},
    "HindIII": {"site": "AAGCTT", "cut_offset": 1}
}

fragments_1 = digest_dna(genome_1, enzymes)
fragments_2 = digest_dna(genome_2, enzymes)
ladder_bands = [1500, 1000, 750, 500, 400, 300, 200, 100]

fig1 = draw_gel(
    title="Gel 1: Digestie Genom 1",
    sample_name="Genom 1 Digest\n(EcoRI + BamHI + HindIII)",
    sample_fragments=fragments_1,
    ladder_bands=ladder_bands,
    band_color="#ff77ff"
)

fig2 = draw_gel(
    title="Gel 2: Digestie Genom 2",
    sample_name="Genom 2 Digest\n(EcoRI + BamHI + HindIII)",
    sample_fragments=fragments_2,
    ladder_bands=ladder_bands,
    band_color="#ffaa00"
)

plt.show()