import gzip
import math
import random
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

LADDER_SIZES = [3000, 2000, 1500, 1200, 1000, 800, 600, 400, 200, 100]


def parse_ncbi_fasta(file_path, target_length=2500):
    """Citește fișiere .fna, .fa, .fasta sau variante arhivate .gz de pe NCBI

    și extrage o secvență continuă de ADN de lungime între 1000 și 3000 bp.
    """
    open_func = (
        gzip.open
        if file_path.endswith(".gz")
        else open
    )

    sequence_chars = []
    with open_func(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            cleaned_line = "".join(ch for ch in line.upper() if ch in "ATCGN")
            sequence_chars.append(cleaned_line)

            if sum(len(s) for s in sequence_chars) >= target_length:
                break

    full_seq = "".join(sequence_chars)
    return full_seq[:target_length]


class GelElectrophoresisApp:

    def __init__(self, root):
        self.root = root
        self.root.title("NCBI Gel Electrophoresis Simulator")
        self.root.geometry("980x680")
        self.root.configure(bg="#1e1e1e")

        self.current_dna_seq = ""
        self.samples = []  

        self.setup_ui()

    def setup_ui(self):
        top_bar = tk.Frame(self.root, bg="#2d2d2d", padx=10, pady=10)
        top_bar.pack(side=tk.TOP, fill=tk.X)

        self.file_btn = ttk.Button(
            top_bar,
            text="Încarcă fișier NCBI (.fna.gz / .fna)",
            command=self.load_file,
        )
        self.file_btn.pack(side=tk.LEFT, padx=5)

        self.run_btn = ttk.Button(
            top_bar,
            text="Generează 10 Fragmente & Rulează Gel",
            command=self.run_simulation,
            state=tk.DISABLED,
        )
        self.run_btn.pack(side=tk.LEFT, padx=5)

        self.status_lbl = tk.Label(
            top_bar,
            text="Selectează fișierul descărcat de pe NCBI...",
            fg="#cccccc",
            bg="#2d2d2d",
            font=("Arial", 10),
        )
        self.status_lbl.pack(side=tk.LEFT, padx=15)

        main_pane = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, bg="#1e1e1e")
        main_pane.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.canvas_width = 520
        self.canvas_height = 560
        self.canvas = tk.Canvas(
            main_pane,
            width=self.canvas_width,
            height=self.canvas_height,
            bg="#0b1118",
            highlightthickness=0,
        )
        main_pane.add(self.canvas)

        right_frame = tk.Frame(main_pane, bg="#252526", padx=10, pady=10)
        main_pane.add(right_frame)

        tk.Label(
            right_frame,
            text="Fragmente Extrase în Array (10 Probe):",
            fg="#4ec9b0",
            bg="#252526",
            font=("Arial", 11, "bold"),
        ).pack(anchor="w", pady=(0, 5))

        cols = ("sample", "length", "start_pos")
        self.tree = ttk.Treeview(
            right_frame, columns=cols, show="headings", height=12
        )
        self.tree.heading("sample", text="Pistă (Lane)")
        self.tree.heading("length", text="Mărime (bp)")
        self.tree.heading("start_pos", text="Start Index")

        self.tree.column("sample", width=90, anchor="center")
        self.tree.column("length", width=90, anchor="center")
        self.tree.column("start_pos", width=90, anchor="center")
        self.tree.pack(fill=tk.BOTH, expand=True)

        legend_text = (
            "Principiul Electrogenezei pe Gel:\n\n"
            "• Godeurile sunt la catod (-), ADN-ul migrează spre anod (+).\n"
            "• Fragmentele scurte întâmpină forțe de frecare reduse\n"
            "  și ajung cel mai departe pe gel.\n"
            "• Fragmentele lungi au frecare mare și migrează puțin."
        )
        tk.Label(
            right_frame,
            text=legend_text,
            fg="#9cdcfe",
            bg="#252526",
            justify=tk.LEFT,
            font=("Consolas", 9),
        ).pack(fill=tk.X, pady=10)

    def load_file(self):
        file_path = filedialog.askopenfilename(
            title="Alege fișierul NCBI",
            filetypes=[
                (
                    "Fișiere NCBI/FASTA",
                    "*.fna.gz *.fasta.gz *.fa.gz *.fna *.fasta *.fa",
                ),
                ("Toate fișierele", "*.*"),
            ],
        )
        if not file_path:
            return

        try:
            target_len = random.randint(1800, 2600)
            self.current_dna_seq = parse_ncbi_fasta(
                file_path, target_length=target_len
            )

            if len(self.current_dna_seq) < 500:
                messagebox.showerror(
                    "Eroare",
                    "Fișierul nu conține suficiente nucleotide valide.",
                )
                return

            self.status_lbl.config(
                text=f"Secvență extrasă: {len(self.current_dna_seq)} nt din fișier.",
                fg="#4ec9b0",
            )
            self.run_btn.config(state=tk.NORMAL)
            self.run_simulation()

        except Exception as e:
            messagebox.showerror(
                "Eroare la citire", f"Nu s-a putut citi fișierul:\n{e}"
            )

    def run_simulation(self):
        self.samples.clear()
        total_len = len(self.current_dna_seq)

        for i in range(10):
            frag_len = random.randint(100, total_len)
            max_start = total_len - frag_len
            start_idx = random.randint(0, max_start) if max_start > 0 else 0

            self.samples.append(
                {
                    "lane": f"Pista {i+1}",
                    "length": frag_len,
                    "start": start_idx,
                }
            )

        self.update_table()
        self.draw_gel()

    def update_table(self):
        for row in self.tree.get_children():
            self.tree.delete(row)
        for s in self.samples:
            self.tree.insert(
                "",
                tk.END,
                values=(s["lane"], f"{s['length']} bp", f"nt {s['start']}"),
            )

    def calculate_y(self, bp, top_y, bottom_y):
        """Calcul pe baza legii migrației: distanța este invers proporțională

        cu logaritmul mărimii (bp).
        """
        min_bp, max_bp = 100, 3000
        log_min = math.log10(min_bp)
        log_max = math.log10(max_bp)
        log_bp = math.log10(max(min_bp, min(max_bp, bp)))

        fraction = (log_max - log_bp) / (log_max - log_min)
        return top_y + fraction * (bottom_y - top_y)

    def draw_gel(self):
        self.canvas.delete("all")

        gx1, gy1 = 45, 30
        gx2, gy2 = self.canvas_width - 25, self.canvas_height - 35
        self.canvas.create_rectangle(
            gx1, gy1, gx2, gy2, fill="#121e2b", outline="#29465b", width=2
        )

        self.canvas.create_text(
            22,
            gy1 + 15,
            text="(-)",
            fill="#e06c75",
            font=("Arial", 11, "bold"),
        )
        self.canvas.create_text(
            22,
            gy2 - 15,
            text="(+)",
            fill="#98c379",
            font=("Arial", 11, "bold"),
        )

        total_lanes = 11  
        lane_step = (gx2 - gx1) / (total_lanes + 1)
        well_y = gy1 + 25
        well_h = 10
        band_w = 22

        for i in range(total_lanes):
            cx = gx1 + (i + 1) * lane_step
            self.canvas.create_rectangle(
                cx - band_w / 2,
                well_y - well_h,
                cx + band_w / 2,
                well_y,
                fill="#000000",
                outline="#3e5366",
            )
            label = "Ladder" if i == 0 else str(i)
            self.canvas.create_text(
                cx,
                well_y - well_h - 10,
                text=label,
                fill="#8fa1b3",
                font=("Arial", 8),
            )

        ladder_x = gx1 + 1 * lane_step
        for size in LADDER_SIZES:
            y = self.calculate_y(size, well_y + 10, gy2 - 25)
            self.canvas.create_line(
                ladder_x - band_w / 2,
                y,
                ladder_x + band_w / 2,
                y,
                fill="#ffd700",
                width=3,
            )
            self.canvas.create_text(
                ladder_x - band_w / 2 - 3,
                y,
                text=str(size),
                anchor="e",
                fill="#6c7a89",
                font=("Arial", 7),
            )

        for i, sample in enumerate(self.samples):
            lane_x = gx1 + (i + 2) * lane_step
            y = self.calculate_y(sample["length"], well_y + 10, gy2 - 25)

            self.canvas.create_line(
                lane_x - band_w / 2,
                y,
                lane_x + band_w / 2,
                y,
                fill="#ff9900",
                width=4,
            )
            self.canvas.create_line(
                lane_x - band_w / 2,
                y,
                lane_x + band_w / 2,
                y,
                fill="#ffffff",
                width=1,
            )


if __name__ == "__main__":
    root = tk.Tk()
    app = GelElectrophoresisApp(root)
    root.mainloop()