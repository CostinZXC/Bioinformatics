import gzip
import itertools
import threading
import tkinter as tk
from collections import defaultdict
from tkinter import filedialog, messagebox, ttk

DEFAULT_NCBI_DNA = (
    "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTG"
    "AGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGG"
    "GCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGA"
    "TGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCT"
    "GGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGT"
    "GGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAA"
    "TTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGT"
    "ATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTA"
    "AACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC"
    "AATGATGTATTTAAATTATTTCTGAATATTTTACTAAAAAGGGAATGTGGGAGGTCAGTGCATTTAAAAC"
    "ATAAAGAAATGAATGTTTTATTTCCTTCCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    "TTTTTTTTTTTTTTTCTTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTT"
    "CTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCT"
    "TTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTT"
    "TTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTC"
    "TTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCC"
    "TTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTT"
    "TCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTC"
    "CCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTT"
    "TTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTT"
    "TCCCTTTCCTTTTTCTTTTTTCCCTTTCTTTTTTCTTTTTTTCCTTTTCTTTTTTCCCTTTCCTTTTTCT"
    "TTTTTCCCTTTCTTT"
)


class DNARepeatFinderApp:

    def __init__(self, root):
        self.root = root
        self.root.title("NCBI DNA Brute-Force Repetition Detector")
        self.root.geometry("1100x720")
        self.root.minsize(900, 600)
        self.root.configure(bg="#1e1e1e")

        self.current_dna_seq = DEFAULT_NCBI_DNA
        self.results_data = []

        self.setup_ui()
        self.update_dna_view()

    def setup_ui(self):
        top_bar = tk.Frame(self.root, bg="#2d2d2d", padx=12, pady=10)
        top_bar.pack(side=tk.TOP, fill=tk.X)

        self.btn_load = ttk.Button(
            top_bar,
            text="Încarcă din fișier (.fna.gz / .fna)",
            command=self.load_ncbi_file,
        )
        self.btn_load.pack(side=tk.LEFT, padx=5)

        self.btn_default = ttk.Button(
            top_bar,
            text="Resetează la Secvență Demo NCBI",
            command=self.use_default_sequence,
        )
        self.btn_default.pack(side=tk.LEFT, padx=5)

        self.btn_scan = ttk.Button(
            top_bar,
            text="Pornește Motorul Brute-Force",
            command=self.start_scan_thread,
        )
        self.btn_scan.pack(side=tk.LEFT, padx=10)

        self.lbl_status = tk.Label(
            top_bar,
            text=f"Secvență încărcată: {len(self.current_dna_seq)} nucleotide",
            fg="#61dafb",
            bg="#2d2d2d",
            font=("Arial", 10, "bold"),
        )
        self.lbl_status.pack(side=tk.LEFT, padx=10)

        self.progress = ttk.Progressbar(top_bar, mode="indeterminate", length=140)

        paned = tk.PanedWindow(
            self.root, orient=tk.VERTICAL, bg="#1e1e1e", sashwidth=5
        )
        paned.pack(fill=tk.BOTH, expand=True, padx=12, pady=10)

        seq_frame = tk.LabelFrame(
            paned,
            text="Secvență ADN NCBI (Selectează un rând din tabel pentru a evidenția repetiția)",
            fg="#9cdcfe",
            bg="#252526",
            font=("Arial", 10, "bold"),
            padx=8,
            pady=8,
        )
        paned.add(seq_frame, height=180)

        self.txt_dna = tk.Text(
            seq_frame,
            wrap=tk.CHAR,
            bg="#181818",
            fg="#d4d4d4",
            insertbackground="white",
            font=("Consolas", 10),
            padx=8,
            pady=8,
        )
        seq_scroll = ttk.Scrollbar(
            seq_frame, orient=tk.VERTICAL, command=self.txt_dna.yview
        )
        self.txt_dna.configure(yscrollcommand=seq_scroll.set)
        seq_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_dna.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.txt_dna.tag_config(
            "highlight", background="#ffd700", foreground="#000000"
        )

        table_frame = tk.LabelFrame(
            paned,
            text="Repetiții Detectate (Combinații brute-force de 3-6 litere, lungime totală 6-10 litere)",
            fg="#4ec9b0",
            bg="#252526",
            font=("Arial", 10, "bold"),
            padx=8,
            pady=8,
        )
        paned.add(table_frame, height=440)

        cols = (
            "pattern",
            "motif",
            "motif_len",
            "repeats",
            "total_len",
            "count",
            "positions",
        )
        self.tree = ttk.Treeview(
            table_frame, columns=cols, show="headings", selectmode="browse"
        )

        self.tree.heading("pattern", text="Tipar Repetat")
        self.tree.heading("motif", text="Motiv de Bază")
        self.tree.heading("motif_len", text="Lung. Motiv (3-6)")
        self.tree.heading("repeats", text="Multiplicator (n)")
        self.tree.heading("total_len", text="Lungime Totală (6-10)")
        self.tree.heading("count", text="Apariții")
        self.tree.heading("positions", text="Poziții (Indici de Start)")

        self.tree.column("pattern", width=120, anchor="center")
        self.tree.column("motif", width=90, anchor="center")
        self.tree.column("motif_len", width=120, anchor="center")
        self.tree.column("repeats", width=110, anchor="center")
        self.tree.column("total_len", width=140, anchor="center")
        self.tree.column("count", width=80, anchor="center")
        self.tree.column("positions", width=340, anchor="w")

        tree_scroll_y = ttk.Scrollbar(
            table_frame, orient=tk.VERTICAL, command=self.tree.yview
        )
        tree_scroll_x = ttk.Scrollbar(
            table_frame, orient=tk.HORIZONTAL, command=self.tree.xview
        )
        self.tree.configure(
            yscrollcommand=tree_scroll_y.set, xscrollcommand=tree_scroll_x.set
        )

        tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.tree.pack(fill=tk.BOTH, expand=True)

        self.tree.bind("<<TreeviewSelect>>", self.on_tree_select)

    def update_dna_view(self):
        self.txt_dna.delete("1.0", tk.END)
        self.txt_dna.insert("1.0", self.current_dna_seq)
        self.lbl_status.config(
            text=f"Secvență încărcată: {len(self.current_dna_seq)} nucleotide"
        )

    def use_default_sequence(self):
        self.current_dna_seq = DEFAULT_NCBI_DNA
        self.update_dna_view()
        self.clear_table()

    def clear_table(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.txt_dna.tag_remove("highlight", "1.0", tk.END)

    def load_ncbi_file(self):
        file_path = filedialog.askopenfilename(
            title="Alege fișierul descărcat de pe NCBI",
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

        open_func = (
            gzip.open
            if file_path.endswith(".gz")
            else open
        )
        chars = []

        try:
            with open_func(file_path, "rt", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith(">"):
                        continue
                    cleaned = "".join(
                        c for c in line.upper() if c in ("A", "C", "G", "T")
                    )
                    chars.append(cleaned)
                    if sum(len(x) for x in chars) >= 2000:
                        break

            extracted_seq = "".join(chars)[:2000]

            if len(extracted_seq) < 1000:
                messagebox.showerror(
                    "Eroare",
                    f"Secvența găsită are doar {len(extracted_seq)} nucleotide. Este nevoie de cel puțin 1000.",
                )
                return

            self.current_dna_seq = extracted_seq
            self.update_dna_view()
            self.clear_table()
            messagebox.showinfo(
                "Succes",
                f"Au fost extrase cu succes {len(extracted_seq)} nucleotide din fișier!",
            )

        except Exception as ex:
            messagebox.showerror(
                "Eroare la citire", f"Nu s-a putut procesa fișierul:\n{ex}"
            )

    def start_scan_thread(self):
        self.btn_scan.config(state=tk.DISABLED)
        self.btn_load.config(state=tk.DISABLED)
        self.progress.pack(side=tk.LEFT, padx=10)
        self.progress.start(10)
        self.lbl_status.config(
            text="Rulare motor brute-force...", fg="#e5c07b"
        )

        threading.Thread(target=self.run_brute_force, daemon=True).start()

    def run_brute_force(self):
        """Generează exhaustiv toate motivele de 3-6 litere și verifică

        repetițiile lor în secvență (lungime totală 6-10).
        """
        alphabet = ("A", "C", "G", "T")
        seq = self.current_dna_seq
        found_data = []

        for m_len in range(3, 7):
            for motif_tuple in itertools.product(alphabet, repeat=m_len):
                motif = "".join(motif_tuple)

                min_rep = max(2, (6 + m_len - 1) // m_len)
                max_rep = 10 // m_len

                if min_rep > max_rep:
                    continue

                for n in range(min_rep, max_rep + 1):
                    repeated_pattern = motif * n
                    pat_len = len(repeated_pattern)

                    if 6 <= pat_len <= 10:
                        indices = []
                        start_pos = 0
                        while True:
                            idx = seq.find(repeated_pattern, start_pos)
                            if idx == -1:
                                break
                            indices.append(idx)
                            start_pos = idx + 1

                        if indices:
                            found_data.append(
                                {
                                    "pattern": repeated_pattern,
                                    "motif": motif,
                                    "motif_len": m_len,
                                    "repeats": n,
                                    "total_len": pat_len,
                                    "count": len(indices),
                                    "positions": indices,
                                }
                            )

        found_data.sort(
            key=lambda x: (x["motif_len"], x["total_len"], x["pattern"])
        )

        self.root.after(0, self.on_scan_finished, found_data)

    def on_scan_finished(self, found_data):
        self.progress.stop()
        self.progress.pack_forget()
        self.btn_scan.config(state=tk.NORMAL)
        self.btn_load.config(state=tk.NORMAL)
        self.lbl_status.config(
            text=f"Scanare finalizată! Găsite: {len(found_data)} repetiții distincte.",
            fg="#98c379",
        )

        self.clear_table()
        self.results_data = found_data

        for item in found_data:
            pos_str = ", ".join(map(str, item["positions"]))
            self.tree.insert(
                "",
                tk.END,
                values=(
                    item["pattern"],
                    item["motif"],
                    item["motif_len"],
                    item["repeats"],
                    item["total_len"],
                    item["count"],
                    pos_str,
                ),
            )

    def on_tree_select(self, event):
        """Când utilizatorul selectează un rând din tabel, evidențiază toate

        aparițiile pe textul secvenței ADN.
        """
        selected_items = self.tree.selection()
        if not selected_items:
            return

        item = self.tree.item(selected_items[0])
        vals = item["values"]
        pattern = vals[0]

        self.txt_dna.tag_remove("highlight", "1.0", tk.END)

        seq = self.current_dna_seq
        first_pos = None
        start = 0
        while True:
            idx = seq.find(pattern, start)
            if idx == -1:
                break
            if first_pos is None:
                first_pos = idx

            start_index = f"1.0 + {idx} chars"
            end_index = f"1.0 + {idx + len(pattern)} chars"
            self.txt_dna.tag_add("highlight", start_index, end_index)
            start = idx + 1

        if first_pos is not None:
            self.txt_dna.see(f"1.0 + {first_pos} chars")


if __name__ == "__main__":
    root = tk.Tk()
    app = DNARepeatFinderApp(root)
    root.mainloop()