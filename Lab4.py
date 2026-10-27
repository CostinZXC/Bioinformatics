import gzip
import os
import random
import threading
import tkinter as tk
from collections import defaultdict
from tkinter import messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import urllib.request


class DNAAssembler:

    def __init__(
        self,
        filepath="GCF_002600855.1_ASM260085v1_genomic.fna.gz",
        accession="NC_001422.1",
    ):
        self.filepath = filepath
        self.accession = accession
        self.genome = ""
        self.reads = []
        self.contigs = []

    def load_genome(self):
        """Încarcă genomul din arhiva locală sau îl preia de pe NCBI."""
        valid_bases = {"A", "C", "G", "T"}
        seq_parts = []

        if os.path.exists(self.filepath):
            open_func = gzip.open if self.filepath.endswith(".gz") else open
            mode = "rt" if self.filepath.endswith(".gz") else "r"
            with open_func(
                self.filepath, mode, encoding="utf-8", errors="ignore"
            ) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith(">"):
                        continue
                    for char in line.upper():
                        if char == "U":
                            char = "T"
                        if char in valid_bases:
                            seq_parts.append(char)
            self.genome = "".join(seq_parts)
            return len(self.genome), f"Fișier local: {self.filepath}"

        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={self.accession}&rettype=fasta&retmode=text"
        req = urllib.request.Request(
            url, headers={"User-Agent": "NCBI-Assembler-App"}
        )
        with urllib.request.urlopen(req) as resp:
            for line in resp.read().decode("utf-8").splitlines():
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                for char in line.upper():
                    if char in valid_bases:
                        seq_parts.append(char)
        self.genome = "".join(seq_parts)
        return len(self.genome), f"Descărcat de pe NCBI ({self.accession})"

    def sample_reads(self, num_reads=2000, min_len=100, max_len=120):
        """Extrage 2000 de secvențe aleatorii (reads) de lungime 100-120 nucleotide."""
        if not self.genome:
            return []
        g_len = len(self.genome)
        self.reads = []
        for _ in range(num_reads):
            r_len = random.randint(min_len, max_len)
            start = random.randint(0, g_len - r_len)
            self.reads.append(self.genome[start : start + r_len])
        return self.reads

    def _overlap(s1, s2, min_ov=20):
        start = 0
        while True:
            pos = s1.find(s2[:min_ov], start)
            if pos == -1:
                return 0
            if s2.startswith(s1[pos:]):
                return len(s1) - pos
            start = pos + 1

    def run_assembly(self, sample_size=300, min_overlap=25, progress_callback=None):
        """Asamblare Greedy pe un subset de reads pentru performanță."""
        pool = list(set(self.reads[:sample_size]))
        total_initial = len(pool)
        step = 0

        while len(pool) > 1:
            best_ov = 0
            best_pair = (None, None)

            for i in range(len(pool)):
                for j in range(len(pool)):
                    if i != j:
                        ov = self._overlap(pool[i], pool[j], min_overlap)
                        if ov > best_ov:
                            best_ov = ov
                            best_pair = (i, j)

            if best_ov < min_overlap:
                break

            i, j = best_pair
            merged = pool[i] + pool[j][best_ov:]
            pool = [read for idx, read in enumerate(pool) if idx not in (i, j)]
            pool.append(merged)
            step += 1

            if progress_callback and step % 15 == 0:
                progress_callback(step, total_initial, len(pool))

        self.contigs = sorted(pool, key=len, reverse=True)
        return self.contigs


class GenomeAssemblerApp(tk.Tk):

    def __init__(self, assembler):
        super().__init__()
        self.assembler = assembler

        self.title("DNA Sequence Reconstructor & Algorithmic Analysis")
        self.geometry("900x720")
        self.minsize(800, 600)

        self.setup_ui()

    def setup_ui(self):
        header_frame = ttk.LabelFrame(self, text="Informații Genom Sursă")
        header_frame.pack(fill="x", padx=15, pady=8)

        self.lbl_genome_info = ttk.Label(
            header_frame,
            text="Genom: Se încarcă...",
            font=("Arial", 10, "bold"),
            foreground="#0055A5",
        )
        self.lbl_genome_info.pack(anchor="w", padx=10, pady=5)

        ctrl_frame = ttk.Frame(self)
        ctrl_frame.pack(fill="x", padx=15, pady=5)

        self.btn_sample = ttk.Button(
            ctrl_frame,
            text="1. Generează 2000 Fragmente (100-120 pb)",
            command=self.on_sample_reads,
        )
        self.btn_sample.grid(row=0, column=0, padx=5, pady=5)

        self.btn_assemble = ttk.Button(
            ctrl_frame,
            text="2. Reconstruiește Secvența (Greedy)",
            command=self.start_assembly_thread,
            state="disabled",
        )
        self.btn_assemble.grid(row=0, column=1, padx=5, pady=5)

        self.btn_critique = ttk.Button(
            ctrl_frame,
            text="3. Probleme Algoritmice (Analiză)",
            command=self.open_analysis_popup,
        )
        self.btn_critique.grid(row=0, column=2, padx=5, pady=5)

        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(
            self, variable=self.progress_var, maximum=100
        )
        self.progress_bar.pack(fill="x", padx=15, pady=5)

        self.lbl_status = ttk.Label(
            self,
            text="Status: Așteaptă inițializarea datelor...",
            font=("Arial", 9, "italic"),
        )
        self.lbl_status.pack(anchor="w", padx=15, pady=(0, 5))

        table_frame = ttk.LabelFrame(
            self, text="Fragmente Extrase (2000 Reads) / Contiguri Asamblate"
        )
        table_frame.pack(fill="both", expand=True, padx=15, pady=5)

        columns = ("id", "type", "length", "sequence")
        self.tree = ttk.Treeview(
            table_frame, columns=columns, show="headings", height=12
        )
        self.tree.heading("id", text="Index")
        self.tree.heading("type", text="Tip Fragment")
        self.tree.heading("length", text="Lungime (pb)")
        self.tree.heading("sequence", text="Secvență (Primele 60 pb)")

        self.tree.column("id", width=70, anchor="center")
        self.tree.column("type", width=120, anchor="center")
        self.tree.column("length", width=100, anchor="e")
        self.tree.column("sequence", width=550, anchor="w")

        scrollbar = ttk.Scrollbar(
            table_frame, orient="vertical", command=self.tree.yview
        )
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        threading.Thread(target=self.init_data, daemon=True).start()

    def init_data(self):
        try:
            length, source = self.assembler.load_genome()
            self.lbl_genome_info.config(
                text=f"Sursă: {source} | Dimensiune Genom: {length:,} nucleotide"
            )
            self.lbl_status.config(
                text="Genom încărcat cu succes. Poți genera cele 2000 de fragmente."
            )
        except Exception as e:
            messagebox.showerror(
                "Eroare", f"Nu s-a putut încărca secvența: {str(e)}"
            )

    def on_sample_reads(self):
        reads = self.assembler.sample_reads(num_reads=2000)
        self.tree.delete(*self.tree.get_children())

        for idx, r in enumerate(reads[:300]):
            self.tree.insert(
                "",
                "end",
                values=(
                    f"Read {idx+1}",
                    "Fragment aleator",
                    f"{len(r)} pb",
                    r[:60] + "...",
                ),
            )

        self.lbl_status.config(
            text=f"Generat cu succes 2000 fragmente (afișate primele 300). Gata de asamblare."
        )
        self.btn_assemble.config(state="normal")

    def start_assembly_thread(self):
        self.btn_assemble.config(state="disabled")
        self.btn_sample.config(state="disabled")
        self.lbl_status.config(
            text="Asamblare în curs (algoritm Greedy Overlap)..."
        )
        self.progress_var.set(10)

        threading.Thread(target=self.run_assembly_task, daemon=True).start()

    def run_assembly_task(self):
        def update_prog(step, initial, current_pool):
            pct = min(95.0, (step / (initial * 0.75)) * 100)
            self.progress_var.set(pct)
            self.lbl_status.config(
                text=f"Îmbinare suprapuneri: Pasul {step} | Secvențe rămase: {current_pool}"
            )

        contigs = self.assembler.run_assembly(
            sample_size=250, min_overlap=25, progress_callback=update_prog
        )

        self.progress_var.set(100)
        self.after(0, self.display_assembly_results, contigs)

    def display_assembly_results(self, contigs):
        self.tree.delete(*self.tree.get_children())
        for idx, c in enumerate(contigs):
            self.tree.insert(
                "",
                "end",
                values=(
                    f"Contig {idx+1}",
                    "Reconstruit",
                    f"{len(c)} pb",
                    c[:60] + "...",
                ),
            )

        largest = len(contigs[0]) if contigs else 0
        self.lbl_status.config(
            text=f"Asamblare finalizată! Obținute {len(contigs)} contiguri. Cel mai mare contig: {largest} pb."
        )
        self.btn_sample.config(state="normal")
        self.btn_assemble.config(state="normal")

    def open_analysis_popup(self):
        """Deschide fereastra cu analiza critică a abordării algoritmice."""
        popup = tk.Toplevel(self)
        popup.title("Analiza Algoritmică: De ce reconstrucția perfectă eșuează")
        popup.geometry("680x520")

        lbl = ttk.Label(
            popup,
            text="Probleme fundamentale identificate conform cerinței:",
            font=("Arial", 11, "bold"),
        )
        lbl.pack(padx=10, pady=(10, 5), anchor="w")

        text_area = ScrolledText(popup, wrap="word", font=("Arial", 10))
        analysis_text = (
            "CELE MAI MARI PROBLEME ALE ABORDĂRII ALGORITMICE (GREEDY / OVERLAP):\n\n"
            "1. Repetițiile Genomice (Repeats):\n"
            "   - Genomul conține frecvent secvențe identice în locații complet diferite.\n"
            "   - Algoritmul nu are cum să diferențieze copia A de copia B folosind fragmente scurte de 100-120 pb.\n"
            "   - Rezultat: Regiunile dintre repetiții sunt omise, iar genomul asamblat este mult mai scurt (colapsare).\n\n"
            "2. Suprapuneri False și Joncțiuni Chimerice (Ambiguous Overlaps):\n"
            "   - Fragmente din zone opuse ale virusului se pot suprapune accidental pe o lungime de 25-30 baze.\n"
            "   - O decizie Greedy greșită unește două bucăți neînrudite, generând o secvență hibridă falsă (chimeric join).\n\n"
            "3. Goluri de Acoperire (Coverage Gaps):\n"
            "   - Deși avem 2000 de fragmente, eșantionarea fiind strict pseudo-aleatorie, anumite regiuni rămân neacoperite (acoperire = 0).\n"
            "   - Din această cauză, genomul nu poate fi închis într-un singur șir continuu, ci rămâne fragmentat în 'contiguri'.\n\n"
            "4. Complexitatea Computațională (NP-Hard):\n"
            "   - Problema reconstrucției prin cel mai scurt superșir comun (Shortest Common Superstring) este NP-hard.\n"
            "   - Compararea fiecărui fragment cu toate celelalte (O(N^2) sau O(N^3)) devine imposibil de lentă la scară mare fără euristici grafice (de Bruijn Graphs).\n\n"
            "5. Pierderea Orientării și Circularității:\n"
            "   - Multe genomuri virale sunt circulare. Reconstrucția liniară rupe continuitatea la capătul circular dacă nu există un pas dedicat de circularizare."
        )
        text_area.insert("1.0", analysis_text)
        text_area.pack(fill="both", expand=True, padx=10, pady=5)

        btn_close = ttk.Button(popup, text="Închide", command=popup.destroy)
        btn_close.pack(pady=5)


if __name__ == "__main__":
    assembler = DNAAssembler()
    app = GenomeAssemblerApp(assembler)
    app.mainloop()