import re
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def expand_sequence(seq: str) -> str:
    """Extinde notatii de tipul (A27) în secventa completa de nucleotide"""
    seq = seq.upper().strip()
    seq = re.sub(r'\(([ACGT])(\d+)\)', lambda m: m.group(1) * int(m.group(2)), seq)
    return re.sub(r'[^ACGT]', '', seq)


def compute_cg(window: str) -> float:
    """Calculeaza procentul C+G pentru o fereastra"""
    if not window:
        return 0.0
    cg_count = window.count('C') + window.count('G')
    return (cg_count / len(window)) * 100.0


def compute_ic(window: str) -> float:
    """
    Calculeaza Information Content (Ic)
    """
    n = len(window)
    if n <= 1:
        return 0.0

    total_sum = 0.0
    for u in range(1, n):
        matches = 0
        limit = n - u
        for i in range(limit):
            if window[i] == window[i + u]:
                matches += 1
        total_sum += (matches / limit) * 100.0

    ic_value = total_sum / (n - 1)
    return ic_value


def analyze_sequence(seq: str, window_size: int = 30):
    """Proceseaza întreaga secventa folosind fereastra glisanta"""
    clean_seq = expand_sequence(seq)
    if len(clean_seq) < window_size:
        raise ValueError(f"Secventa are {len(clean_seq)} baze, dar fereastra necesita minim {window_size}.")

    cg_values = []
    ic_values = []

    for i in range(len(clean_seq) - window_size + 1):
        window = clean_seq[i:i + window_size]
        cg_values.append(compute_cg(window))
        ic_values.append(compute_ic(window))

    mean_cg = sum(cg_values) / len(cg_values)
    mean_ic = sum(ic_values) / len(ic_values)

    return {
        "clean_seq": clean_seq,
        "cg_list": cg_values,
        "ic_list": ic_values,
        "center_x": mean_cg,
        "center_y": mean_ic
    }


class PromKappaApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Pattern Generator (PromKappa Algorithm)")
        self.root.geometry("1150x720")

        self.centers_history = []

        self._setup_ui()
        self.run_analysis()

    def _setup_ui(self):
        control_frame = ttk.LabelFrame(self.root, text="Date de intrare & Parametri", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(control_frame, text="Nume Promotor / Eticheta:").grid(row=0, column=0, sticky=tk.W)
        self.name_entry = ttk.Entry(control_frame, width=25)
        self.name_entry.insert(0, "Secventa Test (Pasul 1)")
        self.name_entry.grid(row=0, column=1, padx=5, pady=2, sticky=tk.W)

        ttk.Label(control_frame, text="Lungime fereastra (baze):").grid(row=0, column=2, padx=(15, 0), sticky=tk.W)
        self.win_entry = ttk.Entry(control_frame, width=8)
        self.win_entry.insert(0, "30")
        self.win_entry.grid(row=0, column=3, padx=5, pady=2, sticky=tk.W)

        ttk.Label(control_frame, text="Secventa ADN:").grid(row=1, column=0, sticky=tk.NW, pady=5)
        self.seq_text = scrolledtext.ScrolledText(control_frame, height=4, width=85, font=("Consolas", 10))
        default_seq = "CGGACTGATCTATCT(A27)CGTAGCATCTATCGATCTATCTCTAGCGATCTA"
        self.seq_text.insert(tk.END, default_seq)
        self.seq_text.grid(row=1, column=1, columnspan=4, padx=5, pady=5, sticky=tk.W)

        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=1, column=5, padx=10, sticky=tk.N)

        ttk.Button(btn_frame, text="Genereaza Pattern", command=self.run_analysis).pack(fill=tk.X, pady=2)
        ttk.Button(btn_frame, text="Încarca Promotor Exemplu", command=self.load_sample_promoter).pack(fill=tk.X, pady=2)
        ttk.Button(btn_frame, text="Reseteaza Centrele (Grafic 2)", command=self.clear_centers).pack(fill=tk.X, pady=2)

        res_frame = ttk.LabelFrame(self.root, text="Rezultate Numerice & Verificare", padding=10)
        res_frame.pack(fill=tk.X, padx=10, pady=5)

        self.lbl_len = ttk.Label(res_frame, text="Lungime secventa: -", font=("Arial", 10))
        self.lbl_len.pack(side=tk.LEFT, padx=15)

        self.lbl_cg = ttk.Label(res_frame, text="CG% mediu (Centru X): -", font=("Arial", 10, "bold"), foreground="#1f77b4")
        self.lbl_cg.pack(side=tk.LEFT, padx=15)

        self.lbl_ic = ttk.Label(res_frame, text="Ic mediu (Centru Y): -", font=("Arial", 10, "bold"), foreground="#d62728")
        self.lbl_ic.pack(side=tk.LEFT, padx=15)

        charts_frame = ttk.Frame(self.root)
        charts_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(11, 4.5), dpi=100)
        self.fig.tight_layout(pad=3.5)

        self.canvas = FigureCanvasTkAgg(self.fig, master=charts_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def run_analysis(self):
        seq = self.seq_text.get("1.0", tk.END).strip()
        name = self.name_entry.get().strip() or "Promotor"
        try:
            win_size = int(self.win_entry.get().strip())
            res = analyze_sequence(seq, win_size)
        except ValueError as e:
            messagebox.showerror("Eroare", str(e))
            return

        self.lbl_len.config(text=f"Lungime procesata: {len(res['clean_seq'])} bp ({len(res['cg_list'])} ferestre)")
        self.lbl_cg.config(text=f"CG% mediu (Centru X): {res['center_x']:.2f}%")
        self.lbl_ic.config(text=f"Ic mediu (Centru Y): {res['center_y']:.2f}%")

        center_entry = (name, res["center_x"], res["center_y"])
        if not self.centers_history or self.centers_history[-1] != center_entry:
            self.centers_history.append(center_entry)

        self.update_charts(res, name)

    def update_charts(self, res, name):
        self.ax1.clear()
        self.ax1.scatter(res["cg_list"], res["ic_list"], c="#2b5c8f", alpha=0.65, edgecolors="none", s=35, label="Ferestre glisante")
        self.ax1.scatter([res["center_x"]], [res["center_y"]], c="red", marker="X", s=140, zorder=5, label=f"Centru ({res['center_x']:.2f}, {res['center_y']:.2f})")

        self.ax1.set_title(f"Pattern ADN: {name}")
        self.ax1.set_xlabel("Continut C+G (%)")
        self.ax1.set_ylabel("Information Content / Ic (%)")
        self.ax1.set_xlim(0, 100)
        self.ax1.set_ylim(0, 100)
        self.ax1.grid(True, linestyle="--", alpha=0.5)
        self.ax1.legend(loc="upper right", fontsize=8)

        self.ax2.clear()
        for idx, (p_name, cx, cy) in enumerate(self.centers_history):
            self.ax2.scatter(cx, cy, s=90, zorder=4)
            self.ax2.annotate(f"{idx+1}. {p_name}", (cx, cy), textcoords="offset points", xytext=(5, 5), fontsize=8)

        self.ax2.set_title("Centrele de greutate (Comparatie Promotori)")
        self.ax2.set_xlabel("Media C+G (%)")
        self.ax2.set_ylabel("Media Ic (%)")
        self.ax2.set_xlim(0, 100)
        self.ax2.set_ylim(0, 100)
        self.ax2.grid(True, linestyle="--", alpha=0.5)

        self.canvas.draw()

    def load_sample_promoter(self):
        """Încarca o secventa reala de promotor"""
        sample_promoter = (
            "GCGCGGGGCCGGGCGCGGCGGGGCGGGGCGGGCGCGCGGGGCCGGGCGCGGCGGGGCGGG"
            "TATAAAAGGCGCGGGGCCGGGCGCGGCGGGGCGGGGCGGGCGCGCGGGGCCGGGCGCGGC"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTACGTACGTACGTACGTACGTACGTAC"
            "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCT"
        )
        self.name_entry.delete(0, tk.END)
        self.name_entry.insert(0, "Promotor Sintetic/Test 2")
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert(tk.END, sample_promoter)
        self.run_analysis()

    def clear_centers(self):
        self.centers_history.clear()
        self.ax2.clear()
        self.ax2.set_title("Centrele de greutate (Comparatie Promotori)")
        self.ax2.set_xlabel("Media C+G (%)")
        self.ax2.set_ylabel("Media Ic (%)")
        self.ax2.set_xlim(0, 100)
        self.ax2.set_ylim(0, 100)
        self.ax2.grid(True, linestyle="--", alpha=0.5)
        self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = PromKappaApp(root)
    root.mainloop()