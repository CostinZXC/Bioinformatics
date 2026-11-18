import re
import tkinter as tk
from tkinter import ttk, messagebox


def clean_dna(seq: str) -> str:
    """Curata secventa si pastreaza doar bazele A, C, G, T."""
    return re.sub(r'[^ACGT]', '', seq.upper().strip())


def needleman_wunsch(s1: str, s2: str, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2):
    """
    Implementarea algoritmului Needleman-Wunsch pentru aliniere globala.
    Returneaza matricea de scoruri, alinierea celor doua secvente, scorul final si coordonatele traseului.
    """
    n = len(s1)
    m = len(s2)

    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(1, m + 1):
        dp[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag_score = dp[i - 1][j - 1] + (match_score if s1[i - 1] == s2[j - 1] else mismatch_penalty)
            up_score = dp[i - 1][j] + gap_penalty
            left_score = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(diag_score, up_score, left_score)

    align1 = []
    align2 = []
    match_line = []
    path = {(n, m)}

    i, j = n, m
    while i > 0 or j > 0:
        current = dp[i][j]
        if i > 0 and j > 0 and current == dp[i - 1][j - 1] + (match_score if s1[i - 1] == s2[j - 1] else mismatch_penalty):
            align1.append(s1[i - 1])
            align2.append(s2[j - 1])
            match_line.append("|" if s1[i - 1] == s2[j - 1] else "·")
            i -= 1
            j -= 1
        elif i > 0 and current == dp[i - 1][j] + gap_penalty:
            align1.append(s1[i - 1])
            align2.append("-")
            match_line.append(" ")
            i -= 1
        else:
            align1.append("-")
            align2.append(s2[j - 1])
            match_line.append(" ")
            j -= 1
        path.add((i, j))

    align1.reverse()
    align2.reverse()
    match_line.reverse()

    aligned_s1 = "".join(align1)
    aligned_s2 = "".join(align2)
    middle_str = "".join(match_line)

    matches_count = sum(1 for a, b in zip(aligned_s1, aligned_s2) if a == b and a != '-')
    identity_pct = (matches_count / len(aligned_s1)) * 100.0 if aligned_s1 else 0.0

    return {
        "dp": dp,
        "aligned_s1": aligned_s1,
        "aligned_s2": aligned_s2,
        "middle_str": middle_str,
        "final_score": dp[n][m],
        "matches": matches_count,
        "length": len(aligned_s1),
        "identity": identity_pct,
        "path": path
    }


class GlobalAlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Global DNA Sequence Alignment (Needleman-Wunsch)")
        self.root.geometry("950x700")

        self._setup_ui()
        self.run_alignment()

    def _setup_ui(self):
        top_frame = ttk.LabelFrame(self.root, text="Secvente ADN & Parametri de Scor", padding=10)
        top_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(top_frame, text="Secventa S1:").grid(row=0, column=0, sticky=tk.W, pady=3)
        self.s1_entry = ttk.Entry(top_frame, width=45, font=("Consolas", 11))
        self.s1_entry.insert(0, "ACCGTGAAGCCAATAC")
        self.s1_entry.grid(row=0, column=1, columnspan=5, sticky=tk.W, padx=5, pady=3)

        ttk.Label(top_frame, text="Secventa S2:").grid(row=1, column=0, sticky=tk.W, pady=3)
        self.s2_entry = ttk.Entry(top_frame, width=45, font=("Consolas", 11))
        self.s2_entry.insert(0, "AGCGTGCAGCCAATAC")
        self.s2_entry.grid(row=1, column=1, columnspan=5, sticky=tk.W, padx=5, pady=3)

        ttk.Label(top_frame, text="Match:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.match_entry = ttk.Entry(top_frame, width=6)
        self.match_entry.insert(0, "1")
        self.match_entry.grid(row=2, column=1, sticky=tk.W, padx=5)

        ttk.Label(top_frame, text="Mismatch:").grid(row=2, column=2, sticky=tk.W, pady=5)
        self.mismatch_entry = ttk.Entry(top_frame, width=6)
        self.mismatch_entry.insert(0, "-1")
        self.mismatch_entry.grid(row=2, column=3, sticky=tk.W, padx=5)

        ttk.Label(top_frame, text="Gap:").grid(row=2, column=4, sticky=tk.W, pady=5)
        self.gap_entry = ttk.Entry(top_frame, width=6)
        self.gap_entry.insert(0, "-2")
        self.gap_entry.grid(row=2, column=5, sticky=tk.W, padx=5)

        ttk.Button(top_frame, text="Aliniaza Secventele", command=self.run_alignment).grid(row=0, column=6, rowspan=3, padx=20, ipadx=10, ipady=10)

        res_frame = ttk.LabelFrame(self.root, text="Rezultat Aliniere Globala", padding=10)
        res_frame.pack(fill=tk.X, padx=10, pady=5)

        self.stats_label = ttk.Label(res_frame, text="", font=("Arial", 10, "bold"), foreground="#1f77b4")
        self.stats_label.pack(anchor=tk.W, pady=(0, 5))

        self.align_text = tk.Text(res_frame, height=4, font=("Consolas", 13, "bold"), bg="#f8f9fa", relief=tk.FLAT)
        self.align_text.pack(fill=tk.X)
        self.align_text.config(state=tk.DISABLED)

        matrix_frame = ttk.LabelFrame(self.root, text="Matricea de Scoruri & Traseul Optim (Traceback)", padding=10)
        matrix_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.canvas = tk.Canvas(matrix_frame, bg="white")
        v_scroll = ttk.Scrollbar(matrix_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        h_scroll = ttk.Scrollbar(matrix_frame, orient=tk.HORIZONTAL, command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=v_scroll.set, xscrollcommand=h_scroll.set)

        v_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        h_scroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    def run_alignment(self):
        s1 = clean_dna(self.s1_entry.get())
        s2 = clean_dna(self.s2_entry.get())

        if not s1 or not s2:
            messagebox.showerror("Eroare", "Ambele secvente trebuie sa contina nucleotide valide (A, C, G, T).")
            return

        try:
            match_s = int(self.match_entry.get())
            mismatch_s = int(self.mismatch_entry.get())
            gap_s = int(self.gap_entry.get())
        except ValueError:
            messagebox.showerror("Eroare", "Valorile pentru Match, Mismatch si Gap trebuie sa fie numere intregi.")
            return

        res = needleman_wunsch(s1, s2, match_s, mismatch_s, gap_s)

        self.stats_label.config(
            text=f"Scor Final: {res['final_score']}   |   "
                 f"Potriviri (Matches): {res['matches']}/{res['length']}   |   "
                 f"Identitate: {res['identity']:.2f}%"
        )

        self.align_text.config(state=tk.NORMAL)
        self.align_text.delete("1.0", tk.END)
        output = (
            f"S1: {res['aligned_s1']}\n"
            f"    {res['middle_str']}\n"
            f"S2: {res['aligned_s2']}"
        )
        self.align_text.insert(tk.END, output)
        self.align_text.config(state=tk.DISABLED)

        self.draw_matrix(s1, s2, res["dp"], res["path"])

    def draw_matrix(self, s1, s2, dp, path):
        self.canvas.delete("all")
        cell_w = 42
        cell_h = 32
        offset_x = 15
        offset_y = 15

        rows = len(s1) + 1
        cols = len(s2) + 1

        for j in range(len(s2)):
            x = offset_x + (j + 2) * cell_w
            y = offset_y
            self.canvas.create_rectangle(x, y, x + cell_w, y + cell_h, fill="#e9ecef", outline="#adb5bd")
            self.canvas.create_text(x + cell_w / 2, y + cell_h / 2, text=s2[j], font=("Consolas", 11, "bold"))

        for i in range(len(s1)):
            x = offset_x
            y = offset_y + (i + 2) * cell_h
            self.canvas.create_rectangle(x, y, x + cell_w, y + cell_h, fill="#e9ecef", outline="#adb5bd")
            self.canvas.create_text(x + cell_w / 2, y + cell_h / 2, text=s1[i], font=("Consolas", 11, "bold"))

        for i in range(rows):
            for j in range(cols):
                x = offset_x + (j + 1) * cell_w
                y = offset_y + (i + 1) * cell_h

                is_path = (i, j) in path
                bg_color = "#ffd166" if is_path else "#ffffff"
                text_font = ("Consolas", 10, "bold") if is_path else ("Consolas", 10)

                self.canvas.create_rectangle(x, y, x + cell_w, y + cell_h, fill=bg_color, outline="#dee2e6")
                self.canvas.create_text(x + cell_w / 2, y + cell_h / 2, text=str(dp[i][j]), font=text_font)

        total_w = offset_x * 2 + (cols + 1) * cell_w
        total_h = offset_y * 2 + (rows + 1) * cell_h
        self.canvas.configure(scrollregion=(0, 0, total_w, total_h))


if __name__ == "__main__":
    root = tk.Tk()
    app = GlobalAlignmentApp(root)
    root.mainloop()