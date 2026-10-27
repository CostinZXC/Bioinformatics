import gzip

import json

import os

import urllib.request

import tkinter as tk

from tkinter import ttk, messagebox, scrolledtext

from collections import defaultdict



CODON_TABLE = {

    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',

    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',

    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',

    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',

    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',

    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',

    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',

    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',

    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',

    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',

    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',

    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',

    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',

    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',

    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',

    'TGC': 'C', 'TGT': 'C', 'TGA': 'Stop', 'TGG': 'W',

}



AMINO_ACID_NAMES = {

    'A': 'Alanină', 'R': 'Arginină', 'N': 'Asparagină', 'D': 'Acid aspartic',

    'C': 'Cisteină', 'E': 'Acid glutamic', 'Q': 'Glutamină', 'G': 'Glicină',

    'H': 'Histidină', 'I': 'Izoleucină', 'L': 'Leucină', 'K': 'Lizină',

    'M': 'Metionină', 'F': 'Fenilalanină', 'P': 'Prolină', 'S': 'Serină',

    'T': 'Treonină', 'W': 'Triptofan', 'Y': 'Tirozină', 'V': 'Valină',

    'Stop': 'Codon Stop'

}





class CodonAnalyzer:

    def __init__(self, window_size=3, step=1):

        self.window_size = window_size

        self.step = step

        self.codon_counts = defaultdict(int)

        self.total_windows = 0

        self.sequence_preview = "" 



    def process_file(self, filepath, max_preview_bases=10000):

        buffer = []

        valid_bases = {'A', 'C', 'G', 'T'}

        preview_chars = []



        open_func = gzip.open if filepath.endswith('.gz') else open

        mode = 'rt' if filepath.endswith('.gz') else 'r'



        with open_func(filepath, mode, encoding='utf-8', errors='ignore') as file:

            for line in file:

                line = line.strip()

                if not line or line.startswith('>'):

                    continue



                for char in line.upper():

                    if char == 'U':

                        char = 'T'

                    if char in valid_bases:

                        if len(preview_chars) < max_preview_bases:

                            preview_chars.append(char)



                        buffer.append(char)

                        if len(buffer) == self.window_size:

                            codon = "".join(buffer)

                            self.codon_counts[codon] += 1

                            self.total_windows += 1

                            buffer.pop(0)



        self.sequence_preview = "".join(preview_chars)



    def get_top_codons(self, top_n=10):

        if self.total_windows == 0:

            return []



        sorted_codons = sorted(

            self.codon_counts.items(), key=lambda item: item[1], reverse=True

        )



        results = []

        for codon, count in sorted_codons[:top_n]:

            freq_percent = (count / self.total_windows) * 100

            aa_code = CODON_TABLE.get(codon, "Unknown")

            aa_name = AMINO_ACID_NAMES.get(aa_code, aa_code)

            results.append((codon, count, freq_percent, aa_name))



        return results





def build_ai_prompt(top_codons):

    summary_lines = []

    unique_aas = set()



    for codon, count, freq, aa_name in top_codons:

        summary_lines.append(f"- Codon {codon}: {freq:.3f}% ({count:,} apariții) -> {aa_name}")

        if aa_name != 'Codon Stop':

            unique_aas.add(aa_name)



    codons_text = "\n".join(summary_lines)

    amino_acids_list = ", ".join(sorted(unique_aas))



    prompt = (

        f"Am scanat un genom complet cu o fereastră glisantă de 3 baze (pas de 1) "

        f"și am obținut următorii top 10 codoni ca frecvență:\n\n{codons_text}\n\n"

        f"Aminoacizii unici asociați sunt: {amino_acids_list}.\n\n"

        f"Te rog să identifici alimentele sau grupele alimentare care sunt lipsite "

        f"sau extrem de sărace în acești aminoacizi specifici și să explici pe scurt de ce biochimic."

    )

    return prompt





def call_gemini_api_native(prompt, api_key):

    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key={api_key}"

    payload = {"contents": [{"parts": [{"text": prompt}]}]}

    req = urllib.request.Request(

        url,

        data=json.dumps(payload).encode("utf-8"),

        headers={"Content-Type": "application/json"}

    )

    with urllib.request.urlopen(req) as response:

        res_data = json.loads(response.read().decode("utf-8"))

        return res_data["candidates"][0]["content"]["parts"][0]["text"]





class GenomeApp(tk.Tk):

    def __init__(self, analyzer, top_codons, prompt_text):

        super().__init__()

        self.analyzer = analyzer

        self.top_codons = top_codons

        self.prompt_text = prompt_text

        self.seq = analyzer.sequence_preview



        self.title("Genome Codon Analyzer & AI Nutrition Prompt")

        self.geometry("880x700")



        self.setup_ui()



    def setup_ui(self):

        lbl_top = ttk.Label(self, text="Top 10 Frecvențe Codoni (Fereastră = 3, Pas = 1)", font=("Arial", 12, "bold"))

        lbl_top.pack(pady=(10, 5))



        columns = ("codon", "aa", "count", "freq")

        self.tree = ttk.Treeview(self, columns=columns, show="headings", height=8)

        self.tree.heading("codon", text="Codon")

        self.tree.heading("aa", text="Aminoacid")

        self.tree.heading("count", text="Număr Apariții")

        self.tree.heading("freq", text="Frecvență (%)")



        self.tree.column("codon", width=100, anchor="center")

        self.tree.column("aa", width=180, anchor="w")

        self.tree.column("count", width=140, anchor="e")

        self.tree.column("freq", width=120, anchor="e")



        for codon, count, freq, aa_name in self.top_codons:

            self.tree.insert("", "end", values=(codon, aa_name, f"{count:,}", f"{freq:.3f}%"))



        self.tree.pack(fill="x", padx=15, pady=5)



        frame_slider = ttk.LabelFrame(self, text="Scanare Secvență (Sliding Window - Vizualizator)")

        frame_slider.pack(fill="x", padx=15, pady=10)



        self.lbl_codon_info = ttk.Label(

            frame_slider, 

            text="Poziție: 0 | Codon Curent: --- | Traducere: ---", 

            font=("Arial", 10, "bold"),

            foreground="#0055A5"

        )

        self.lbl_codon_info.pack(anchor="w", padx=10, pady=(5, 2))



        max_idx = max(0, len(self.seq) - 3)

        self.slider = ttk.Scale(frame_slider, from_=0, to=max_idx, orient="horizontal", command=self.update_highlight)

        self.slider.pack(fill="x", padx=10, pady=5)



        self.txt_seq = tk.Text(frame_slider, height=4, font=("Courier", 11), wrap="char")

        self.txt_seq.pack(fill="x", padx=10, pady=(0, 10))

        self.txt_seq.tag_configure("highlight", background="yellow", foreground="black", font=("Courier", 11, "bold"))



        self.update_highlight(0)



        btn_frame = ttk.Frame(self)

        btn_frame.pack(pady=10)



        btn_popup_prompt = ttk.Button(btn_frame, text="Afișează Popup Prompt AI", command=self.open_prompt_popup)

        btn_popup_prompt.grid(row=0, column=0, padx=8)



        btn_call_ai = ttk.Button(btn_frame, text="Interoghează AI (Popup Răspuns)", command=self.open_ai_popup)

        btn_call_ai.grid(row=0, column=1, padx=8)



    def update_highlight(self, value):

        if not self.seq:

            return



        idx = int(float(value))

        current_codon = self.seq[idx:idx + 3]

        aa_code = CODON_TABLE.get(current_codon, "Unknown")

        aa_name = AMINO_ACID_NAMES.get(aa_code, aa_code)



        self.lbl_codon_info.config(

            text=f"Poziție baze: {idx} - {idx+2} | Codon: {current_codon} | Aminoacid: {aa_name}"

        )



        start = max(0, idx - 40)

        end = min(len(self.seq), idx + 80)

        displayed_text = self.seq[start:end]



        self.txt_seq.delete("1.0", tk.END)

        self.txt_seq.insert("1.0", displayed_text)



        highlight_start = f"1.{idx - start}"

        highlight_end = f"1.{idx - start + 3}"

        self.txt_seq.tag_add("highlight", highlight_start, highlight_end)



    def open_prompt_popup(self):

        """Deschide un Popup Toplevel cu textul promptului pregătit."""

        popup = tk.Toplevel(self)

        popup.title("Prompt Generat pentru AI")

        popup.geometry("620x450")



        lbl = ttk.Label(popup, text="Prompt gata asamblat pe baza frecvențelor genomice:", font=("Arial", 10, "bold"))

        lbl.pack(padx=10, pady=(10, 5), anchor="w")



        text_area = scrolledtext.ScrolledText(popup, wrap="word", font=("Arial", 10))

        text_area.insert("1.0", self.prompt_text)

        text_area.pack(fill="both", expand=True, padx=10, pady=5)



        btn_close = ttk.Button(popup, text="Închide", command=popup.destroy)

        btn_close.pack(pady=5)



    def open_ai_popup(self):

        """Apelează API-ul sau oferă direct răspunsul biochimic într-un popup."""

        api_key = os.environ.get("GEMINI_API_KEY")



        popup = tk.Toplevel(self)

        popup.title("Rezultat / Recomandare Nutrițională AI")

        popup.geometry("650x500")



        lbl = ttk.Label(popup, text="Analiză Nutrițională (Alimente lipsite de acești aminoacizi):", font=("Arial", 10, "bold"))

        lbl.pack(padx=10, pady=(10, 5), anchor="w")



        text_area = scrolledtext.ScrolledText(popup, wrap="word", font=("Arial", 10))

        text_area.pack(fill="both", expand=True, padx=10, pady=5)



        if api_key:

            try:

                text_area.insert(tk.END, "Se contactează Gemini API...\n\n")

                self.update_idletasks()

                response = call_gemini_api_native(self.prompt_text, api_key)

                text_area.delete("1.0", tk.END)

                text_area.insert(tk.END, response)

            except Exception as e:

                text_area.insert(tk.END, f"Eroare la apelul API: {str(e)}\n\nPromptul a rămas disponibil în prima fereastră.")

        else:

            default_ai_answer = (

                "[Cheia GEMINI_API_KEY nu este setată. Răspuns biochimic de referință:]\n\n"

                "1. Lipide și Grăsimi Pure (Zero Aminoacizi):\n"

                "   - Ulei de măsline, floarea-soarelui, unt clarifiat (ghee), untură pură.\n"

                "   - Motiv: Sunt compuse exclusiv din trigliceride și acizi grași, fără fracțiuni proteice.\n\n"

                "2. Zaharuri Rafinate și Amidon Pur:\n"

                "   - Sirop de glucoză, zahăr alb rafinat, amidon de porumb, tapioca.\n"

                "   - Motiv: Conțin polizaharide sau zaharuri simple și sunt lipsite de legături peptidice.\n\n"

                "3. Diete Medicale Speciale (Formule Hipoproteice):\n"

                "   - Preparate medicale concepute pentru afecțiuni metabolice ereditare (ex. formule fără fenilalanină/tirozină).\n"

            )

            text_area.insert(tk.END, default_ai_answer)



        btn_close = ttk.Button(popup, text="Închide", command=popup.destroy)

        btn_close.pack(pady=5)





if __name__ == "__main__":

    analyzer = CodonAnalyzer(window_size=3, step=1)

    genome_file = "GCF_002600855.1_ASM260085v1_genomic.fna.gz"



    if not os.path.exists(genome_file):

        genome_file = "test_seq.fna"

        with open(genome_file, "w") as f:

            f.write(">test\nATGGCGTTTCCCAAAATGATGATGAAATTTAAA")



    analyzer.process_file(genome_file)

    top_10 = analyzer.get_top_codons(10)

    prompt = build_ai_prompt(top_10)


    app = GenomeApp(analyzer, top_10, prompt)

    app.mainloop() 