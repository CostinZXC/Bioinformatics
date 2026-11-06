import tkinter as tk
from tkinter import ttk, messagebox
import random

ELEMENTE_DISPONIBILE = {
    "IS101": "TAGCTCGATTACGCTAGCTA",
    "IS1":   "GGTAATGACCCGGTACCGG",
    "IS5":   "CTAGTTACCGGATCGATAG"
}

def genereaza_adn(lungime):
    return "".join(random.choice("ACGT") for _ in range(lungime))

def insereaza_transpozoni(gazda, transpozon_seq, nr_inserari=3):
    secventa = gazda
    pozitii_reale = []
    
    for _ in range(nr_inserari):
        poz = random.randint(10, len(secventa) - 10)
        dr = secventa[poz:poz+3]
        bloc = dr + transpozon_seq + dr
        
        secventa = secventa[:poz] + bloc + secventa[poz:]
        pozitii_reale.append((poz, poz + len(bloc)))
        
    return secventa, sorted(pozitii_reale)

class TransposonDetectorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Detector Elemente Transpozabile (IS Finder)")
        self.root.geometry("780x680")
        self.root.configure(bg="#f4f6f8")

        self.secventa_generata = ""
        self.creeaza_interfata()

    def creeaza_interfata(self):
        control_frame = ttk.LabelFrame(self.root, text=" Parametri Simulare ", padding=12)
        control_frame.pack(fill="x", padx=15, pady=10)

        ttk.Label(control_frame, text="Element Transpozabil:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.combo_element = ttk.Combobox(control_frame, values=list(ELEMENTE_DISPONIBILE.keys()), state="readonly", width=12)
        self.combo_element.current(0)
        self.combo_element.grid(row=0, column=1, sticky="w", padx=5, pady=5)

        ttk.Label(control_frame, text="Nr. Inserari (3-4):").grid(row=0, column=2, sticky="w", padx=15, pady=5)
        self.spin_inserari = ttk.Spinbox(control_frame, from_=3, to=4, width=5)
        self.spin_inserari.set(3)
        self.spin_inserari.grid(row=0, column=3, sticky="w", padx=5, pady=5)

        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=1, column=0, columnspan=4, pady=10)

        self.btn_genereaza = ttk.Button(btn_frame, text="1. Genereaza Secventa ADN", command=self.action_genereaza)
        self.btn_genereaza.pack(side="left", padx=8)

        self.btn_detecteaza = ttk.Button(btn_frame, text="2. Detecteaza Pozitiile", command=self.action_detecteaza, state="disabled")
        self.btn_detecteaza.pack(side="left", padx=8)

        seq_frame = ttk.LabelFrame(self.root, text=" Secventa ADN Gazda ", padding=10)
        seq_frame.pack(fill="both", expand=True, padx=15, pady=5)

        self.text_adn = tk.Text(seq_frame, wrap="char", height=10, font=("Consolas", 11), bg="#ffffff")
        self.text_adn.pack(side="left", fill="both", expand=True)
        
        scroll_seq = ttk.Scrollbar(seq_frame, command=self.text_adn.yview)
        scroll_seq.pack(side="right", fill="y")
        self.text_adn.config(yscrollcommand=scroll_seq.set)

        self.text_adn.tag_config("highlight", background="#ffe066", foreground="#d90429", font=("Consolas", 11, "bold"))

        res_frame = ttk.LabelFrame(self.root, text=" Raport Detectie Pozitii ", padding=10)
        res_frame.pack(fill="both", expand=True, padx=15, pady=10)

        self.tree = ttk.Treeview(res_frame, columns=("Index", "Start", "End", "Lungime", "Fragment"), show="headings", height=5)
        self.tree.heading("Index", text="#")
        self.tree.heading("Start", text="Index Start")
        self.tree.heading("End", text="Index Sfarsit")
        self.tree.heading("Lungime", text="Lungime (bp)")
        self.tree.heading("Fragment", text="Secventa Detectata")

        self.tree.column("Index", width=40, anchor="center")
        self.tree.column("Start", width=90, anchor="center")
        self.tree.column("End", width=90, anchor="center")
        self.tree.column("Lungime", width=90, anchor="center")
        self.tree.column("Fragment", width=380, anchor="w")

        self.tree.pack(side="left", fill="both", expand=True)

        scroll_tree = ttk.Scrollbar(res_frame, command=self.tree.yview)
        scroll_tree.pack(side="right", fill="y")
        self.tree.config(yscrollcommand=scroll_tree.set)

    def action_genereaza(self):
        lungime_gazda = random.randint(200, 400)
        gazda = genereaza_adn(lungime_gazda)
        
        element_selectat = self.combo_element.get()
        seq_element = ELEMENTE_DISPONIBILE[element_selectat]
        nr_ins = int(self.spin_inserari.get())

        self.secventa_generata, _ = insereaza_transpozoni(gazda, seq_element, nr_inserari=nr_ins)

        self.text_adn.delete("1.0", tk.END)
        self.text_adn.insert("1.0", self.secventa_generata)
        for item in self.tree.get_children():
            self.tree.delete(item)

        self.btn_detecteaza.config(state="normal")

    def action_detecteaza(self):
        if not self.secventa_generata:
            messagebox.showwarning("Atentie", "Generati mai întai secventa ADN!")
            return

        element_selectat = self.combo_element.get()
        seq_element = ELEMENTE_DISPONIBILE[element_selectat]

        for item in self.tree.get_children():
            self.tree.delete(item)
        self.text_adn.tag_remove("highlight", "1.0", tk.END)

        start_idx = 0
        gasite = 0

        while True:
            pos = self.secventa_generata.find(seq_element, start_idx)
            if pos == -1:
                break
            
            s_pos = max(0, pos - 3)
            e_pos = min(len(self.secventa_generata), pos + len(seq_element) + 3)
            segment = self.secventa_generata[s_pos:e_pos]
            gasite += 1

            self.tree.insert("", "end", values=(gasite, s_pos, e_pos, len(segment), segment))

            tk_start = f"1.0 + {s_pos} chars"
            tk_end = f"1.0 + {e_pos} chars"
            self.text_adn.tag_add("highlight", tk_start, tk_end)

            start_idx = pos + len(seq_element)

        if gasite == 0:
            messagebox.showinfo("Rezultat", "Nu a fost gasit niciun element transpozabil.")

if __name__ == "__main__":
    root = tk.Tk()
    app = TransposonDetectorApp(root)
    root.mainloop()