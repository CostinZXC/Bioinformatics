def digest_dna(sequence, enzymes):
    """
    Simuleaza digestia enzimatica a unei secvente liniare de ADN.
    
    :param sequence: Secventa de ADN ca string 
    :param enzymes: Dictionar cu numele enzimei, tiparul de recunoastere si pozitia de taiere relativa
    :return: Lista fragmentelor rezultate cu lungimile lor
    """
    cut_positions = set()
    found_sites = []

    for name, info in enzymes.items():
        pattern = info["site"]
        cut_offset = info["cut_offset"]
        
        start = 0
        while True:
            idx = sequence.find(pattern, start)
            if idx == -1:
                break
            
            cut_pos = idx + cut_offset
            cut_positions.add(cut_pos)
            found_sites.append((cut_pos, name, idx + 1))
            start = idx + 1

    sorted_cuts = sorted(list(cut_positions))
    
    fragments = []
    prev = 0
    for cut in sorted_cuts:
        frag = sequence[prev:cut]
        fragments.append(frag)
        prev = cut
    fragments.append(sequence[prev:])

    return found_sites, fragments


dna_seq = (
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

enzymes = {
    "EcoRI":   {"site": "GAATTC", "cut_offset": 1},
    "BamHI":   {"site": "GGATCC", "cut_offset": 1},
    "HindIII": {"site": "AAGCTT", "cut_offset": 1}
}

sites, fragments = digest_dna(dna_seq, enzymes)

print(f"Lungime totala secventa: {len(dna_seq)} pb\n")
print("Situsuri de taiere identificate:")
for cut_pos, name, start_idx in sorted(sites):
    print(f"  - {name}: recunoaste la baza {start_idx}, taie dupa baza {cut_pos}")

print(f"\nNumar fragmente rezultate: {len(fragments)}")
for i, frag in enumerate(fragments, 1):
    print(f"  Fragment {i}: {len(frag)} pb")