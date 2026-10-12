import math

"""def calculate_tm(sequence):
    seq = sequence.upper().strip()

    count_a = seq.count("A")
    count_t = seq.count("T")
    count_c = seq.count("C")
    count_g = seq.count("G")

    tm = 4 * (count_g + count_c) + 2 * (count_a + count_t)
    return tm


dna = input("DNA sequence: ")
result = calculate_tm(dna)
print(f"Melting temperature (Tm) = {result}°C")"""


def calculate_tm_both(sequence, na_conc=0.05):
    seq = sequence.upper().strip()
    length = len(seq)

    if length == 0:
        return 0, 0

    count_a = seq.count("A")
    count_t = seq.count("T")
    count_c = seq.count("C")
    count_g = seq.count("G")

    tm_simple = 4 * (count_g + count_c) + 2 * (count_a + count_t)

    gc_fraction = (count_g + count_c) / length
    tm_advanced = (
        81.5 + 16.6 * math.log10(na_conc) + 41 * gc_fraction - (600 / length)
    )

    return tm_simple, tm_advanced

dna = "ACGCGTGCCA"

tm1, tm2 = calculate_tm_both(dna, na_conc=0.05)

print(f"Sequencee: {dna}")
print(f"1. Tm (simple formula): {tm1:.2f}°C")
print(f"2. Tm (second task): {tm2:.2f}°C")