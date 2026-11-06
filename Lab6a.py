import random

def genereaza_adn(lungime):
    """Genereaza o secventa aleatoare de ADN de lungimea specificata"""
    return "".join(random.choice("ACGT") for _ in range(lungime))

def simuleaza_transpozonul(secventa_gazda, transpozon, nr_inserari=3):
    """
    Insereaza transpozonul in secventa gazda in mod aleator, 
    simuland si o scurta repetare directa de 3 baze.
    """
    secventa_modificata = secventa_gazda
    pozitii_reale = []
    
    for i in range(nr_inserari):
        poz = random.randint(10, len(secventa_modificata) - 10)
        
        dr = secventa_modificata[poz:poz+3]
        
        element_inserat = dr + transpozon + dr
        
        secventa_modificata = secventa_modificata[:poz] + element_inserat + secventa_modificata[poz:]
        pozitii_reale.append((poz, poz + len(element_inserat)))
        
    return secventa_modificata, pozitii_reale

def detecteaza_transpozoni(secventa_finala, transpozon):
    """
    Detecteaza pozitiile la care se afla transpozonul în secventa ADN.
    """
    pozitii_detectate = []
    start = 0
    while True:
        idx = secventa_finala.find(transpozon, start)
        if idx == -1:
            break
        start_complet = idx - 3 if idx >= 3 else 0
        end_complet = idx + len(transpozon) + 3
        pozitii_detectate.append((start_complet, end_complet, secventa_finala[start_complet:end_complet]))
        start = idx + 1
        
    return pozitii_detectate

lungime_gazda = 300
gazda = genereaza_adn(lungime_gazda)
print(f"Secventa gazda initiala (lungime {len(gazda)}):")
print(gazda)
print("-" * 60)

transpozon = "ATCGTACG"

secventa_noua, pozitii_reale = simuleaza_transpozonul(gazda, transpozon, nr_inserari=3)
print(f"Secvenaa ADN finala creata dupa inserari (lungime {len(secventa_noua)}):")
print(secventa_noua)
print("-" * 60)

rezultate_detectie = detecteaza_transpozoni(secventa_noua, transpozon)
print("Pozitiile si elementele detectate:")
for i, (start, end, segment) in enumerate(rezultate_detectie, 1):
    print(f"Elementul {i}: Interval [ {start} : {end} ], Secventa segment: {segment}")