# flexible_lock_checker.py
import os

complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
def get_complement(seq):
    return "".join([complement.get(b, b) for b in seq])

with open('blocker_key.txt', 'r') as f:
    lock_key = f.readlines()[1].strip()

print(f"{'Файл':<25} | {'Макс. совпадение (%)'}")
print("-" * 50)

for file in [f for f in os.listdir('.') if f.endswith('.fasta')]:
    with open(file, 'r') as f:
        seq = "".join([l.strip().upper() for l in f if not l.startswith('>')])
    
    # Сканируем весь файл в поисках "прилипания" ключа (окнами по 56)
    best_match = 0
    for i in range(len(seq) - 56):
        chunk = seq[i:i+56]
        matches = sum(1 for a, b in zip(lock_key, get_complement(chunk)) if a == b)
        best_match = max(best_match, (matches / 56) * 100)
        
    print(f"{file:<25} | {best_match:.2f}%")