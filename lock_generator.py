# lock_generator.py
def get_complement(sequence):
    complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return "".join([complement.get(base, base) for base in sequence])

with open('NC_004102.fasta', 'r') as f:
    seq = "".join([l.strip().upper() for l in f if not l.startswith('>')])
    node_seq = seq[9080:9135]

lock_key = get_complement(node_seq)

with open('blocker_key.txt', 'w') as f:
    f.write(f"LOCK_KEY_HCV_9080_9135:\n{lock_key}")

print("--- Ключ-блокатор сгенерирован и сохранен в blocker_key.txt ---")