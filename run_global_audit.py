import os
import numpy as np
from analysis_logic import analyze_structural_divergence

def get_total_tension(file_path):
    with open(file_path, 'r') as f:
        seq = "".join([l.strip().upper() for l in f if not l.startswith('>')])
    # Если файл пустой или слишком маленький
    if len(seq) < 100: return 0
    div = analyze_structural_divergence(seq)
    return np.sum(div)

files = [f for f in os.listdir('.') if f.endswith('.fasta')]
etalon_tension = get_total_tension('NC_004102.fasta')

print(f"{'Файл':<25} | {'Напряжение':<15} | {'Сходство с NC_004102':<20}")
print("-" * 65)
# Проверка: есть ли хоть один файл для анализа
if not files:
    print("Ошибка: Не найдено файлов .fasta. Положи их в папку со скриптом.")
    return
for file in files:
    tension = get_total_tension(file)
    # Относительное сходство
    similarity = (tension / etalon_tension) * 100 if etalon_tension > 0 else 0
    print(f"{file:<25} | {tension:<15.2f} | {similarity:<20.2f}%")
