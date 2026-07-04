import os
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from analysis_logic import analyze_structural_divergence

def run_batch_analysis(directory='.'):
    results = []
    # Ищем все файлы .fasta, кроме эталонного
    files = [f for f in os.listdir(directory) if f.endswith('.fasta') and 'NC_009828' not in f]
    
    print(f"Найдено геномов для анализа: {len(files)}")
    
    for file in files:
        with open(file, 'r') as f:
            seq = "".join([l.strip().upper() for l in f if not l.startswith('>')])
        
        # Запускаем наш движок
        divergence = analyze_structural_divergence(seq)
        anomaly_count = len(divergence[divergence > (divergence.mean() + 5 * divergence.std())])
        
        results.append({'genome': file, 'anomalies': anomaly_count})
        print(f"Анализ {file}: {anomaly_count} узлов напряжения.")
    
    # Сохраняем итоговый отчет
    pd.DataFrame(results).to_csv('batch_analysis_summary.csv', index=False)
    print("\n--- Отчет готов: batch_analysis_summary.csv ---")

if __name__ == "__main__":
    run_batch_analysis()