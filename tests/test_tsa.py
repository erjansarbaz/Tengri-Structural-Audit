import numpy as np
from analysis_logic import analyze_structural_divergence

def test_analyze_structural_divergence():
    # Создаем простую искусственную нуклеотидную последовательность
    test_seq = "ATGC" * 500  # 2000 нуклеотидов
    window_size = 100
    
    # Запускаем функцию
    divergence = analyze_structural_divergence(test_seq, window_size=window_size)
    
    # Проверяем, что результат — это numpy массив
    assert isinstance(divergence, np.ndarray), "Result must be a numpy ndarray"
    
    # Проверяем, что длина массива совпадает с длиной исходной последовательности
    assert len(divergence) == len(test_seq), "Divergence array length must match sequence length"
    
    # Проверяем, что значения рассчитаны (не все нули)
    assert not np.all(divergence == 0), "Divergence values should not be all zeros"
