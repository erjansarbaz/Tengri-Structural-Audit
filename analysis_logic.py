import numpy as np

def analyze_structural_divergence(sequence, window_size=1000):
    # Кодирование
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4}
    arr = np.array([mapping.get(n, 4) for n in sequence], dtype=int)
    
    # Веса связей
    weights = np.array([-4.5, 1.5, 0.8, 4.2, 0.1, 4.6, -6.1, 2.2, 1.0, 0.1, 
                        0.5, 3.1, -5.2, 1.8, 0.1, 2.5, 0.9, 3.8, -4.8, 0.1, 
                        0.1, 0.1, 0.1, 0.1, 0.1], dtype=np.float32)
    
    divergence = np.zeros(len(arr))
    
    # Скользящее окно (убирает глобальный дрейф)
    for i in range(len(arr) - window_size):
        window = arr[i : i + window_size]
        traj = window[:-1] * 5 + window[1:]
        energy = weights[traj]
        
        s_f = np.cumsum(energy)
        s_b = np.cumsum(energy[::-1])[::-1]
        divergence[i + window_size // 2] = np.mean(np.abs(s_f - s_b))
        
    return divergence