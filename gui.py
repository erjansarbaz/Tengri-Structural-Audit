"""
Tengri Structural Audit — Streamlit Frontend
Запуск: streamlit run gui.py
"""

import streamlit as st
import requests
import base64
import json
import os
from io import BytesIO

# ── Конфиг ───────────────────────────────────────────────────────────────────
API_URL = os.getenv("API_URL", "http://api:8000")   # Docker: api:8000, локально: localhost:8000

st.set_page_config(
    page_title="Tengri Structural Audit",
    page_icon="🧬",
    layout="wide",
)

# ── Стиль ────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main { background-color: #F8FAFD; }
    .stButton > button {
        background-color: #1565C0; color: white;
        border-radius: 8px; border: none;
        padding: 0.5rem 2rem; font-size: 1rem;
    }
    .metric-card {
        background: #0D1B30; border-radius: 10px;
        padding: 1rem; text-align: center; color: white;
    }
</style>
""", unsafe_allow_html=True)

# ── Заголовок ─────────────────────────────────────────────────────────────────
col_logo, col_title = st.columns([1, 6])
with col_logo:
    st.markdown("# 🧬")
with col_title:
    st.markdown("## Tengri Structural Audit (TSA)")
    st.caption("Windowed Thermodynamic Approximation · DOI: 10.5281/zenodo.21203626")

st.divider()

# ── Режим ────────────────────────────────────────────────────────────────────
mode = st.radio(
    "Режим анализа:",
    ["🔬 Одиночный анализ", "⚖️ Сравнение: Референс vs Образец"],
    horizontal=True
)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# РЕЖИМ 1: ОДИНОЧНЫЙ АНАЛИЗ
# ══════════════════════════════════════════════════════════════════════════════
if mode == "🔬 Одиночный анализ":
    st.markdown("### Загрузите FASTA-файл")
    st.info("Поддерживаются файлы .fasta / .fna / .fa · Максимум 200 МБ · **2 бесплатных анализа**")

    uploaded = st.file_uploader("Выберите файл:", type=["fasta","fna","fa","txt"])

    if uploaded:
        st.success(f"✅ Файл загружен: **{uploaded.name}** ({uploaded.size/1024:.1f} КБ)")

        if st.button("🚀 Запустить анализ"):
            with st.spinner("Анализируем последовательность..."):
                try:
                    resp = requests.post(
                        f"{API_URL}/analyze",
                        files={"file": (uploaded.name, uploaded.getvalue(), "text/plain")},
                        timeout=120
                    )
                    if resp.status_code == 429:
                        st.error("⛔ Бесплатный лимит исчерпан (2 анализа). "
                                 "Напишите нам для расширения доступа.")
                    elif resp.status_code != 200:
                        st.error(f"Ошибка API: {resp.status_code} — {resp.text}")
                    else:
                        data = resp.json()
                        st.success("✅ Анализ завершён!")

                        # Метрики
                        s = data["stats"]
                        c1,c2,c3,c4 = st.columns(4)
                        c1.metric("Длина", f"{data['seq_len']:,} нт")
                        c2.metric("Max DI", s["max_di"])
                        c3.metric("SCZ зон", s["scz_count"])
                        c4.metric("Threshold", s["threshold"])

                        remaining = data.get("remaining", 0)
                        if remaining == 0:
                            st.warning("⚠️ Это был ваш последний бесплатный анализ.")
                        else:
                            st.info(f"Осталось бесплатных анализов: **{remaining}**")

                        # График
                        st.markdown("### 📊 Профиль структурного конфликта")
                        img_bytes = base64.b64decode(data["chart_b64"])
                        st.image(img_bytes, use_column_width=True)

                        # Скачать
                        col_dl1, col_dl2 = st.columns(2)
                        with col_dl1:
                            st.download_button(
                                "⬇️ Скачать PNG",
                                data=img_bytes,
                                file_name=f"TSA_{uploaded.name}.png",
                                mime="image/png"
                            )
                        with col_dl2:
                            csv_bytes = base64.b64decode(data["csv_b64"])
                            st.download_button(
                                "⬇️ Скачать CSV",
                                data=csv_bytes,
                                file_name=f"TSA_{uploaded.name}.csv",
                                mime="text/csv"
                            )

                except requests.exceptions.ConnectionError:
                    st.error("❌ Не удалось подключиться к API. Убедитесь что сервер запущен.")

# ══════════════════════════════════════════════════════════════════════════════
# РЕЖИМ 2: СРАВНЕНИЕ
# ══════════════════════════════════════════════════════════════════════════════
else:
    st.markdown("### Сравнение: Референс vs Ваш образец")
    st.info("""
    **Как это работает:**
    1. Загрузите **референсный геном** (например, стандартный E. coli K-12) → движок выставит его как стандарт
    2. Загрузите **ваш образец** (ваш штамм E. coli) → движок сравнит его с референсом
    3. Получите сравнительный график — где ваш штамм отклоняется от нормы
    """)

    col_ref, col_smp = st.columns(2)
    with col_ref:
        st.markdown("#### 📌 Референс (стандарт)")
        ref_file = st.file_uploader("Референсный геном:", type=["fasta","fna","fa","txt"],
                                     key="ref")
        if ref_file:
            st.success(f"✅ {ref_file.name} ({ref_file.size/1024:.1f} КБ)")

    with col_smp:
        st.markdown("#### 🔬 Ваш образец")
        smp_file = st.file_uploader("Ваш образец для анализа:", type=["fasta","fna","fa","txt"],
                                     key="smp")
        if smp_file:
            st.success(f"✅ {smp_file.name} ({smp_file.size/1024:.1f} КБ)")

    if ref_file and smp_file:
        if st.button("⚖️ Запустить сравнительный анализ"):
            with st.spinner("Сравниваем профили..."):
                try:
                    resp = requests.post(
                        f"{API_URL}/compare",
                        files={
                            "reference": (ref_file.name, ref_file.getvalue(), "text/plain"),
                            "sample":    (smp_file.name, smp_file.getvalue(), "text/plain"),
                        },
                        timeout=180
                    )
                    if resp.status_code == 429:
                        st.error("⛔ Бесплатный лимит исчерпан. Напишите нам для продолжения.")
                    elif resp.status_code != 200:
                        st.error(f"Ошибка: {resp.status_code} — {resp.text}")
                    else:
                        data = resp.json()
                        st.success("✅ Сравнение завершено!")

                        # Метрики сравнения
                        st.markdown("### 📊 Результаты сравнения")
                        c1,c2,c3,c4 = st.columns(4)
                        c1.metric("Ref Max DI",    data["reference"]["stats"]["max_di"])
                        c2.metric("Sample Max DI", data["sample"]["stats"]["max_di"],
                                  delta=str(data["delta_max_di"]))
                        c3.metric("Ref SCZ",    data["reference"]["stats"]["scz_count"])
                        c4.metric("Sample SCZ", data["sample"]["stats"]["scz_count"],
                                  delta=str(data["delta_scz"]))

                        remaining = data.get("remaining", 0)
                        if remaining == 0:
                            st.warning("⚠️ Это был ваш последний бесплатный анализ.")
                        else:
                            st.info(f"Осталось бесплатных анализов: **{remaining}**")

                        # График
                        img_bytes = base64.b64decode(data["chart_b64"])
                        st.image(img_bytes, use_column_width=True)

                        # Скачать
                        col_dl1, col_dl2 = st.columns(2)
                        with col_dl1:
                            st.download_button(
                                "⬇️ Скачать PNG",
                                data=img_bytes,
                                file_name="TSA_comparison.png",
                                mime="image/png"
                            )
                        with col_dl2:
                            csv_bytes = base64.b64decode(data["csv_b64"])
                            st.download_button(
                                "⬇️ Скачать CSV (образец)",
                                data=csv_bytes,
                                file_name="TSA_sample.csv",
                                mime="text/csv"
                            )

                except requests.exceptions.ConnectionError:
                    st.error("❌ Не удалось подключиться к API.")
    else:
        st.warning("⬆️ Загрузите оба файла для запуска сравнения")

# ── Футер ────────────────────────────────────────────────────────────────────
st.divider()
st.markdown(
    "**Tengri Lab** · Windowed Thermodynamic Approximation (WTA) · "
    "[DOI: 10.5281/zenodo.21203626](https://doi.org/10.5281/zenodo.21203626) · "
    "[GitHub](https://github.com/erjansarbaz/Tengri-Structural-Audit) · "
    "По вопросам расширения доступа: ybaynazarovarchitect@gmail.com",
    unsafe_allow_html=False
)
