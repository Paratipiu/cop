# Simulador de pH — Potencial Eléctrico ⚗️

Simulación teórica de un electrodo de pH básico basado en la **ecuación de Nernst**.

**Universidad Icesi · Física para Ciencias Aplicadas · 2026**  
Sofía García Guevara · Mariana Suárez · Alejandra Anaya

---

## ¿Qué hace?

- Relaciona el voltaje (ΔE) con la concentración de iones H⁺
- Slider interactivo de pH (0–14) y temperatura (0–60 °C)
- 4 gráficas interactivas con Plotly:
  - ΔE vs pH (curva de calibración)
  - ΔE vs log₁₀\[H⁺\]
  - ΔE vs \[H⁺\] escala lineal
  - Pendiente Nernst vs temperatura
- Tabla de valores teóricos simulados
- Marco teórico completo

## Cómo correr localmente

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Cómo desplegarlo en Streamlit Cloud

1. Sube este repositorio a GitHub
2. Ve a [share.streamlit.io](https://share.streamlit.io)
3. Conecta tu cuenta de GitHub
4. Selecciona este repositorio y el archivo `app.py`
5. Click en **Deploy** ✅

## Ecuación principal

```
ΔE = −59.16 × (pH − 7)   [mV, a 25°C]
```

Pendiente nernstiana: **59.16 mV/pH** a 25 °C para especie monovalente (H⁺).
# cop
