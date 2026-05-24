"""
Simulador de pH mediante Potencial Eléctrico
Basado en la ecuación de Nernst
Sofía García Guevara · Mariana Suárez · Alejandra Anaya
Universidad Icesi — Física para Ciencias Aplicadas 2026
"""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CONFIGURACIÓN DE PÁGINA ─────────────────────────────────────────
st.set_page_config(
    page_title="Simulador pH — Potencial Eléctrico",
    page_icon="⚗️",
    layout="wide",
)

# ── ESTILOS CSS ─────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

html, body, [class*="css"] { font-family: 'Space Grotesk', sans-serif; }

.main-title {
    font-size: 2.2rem;
    font-weight: 700;
    background: linear-gradient(135deg, #e2e8f0 30%, #63b3ed 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    margin-bottom: 0.2rem;
    letter-spacing: -0.02em;
}
.subtitle {
    color: #718096;
    font-size: 0.95rem;
    margin-bottom: 0.2rem;
}
.meta { color: #4a5568; font-size: 0.82rem; margin-bottom: 1.5rem; }
.badge {
    display: inline-block;
    background: rgba(99,179,237,0.1);
    border: 1px solid rgba(99,179,237,0.35);
    color: #63b3ed;
    padding: 2px 12px;
    border-radius: 20px;
    font-size: 0.75rem;
    font-weight: 600;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    margin-bottom: 0.8rem;
}
.metric-card {
    background: #111827;
    border: 1px solid #243049;
    border-radius: 10px;
    padding: 1rem 1.2rem 0.8rem;
    margin-bottom: 0.6rem;
}
.metric-label {
    font-size: 0.72rem;
    color: #718096;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    font-weight: 600;
    margin-bottom: 4px;
}
.metric-value {
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.4rem;
    font-weight: 600;
    line-height: 1.2;
}
.metric-unit {
    font-size: 0.7rem;
    color: #4a5568;
    margin-top: 2px;
}
.equation-box {
    background: #111827;
    border: 1px solid #243049;
    border-left: 3px solid #63b3ed;
    border-radius: 8px;
    padding: 0.9rem 1.1rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.82rem;
    color: #a0aec0;
    margin-bottom: 1rem;
    line-height: 1.9;
}
.eq-live {
    font-size: 0.95rem;
    font-weight: 600;
}
.theory-box {
    background: #111827;
    border: 1px solid #243049;
    border-radius: 10px;
    padding: 1.4rem;
    margin-top: 0.5rem;
}
.theory-box h4 { color: #63b3ed; font-size: 0.95rem; margin-bottom: 0.5rem; }
.theory-box p  { color: #a0aec0; font-size: 0.85rem; line-height: 1.7; }
.formula {
    background: #1a2235;
    border-left: 3px solid #f6ad55;
    border-radius: 6px;
    padding: 0.7rem 1rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.82rem;
    color: #f6ad55;
    margin: 0.6rem 0;
    line-height: 1.9;
}
.stSlider > div > div > div > div { background: #63b3ed !important; }
div[data-testid="stMetric"] {
    background: #111827;
    border: 1px solid #243049;
    border-radius: 10px;
    padding: 0.8rem 1rem;
}
</style>
""", unsafe_allow_html=True)

# ── CONSTANTES FÍSICAS ───────────────────────────────────────────────
R      = 8.314
F_CONST = 96485
PH_REF  = 7

# ── FUNCIONES ────────────────────────────────────────────────────────
def pendiente_nernst(T_C):
    return (2.303 * R * (T_C + 273.15) / F_CONST) * 1000

def delta_E(pH, T_C=25):
    return -pendiente_nernst(T_C) * (pH - PH_REF)

def h_conc(pH):
    return 10 ** (-pH)

def estado_solucion(pH):
    if pH < 4:    return "Muy ácido", "#fc8181"
    elif pH < 7:  return "Ácido",     "#f6ad55"
    elif pH == 7: return "Neutro",    "#68d391"
    elif pH <= 10: return "Básico",   "#63b3ed"
    else:          return "Muy básico", "#b794f4"

def color_delta(dE):
    if dE > 0.5:   return "#fc8181"
    elif dE < -0.5: return "#63b3ed"
    return "#68d391"

# ── HEADER ───────────────────────────────────────────────────────────
st.markdown('<div class="badge">Simulación Teórica · Física para Ciencias Aplicadas</div>', unsafe_allow_html=True)
st.markdown('<div class="main-title">Medición de pH mediante Potencial Eléctrico</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitle">Simulación de un electrodo de pH básico basado en la ecuación de Nernst</div>', unsafe_allow_html=True)
st.markdown('<div class="meta">Sofía García Guevara · Mariana Suárez · Alejandra Anaya — Universidad Icesi · 2026</div>', unsafe_allow_html=True)

st.divider()

# ── LAYOUT: columnas ─────────────────────────────────────────────────
col_ctrl, col_main = st.columns([1, 2.6], gap="large")

with col_ctrl:
    st.markdown("#### Control del simulador")

    pH = st.slider("**pH de la solución**", 0.0, 14.0, 7.0, 0.1,
                   help="Desliza para cambiar el pH de 0 (muy ácido) a 14 (muy básico)")

    T_C = st.slider("**Temperatura (°C)**", 0, 60, 25, 1,
                    help="Cambia la temperatura y observa cómo varía la pendiente de Nernst")

    # Cálculos
    slope  = pendiente_nernst(T_C)
    dE     = delta_E(pH, T_C)
    hc     = h_conc(pH)
    lh     = -pH
    estado, clr_estado = estado_solucion(pH)
    sign   = "+" if dE > 0 else ""
    clr_dE = color_delta(dE)

    st.markdown(f"""
    <div style="margin: 0.8rem 0 1rem;">
        <span style="background:{'rgba(252,129,129,0.15)' if pH<7 else 'rgba(104,211,145,0.15)' if pH==7 else 'rgba(99,179,237,0.15)'};
                     border:1px solid {clr_estado}40;
                     color:{clr_estado};
                     padding:4px 14px;
                     border-radius:20px;
                     font-size:0.8rem;
                     font-weight:600;">
            ● {estado}
        </span>
    </div>
    """, unsafe_allow_html=True)

    # Métricas
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-label">ΔE — Potencial</div>
        <div class="metric-value" style="color:{clr_dE}">{sign}{dE:.2f}</div>
        <div class="metric-unit">mV</div>
    </div>
    <div class="metric-card">
        <div class="metric-label">[H⁺] Concentración</div>
        <div class="metric-value" style="color:#f6ad55">10⁻{pH:.1f}</div>
        <div class="metric-unit">mol/L</div>
    </div>
    <div class="metric-card">
        <div class="metric-label">log₁₀[H⁺]</div>
        <div class="metric-value" style="color:#b794f4">−{pH:.1f}</div>
        <div class="metric-unit">adimensional</div>
    </div>
    <div class="metric-card">
        <div class="metric-label">Pendiente Nernst</div>
        <div class="metric-value" style="color:#63b3ed">{slope:.2f}</div>
        <div class="metric-unit">mV/pH · a {T_C}°C</div>
    </div>
    """, unsafe_allow_html=True)

    # Ecuación activa
    st.markdown(f"""
    <div class="equation-box">
        <div style="color:#4a5568;font-size:0.7rem;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:4px">
            Ecuación activa
        </div>
        <div>ΔE = −{slope:.2f} × (pH − 7)</div>
        <div class="eq-live" style="color:{clr_dE};margin-top:6px">
            ΔE = −{slope:.2f} × ({pH:.1f} − 7) = {sign}{dE:.2f} mV
        </div>
    </div>
    """, unsafe_allow_html=True)

# ── GRÁFICAS ─────────────────────────────────────────────────────────
with col_main:
    tabs = st.tabs(["📈 Gráficas", "📋 Tabla de datos", "📖 Marco teórico"])

    # ── TAB 1: GRÁFICAS ──────────────────────────────────────────────
    with tabs[0]:
        pH_range  = np.linspace(0, 14, 500)
        dE_range  = delta_E(pH_range, T_C)
        logH_range = -pH_range
        temps_range = np.linspace(0, 60, 200)
        slopes_range = [pendiente_nernst(t) for t in temps_range]

        # Gráfica 1: ΔE vs pH
        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(
            x=pH_range, y=dE_range,
            mode='lines', name='Curva teórica',
            line=dict(color='#63b3ed', width=2.5),
            fill='tozeroy',
            fillcolor='rgba(99,179,237,0.07)'
        ))
        fig1.add_trace(go.Scatter(
            x=pH_range[pH_range < 7], y=dE_range[pH_range < 7],
            fill='tozeroy', fillcolor='rgba(252,129,129,0.1)',
            mode='none', name='Zona ácida (ΔE > 0)', showlegend=True
        ))
        fig1.add_trace(go.Scatter(
            x=pH_range[pH_range > 7], y=dE_range[pH_range > 7],
            fill='tozeroy', fillcolor='rgba(99,179,237,0.1)',
            mode='none', name='Zona básica (ΔE < 0)', showlegend=True
        ))
        fig1.add_trace(go.Scatter(
            x=[pH], y=[dE],
            mode='markers', name=f'pH actual = {pH:.1f}',
            marker=dict(color=clr_dE, size=13, symbol='circle',
                        line=dict(color='white', width=2))
        ))
        fig1.add_hline(y=0, line_dash='dash', line_color='#68d391', line_width=1.2, opacity=0.7)
        fig1.add_vline(x=7, line_dash='dash', line_color='#68d391', line_width=1.2, opacity=0.5)
        fig1.update_layout(
            title=dict(text='ΔE vs pH — Curva de calibración teórica', font_size=13, font_color='#a0aec0'),
            xaxis=dict(title='pH', range=[0,14], dtick=1, gridcolor='#1e2d45', color='#718096'),
            yaxis=dict(title='ΔE (mV)', range=[-460,460], gridcolor='#1e2d45', color='#718096'),
            paper_bgcolor='#0b0f1a', plot_bgcolor='#111827',
            legend=dict(font_size=10, bgcolor='#1a2235', bordercolor='#243049'),
            margin=dict(l=50,r=20,t=40,b=40), height=280,
            font=dict(family='Space Grotesk')
        )
        st.plotly_chart(fig1, use_container_width=True)

        c1, c2 = st.columns(2)

        # Gráfica 2: ΔE vs log10[H+]
        with c1:
            fig2 = go.Figure()
            fig2.add_trace(go.Scatter(
                x=logH_range, y=dE_range,
                mode='lines', name='ΔE vs log₁₀[H⁺]',
                line=dict(color='#b794f4', width=2.5),
                fill='tozeroy', fillcolor='rgba(183,148,244,0.07)'
            ))
            fig2.add_trace(go.Scatter(
                x=[lh], y=[dE],
                mode='markers', name=f'log₁₀[H⁺] = {lh:.1f}',
                marker=dict(color='#b794f4', size=12, symbol='circle',
                            line=dict(color='white', width=2))
            ))
            fig2.add_hline(y=0, line_dash='dash', line_color='#68d391', line_width=1, opacity=0.6)
            fig2.add_vline(x=-7, line_dash='dash', line_color='#68d391', line_width=1, opacity=0.5)
            fig2.update_layout(
                title=dict(text='ΔE vs log₁₀[H⁺]', font_size=12, font_color='#a0aec0'),
                xaxis=dict(title='log₁₀[H⁺]', range=[-14,0], dtick=2, gridcolor='#1e2d45', color='#718096'),
                yaxis=dict(title='ΔE (mV)', range=[-460,460], gridcolor='#1e2d45', color='#718096'),
                paper_bgcolor='#0b0f1a', plot_bgcolor='#111827',
                showlegend=False, margin=dict(l=50,r=10,t=40,b=40), height=240,
                font=dict(family='Space Grotesk')
            )
            st.plotly_chart(fig2, use_container_width=True)

        # Gráfica 3: ΔE vs [H+] lineal (pH 5-14)
        with c2:
            ph_lin    = np.linspace(5, 14, 400)
            hconc_lin = h_conc(ph_lin)
            dE_lin    = delta_E(ph_lin, T_C)

            fig3 = go.Figure()
            fig3.add_trace(go.Scatter(
                x=hconc_lin, y=dE_lin,
                mode='lines', name='ΔE vs [H⁺]',
                line=dict(color='#f6ad55', width=2.5),
                fill='tozeroy', fillcolor='rgba(246,173,85,0.07)'
            ))
            if 5 <= pH <= 14:
                fig3.add_trace(go.Scatter(
                    x=[hc], y=[dE],
                    mode='markers', name=f'[H⁺] = 10⁻{pH:.1f}',
                    marker=dict(color='#f6ad55', size=12, symbol='circle',
                                line=dict(color='white', width=2))
                ))
            fig3.add_hline(y=0, line_dash='dash', line_color='#68d391', line_width=1, opacity=0.6)
            fig3.update_layout(
                title=dict(text='ΔE vs [H⁺] escala lineal (pH 5–14)', font_size=12, font_color='#a0aec0'),
                xaxis=dict(title='[H⁺] (mol/L)', gridcolor='#1e2d45', color='#718096',
                           tickformat='.1e'),
                yaxis=dict(title='ΔE (mV)', range=[-460,460], gridcolor='#1e2d45', color='#718096'),
                paper_bgcolor='#0b0f1a', plot_bgcolor='#111827',
                showlegend=False, margin=dict(l=50,r=10,t=40,b=40), height=240,
                font=dict(family='Space Grotesk')
            )
            st.plotly_chart(fig3, use_container_width=True)

        # Gráfica 4: Pendiente Nernst vs T
        fig4 = go.Figure()
        fig4.add_trace(go.Scatter(
            x=list(temps_range), y=slopes_range,
            mode='lines', name='Pendiente Nernst',
            line=dict(color='#fc8181', width=2.5),
            fill='tozeroy', fillcolor='rgba(252,129,129,0.07)'
        ))
        fig4.add_trace(go.Scatter(
            x=[T_C], y=[slope],
            mode='markers', name=f'T = {T_C}°C → {slope:.2f} mV/pH',
            marker=dict(color='#fc8181', size=12, symbol='circle',
                        line=dict(color='white', width=2))
        ))
        fig4.add_annotation(
            x=T_C, y=slope,
            text=f"  {slope:.2f} mV/pH",
            showarrow=False, xanchor='left',
            font=dict(color='#68d391', size=11)
        )
        fig4.update_layout(
            title=dict(text='Pendiente de Nernst vs Temperatura', font_size=12, font_color='#a0aec0'),
            xaxis=dict(title='Temperatura (°C)', range=[0,60], gridcolor='#1e2d45', color='#718096'),
            yaxis=dict(title='Pendiente (mV/pH)', gridcolor='#1e2d45', color='#718096'),
            paper_bgcolor='#0b0f1a', plot_bgcolor='#111827',
            legend=dict(font_size=10, bgcolor='#1a2235', bordercolor='#243049'),
            margin=dict(l=50,r=20,t=40,b=40), height=240,
            font=dict(family='Space Grotesk')
        )
        st.plotly_chart(fig4, use_container_width=True)

    # ── TAB 2: TABLA ─────────────────────────────────────────────────
    with tabs[1]:
        st.markdown("##### Valores teóricos simulados (T = 25 °C, referencia: pH 7)")
        table_ph = list(range(0, 15))
        rows = []
        for ph_t in table_ph:
            dE_t   = delta_E(ph_t, 25)
            sign_t = "+" if dE_t > 0 else ""
            est_t, _ = estado_solucion(ph_t)
            rows.append({
                "pH": ph_t,
                "[H⁺] (mol/L)": f"1.0 × 10⁻{ph_t}",
                "log₁₀[H⁺]": f"−{ph_t}",
                "ΔE (mV)": f"{sign_t}{dE_t:.2f}",
                "Estado": est_t
            })
        import pandas as pd
        df = pd.DataFrame(rows)

        def highlight_row(row):
            ph_r = row['pH']
            dE_r = delta_E(ph_r, 25)
            base = [''] * len(row)
            if abs(ph_r - pH) < 0.5:
                return ['background-color: #1a2e4a; font-weight: bold'] * len(row)
            if dE_r > 0:
                return ['color: #fc8181'] * len(row)
            elif dE_r < 0:
                return ['color: #63b3ed'] * len(row)
            return ['color: #68d391; font-weight: bold'] * len(row)

        st.dataframe(
            df.style.apply(highlight_row, axis=1),
            use_container_width=True,
            height=500
        )

    # ── TAB 3: MARCO TEÓRICO ─────────────────────────────────────────
    with tabs[2]:
        st.markdown("""
        <div class="theory-box">
            <h4>1. Definición de pH y concentración de H⁺</h4>
            <p>El pH se define como el logaritmo decimal negativo de la actividad del ion hidrógeno.
            En soluciones diluidas se aproxima la actividad por concentración molar:</p>
            <div class="formula">
                pH = −log₁₀(aH⁺)<br>
                [H⁺] ≈ 10⁻ᵖᴴ  (mol/L)
            </div>
            <p>Una diferencia de una unidad de pH equivale a un cambio de <strong>10 veces</strong>
            en la concentración de H⁺.</p>

            <h4 style="margin-top:1.2rem">2. Potencial eléctrico y diferencia de potencial</h4>
            <p>La diferencia de potencial entre dos electrodos (sensible y de referencia) depende
            de los equilibrios de carga en la interfaz electrodo-solución:</p>
            <div class="formula">ΔV = V_A − V_B = W / q</div>

            <h4 style="margin-top:1.2rem">3. Ecuación de Nernst</h4>
            <p>Relaciona el potencial del electrodo con la actividad de las especies iónicas.
            Para H⁺ (especie monovalente, n = 1) a temperatura T:</p>
            <div class="formula">
                E = E° − (RT/nF) · ln Q<br><br>
                Para H⁺:  E = K − (2.303RT/F) · pH<br><br>
                A 25 °C:  E(mV) = K − <strong>59.16</strong> · pH
            </div>

            <h4 style="margin-top:1.2rem">4. Modelo de simulación</h4>
            <p>Tomando pH 7 como referencia (ΔE = 0 mV):</p>
            <div class="formula">
                ΔE = −59.16 × (pH − 7)  [mV]<br><br>
                Equivalentemente:<br>
                ΔE = 59.16 × log₁₀([H⁺] / 10⁻⁷)  [mV]
            </div>

            <h4 style="margin-top:1.2rem">5. Conclusión clave</h4>
            <p>El voltaje es <strong>lineal con el pH</strong> y con <strong>log₁₀[H⁺]</strong>,
            pero <em>no</em> con la concentración directa [H⁺] en escala lineal. Por eso un pHmetro
            convierte la señal eléctrica directamente a pH mediante un modelo electroquímico,
            no a concentración molar.</p>
        </div>
        """, unsafe_allow_html=True)

st.divider()
st.markdown(
    "<div style='text-align:center;color:#4a5568;font-size:0.75rem'>"
    "Simulación teórica · datos calculados con la ecuación de Nernst · "
    "Universidad Icesi · Ingeniería Bioquímica · 2026"
    "</div>",
    unsafe_allow_html=True
)
