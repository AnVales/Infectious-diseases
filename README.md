# 🦠 Infectious-Diseases Simulation 💉

Simulaciones de epidemias usando **modelos SIR** en Python, con y sin vacunación.  
¡Observa cómo evoluciona la población en diferentes escenarios! 📊

---

## 🔹 Población
- 👥 **Susceptibles (S)**  
- 🤒 **Infectados (I)**  
- 💪 **Recuperados con inmunidad (R)**  
- 🛡️ **Vacunados (V)**

---

## 🔹 Parámetros
- ⚡ **Beta (β):** tasa de transmisión  
- 💚 **Gamma (γ):** tasa de recuperación  
- 💉 **p:** porcentaje de vacunación diaria  

---

## 🔹 Escenarios
1. **SIR clásico** (sin vacuna) 🏥  
2. **SIR + vacuna constante** 💉  
3. **SIR + vacunación desde cierto día** 📅  
4. **SIR + vacunación proporcional diaria** 📈  

---

## 🔹 Archivos principales
| Archivo | Descripción |
|---------|-------------|
| `SIR_vaccine.py` / `SIR_vaccine_2.py` 🧪 | Simulación SIR con vacuna |
| `SIRdays.py` 📈 | Evolución temporal de la epidemia |
| `npopulation.py` 🧮 | Diferentes estrategias de vacunación |
| `beta.py`, `gamma.py`, `gamma_beta.py` ⚙️ | Ajustes de parámetros |
| `im_leyend.py`, `im_leyend_vaccine.py` 🎨 | Visualización de resultados |

---

## 🔹 Cómo ejecutar
1. Instala las dependencias:  
```bash
pip install numpy scipy matplotlib
