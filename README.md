# 🛢️ Advanced Rabia Well Planner

An automated **Directional Well Planning Tool** developed in Python using Streamlit. This application integrates **H. Rabia’s geometric design formulas** with the industry-standard **Minimum Curvature Method** to generate, visualize, and optimize precise well trajectories.

## 🚀 Features

* **Multi-Profile Support:** Instantly designs Type I (J-curve), Type II (S-curve), and Type III (Deep Kickoff) wells.
* **Physics-Informed:** Uses minimum curvature calculations for accurate Dogleg Severity (DLS) and survey interpolation.
* **3D Visualization:** Interactive 3D and 2D vertical section plotting using Plotly.
* **Risk & Economics:** Includes a drilling cost estimator and offset well collision risk analysis.
* **Data Export:** Download full survey data (MD, Inc, Azi, TVD, N, E) as CSV.

## 🛠️ Installation & Usage

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/rabia-well-planner.git](https://github.com/YOUR_USERNAME/rabia-well-planner.git)
    cd rabia-well-planner
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the application:**
    ```bash
    streamlit run app.py
    ```

## 📐 Methodology

The tool uses a two-step physics engine:
1.  **Geometric Solver:** Calculates critical tie-in points (KOP, EOB, EOD) based on H. Rabia's *Oilwell Drilling Engineering*.
2.  **Survey Engine:** Interpolates the path using the Minimum Curvature method for realistic wellbore representation.

## 👨‍💻 Team

Developed by **Mohit Singh**:
* **Developer 1**


---
*Built with Python 🐍 and Streamlit 🎈*
