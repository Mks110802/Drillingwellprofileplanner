import streamlit as st
import numpy as np
import pandas as pd
import math
import plotly.express as px
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Configuration & Styling
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Advanced Rabia Well Planner",
    page_icon="🛢️",
    layout="wide"
)

# Custom CSS for a professional look
st.markdown("""
<style>
    .metric-card {
        background-color: #f0f2f6;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
        text-align: center;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #f0f2f6;
        border-radius: 4px 4px 0px 0px;
        gap: 1px;
        padding-top: 10px;
        padding-bottom: 10px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #ffffff;
        border-bottom: 2px solid #ff4b4b;
    }
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# 1. Core Physics & Surveying Engine (Minimum Curvature)
# ---------------------------------------------------------------------------

DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi


def calculate_radius_of_curvature(rate_per_100):
    """
    H. Rabia Eq 11.6: R = (360 * 100) / (2 * pi * Rate)
    """
    if rate_per_100 <= 0: return float('inf')
    return (360.0 * 100.0) / (2.0 * math.pi * rate_per_100)


def min_curvature_increment(dmd, inc1_rad, inc2_rad, azi1_rad, azi2_rad):
    """
    Calculates the spatial change between two survey points using the Minimum Curvature Method.
    """
    if dmd <= 0.0: return 0.0, 0.0, 0.0

    cos_dl = (math.cos(inc1_rad) * math.cos(inc2_rad) +
              math.sin(inc1_rad) * math.sin(inc2_rad) * math.cos(azi2_rad - azi1_rad))

    # Numerical stability clamp
    cos_dl = max(-1.0, min(1.0, cos_dl))
    dogleg = math.acos(cos_dl)

    if abs(dogleg) < 1e-6:
        rf = 1.0
    else:
        rf = (2.0 / dogleg) * math.tan(dogleg / 2.0)

    dTVD = (dmd / 2.0) * (math.cos(inc1_rad) + math.cos(inc2_rad)) * rf
    dN = (dmd / 2.0) * (math.sin(inc1_rad) * math.cos(azi1_rad) + math.sin(inc2_rad) * math.cos(azi2_rad)) * rf
    dE = (dmd / 2.0) * (math.sin(inc1_rad) * math.sin(azi1_rad) + math.sin(inc2_rad) * math.sin(azi2_rad)) * rf

    return dTVD, dN, dE


# ---------------------------------------------------------------------------
# 2. Geometric Design Solvers (Based on H. Rabia Text)
# ---------------------------------------------------------------------------

def solve_type_I_rabia(kop, tvd_target, dep_target, bur):
    """
    Solves Type I (Build & Hold) Geometry.
    """
    R = calculate_radius_of_curvature(bur)
    V3_V1 = tvd_target - kop
    D2 = dep_target

    if V3_V1 <= 0: return None, "Target TVD must be deeper than KOP"

    # --- Case 1: R > Target Displacement ---
    if R > D2:
        try:
            dx = D2 - R
            dz = V3_V1
            dist_to_center = math.sqrt(dx ** 2 + dz ** 2)

            if dist_to_center < R:
                return None, "Geometry Impossible: Target is inside the minimum turning radius."

            theta_to_target = math.atan2(dz, dx)
            theta_offset = math.acos(R / dist_to_center)
            alpha_rad = theta_to_target - theta_offset

        except Exception as e:
            return None, f"Calculation Error Case 1: {str(e)}"

    # --- Case 2: R < Target Displacement ---
    else:
        try:
            # Simple geometric fit for large departures
            dx = D2 - R
            dz = V3_V1
            hyp = math.sqrt(dx ** 2 + dz ** 2)
            theta_center = math.atan2(dz, dx)
            theta_off = math.asin(R / hyp)
            alpha_rad = math.pi - (theta_center + (math.pi / 2 - theta_off))

        except:
            return None, "Calculation Error Case 2"

    alpha_deg = alpha_rad * RAD2DEG
    md_build = (alpha_deg * 100) / bur
    tvd_build = R * math.sin(alpha_rad)
    dep_build = R * (1 - math.cos(alpha_rad))
    md_tangent = (V3_V1 - tvd_build) / math.cos(alpha_rad)

    return {
        "alpha": alpha_deg,
        "R": R,
        "md_build": md_build,
        "md_tangent": md_tangent,
        "tvd_eob": kop + tvd_build,
        "dep_eob": dep_build
    }, None


def solve_type_II_rabia(kop, tvd_target, dep_target, bur, dor, end_drop_tvd, final_inc):
    """
    Solves Type II (S-Shape) Geometry.
    """
    R1 = calculate_radius_of_curvature(bur)
    R2 = calculate_radius_of_curvature(dor)

    V1 = kop
    V4 = end_drop_tvd
    D3 = dep_target

    if V4 >= tvd_target: V4 = tvd_target

    def geometry_error(alpha_guess_deg):
        rad = alpha_guess_deg * DEG2RAD
        final_rad = final_inc * DEG2RAD

        # 1. Build
        dz1 = R1 * math.sin(rad)
        dx1 = R1 * (1 - math.cos(rad))

        # 2. Drop
        dz3 = R2 * (math.sin(rad) - math.sin(final_rad))
        dx3 = R2 * (math.cos(final_rad) - math.cos(rad))

        # 3. Tangent
        z_tangent = (V4 - V1) - dz1 - dz3

        if z_tangent < 0: return 1e6

        x_tangent = z_tangent * math.tan(rad)
        total_x = dx1 + x_tangent + dx3

        return abs(total_x - D3)

    # Optimization Loop
    best_alpha = 0
    min_err = 1e9
    scan_range = np.linspace(0.1, 89.0, 500)
    for a in scan_range:
        err = geometry_error(a)
        if err < min_err:
            min_err = err
            best_alpha = a

    if min_err > 10:
        scan_range = np.linspace(best_alpha - 2, best_alpha + 2, 100)
        for a in scan_range:
            err = geometry_error(a)
            if err < min_err:
                min_err = err
                best_alpha = a

    alpha_rad = best_alpha * DEG2RAD
    final_rad = final_inc * DEG2RAD

    dz1 = R1 * math.sin(alpha_rad)
    md_build = (best_alpha * 100) / bur
    dz3 = R2 * (math.sin(alpha_rad) - math.sin(final_rad))
    md_drop = ((best_alpha - final_inc) * 100) / dor
    z_tan = (V4 - V1) - dz1 - dz3
    md_tangent = z_tan / math.cos(alpha_rad)
    z_final = tvd_target - V4
    md_final = z_final / math.cos(final_rad) if z_final > 0 else 0

    return {
        "alpha": best_alpha,
        "R1": R1,
        "R2": R2,
        "md_build": md_build,
        "md_tangent": md_tangent,
        "md_drop": md_drop,
        "md_final": md_final,
        "tvd_eob": kop + dz1,
        "tvd_eod": V4
    }, None


def solve_type_III_rabia(kop, tvd_target, dep_target, bur):
    """
    Type III (Deep Kickoff / J-Type).
    """
    avail_z = tvd_target - kop
    if avail_z <= 0: return None, "KOP is below Target TVD"

    D = dep_target
    R_required = (D ** 2 + avail_z ** 2) / (2 * D)
    bur_required = (360 * 100) / (2 * math.pi * R_required)
    alpha_target_rad = 2 * math.atan(D / avail_z)
    alpha_target_deg = alpha_target_rad * RAD2DEG
    md_build = (alpha_target_deg * 100) / bur_required

    return {
        "R": R_required,
        "bur_actual": bur_required,
        "alpha_target": alpha_target_deg,
        "md_build": md_build
    }, None


# ---------------------------------------------------------------------------
# 3. Trajectory Generation
# ---------------------------------------------------------------------------

def generate_full_survey(inputs, results, type_name):
    surveys = [{
        'md': 0, 'tvd': 0, 'n': 0, 'e': 0,
        'inc': 0, 'azi': inputs['azimuth'],
        'section': 'Vertical', 'dl': 0
    }]

    current_md = 0
    current_inc = 0
    current_azi = inputs['azimuth'] * DEG2RAD

    def march(distance, target_inc_deg, section_label):
        nonlocal current_md, current_inc
        if distance <= 0.1: return

        step_size = 30
        n_steps = int(math.ceil(distance / step_size))
        actual_step = distance / n_steps

        target_inc_rad = target_inc_deg * DEG2RAD
        inc_increment = (target_inc_rad - current_inc) / n_steps

        for _ in range(n_steps):
            next_inc = current_inc + inc_increment
            dtvd, dn, de = min_curvature_increment(
                actual_step, current_inc, next_inc, current_azi, current_azi
            )
            last = surveys[-1]
            dl_deg = (inc_increment / actual_step) * 100 * RAD2DEG if actual_step > 0 else 0

            surveys.append({
                'md': last['md'] + actual_step,
                'tvd': last['tvd'] + dtvd,
                'n': last['n'] + dn,
                'e': last['e'] + de,
                'inc': next_inc * RAD2DEG,
                'azi': inputs['azimuth'],
                'section': section_label,
                'dl': abs(dl_deg)
            })
            current_md += actual_step
            current_inc = next_inc

    kop = inputs['kop']
    march(kop, 0, 'Vertical')

    if type_name == "Type I":
        march(results['md_build'], results['alpha'], 'Build')
        march(results['md_tangent'], results['alpha'], 'Hold (Tangent)')

    elif type_name == "Type II":
        march(results['md_build'], results['alpha'], 'Build')
        march(results['md_tangent'], results['alpha'], 'Hold (Tangent)')
        march(results['md_drop'], inputs['final_inc'], 'Drop')
        if results['md_final'] > 0:
            march(results['md_final'], inputs['final_inc'], 'Final Hold')

    elif type_name == "Type III":
        march(results['md_build'], results['alpha_target'], 'Continuous Build')

    return pd.DataFrame(surveys)


def calculate_cost_estimate(df_survey):
    total_md = df_survey['md'].iloc[-1]
    max_inc = df_survey['inc'].max()
    avg_dl = df_survey['dl'].mean()

    base_cost_per_ft = 350
    inc_factor = 1.0 + (max_inc / 90.0) * 0.4
    tortuosity_factor = 1.0 + (avg_dl / 3.0) * 0.2
    depth_factor = 1.0 + (total_md / 15000.0) * 0.3

    estimated_cost = total_md * base_cost_per_ft * inc_factor * tortuosity_factor * depth_factor
    return estimated_cost


def collision_risk_check(df_survey):
    offset_n, offset_e = 200, 200
    df_survey['dist_to_offset'] = np.sqrt(
        (df_survey['n'] - offset_n) ** 2 + (df_survey['e'] - offset_e) ** 2
    )
    min_dist = df_survey['dist_to_offset'].min()
    return min_dist, offset_n, offset_e


# ---------------------------------------------------------------------------
# 5. Main UI Application
# ---------------------------------------------------------------------------

def main():
    st.title("🛢️ Advanced Rabia Well Planner & Optimizer")
    st.markdown("Combines **H. Rabia's Geometric Design Formulas** with **Minimum Curvature Surveying**.")

    # --- Sidebar Inputs ---
    with st.sidebar:
        st.header("1. Profile Strategy")
        profile_choice = st.selectbox(
            "Select Well Profile",
            ["Type I (Build & Hold)", "Type II (S-Shape)", "Type III (Deep Kickoff)"]
        )

        st.header("2. Geometry Inputs")
        st.markdown("**Enter values to proceed:**")

        # User MUST input values (Defaults are None)
        azimuth = st.number_input("Target Azimuth (°)", min_value=0.0, max_value=360.0, value=None,
                                  placeholder="e.g. 45.0")
        kop = st.number_input("Kick-Off Point (ft)", min_value=0.0, max_value=25000.0, value=None,
                              placeholder="e.g. 2000.0")
        target_tvd = st.number_input("Target TVD (ft)", min_value=100.0, max_value=35000.0, value=None,
                                     placeholder="e.g. 8000.0")
        target_dep = st.number_input("Target Departure (ft)", min_value=0.0, max_value=35000.0, value=None,
                                     placeholder="e.g. 3000.0")

        st.header("3. Operational Constraints")
        bur = st.number_input("Build-Up Rate (°/100ft)", min_value=0.1, max_value=25.0, value=None,
                              placeholder="e.g. 3.0")

        dor = None
        end_drop_tvd = None
        final_inc = 0.0  # Default for Type I/III

        if profile_choice == "Type II (S-Shape)":
            dor = st.number_input("Drop-Off Rate (°/100ft)", min_value=0.1, max_value=25.0, value=None,
                                  placeholder="e.g. 2.0")
            end_drop_tvd = st.number_input("TVD at End of Drop (ft)", min_value=0.0, max_value=35000.0, value=None,
                                           placeholder="Depper than KOP")
            final_inc = st.number_input("Final Inclination (°)", min_value=0.0, max_value=90.0, value=None,
                                        placeholder="e.g. 0.0")

        # Check if inputs are ready
        inputs_ready = (
                azimuth is not None and kop is not None and
                target_tvd is not None and target_dep is not None and
                bur is not None
        )

        if profile_choice == "Type II (S-Shape)":
            inputs_ready = inputs_ready and (dor is not None and end_drop_tvd is not None and final_inc is not None)

        # -----------------------------------------------------------------------
        # NEW SECTION: Credits and About
        # -----------------------------------------------------------------------
        st.markdown("---")
        st.subheader("ℹ️ About & Credits")

        # Project Details
        st.info(
            "**Purpose:** This application automates directional well planning  "
            "using standard geometric methods and Minimum Curvature interpolation.\n\n"
            "**Key Features:**\n"
            "- Trajectory Design (Type I, II, III)\n"
            "- 3D Visualization\n"
            "- Collision Risk Analysis"
        )

        # Team Details
        st.markdown("#### 👨‍💻 Developed By")
        st.markdown("""
        * **Mohit Singh** (Roll No: 23PE3064)
        * **Srijan Srivastava** (Roll No: 23PE3065)
        * **Aman Singh** (Roll No: 23PE3063)
        * **Rajdip Choudhuri** (Roll No: 23PE3066)
        * **Yashashvita Singh** (Roll No: 23PE3061)
        * **Vinay Songara** (Roll No: 23PE3062)
        * **Bedanta Borah** (Roll No: 21ea4pe09)
        """)
        # -----------------------------------------------------------------------

    # --- Calculation Trigger ---
    if not inputs_ready:
        st.info("👈 Please enter all survey parameters in the sidebar to generate the well plan.")
        st.stop()  # Stop execution until data is present

    # --- Calculation Engine ---
    inputs = {'kop': kop, 'azimuth': azimuth, 'target_tvd': target_tvd, 'final_inc': final_inc}
    results = None
    error_msg = None
    type_key = ""

    if profile_choice == "Type I (Build & Hold)":
        results, error_msg = solve_type_I_rabia(kop, target_tvd, target_dep, bur)
        type_key = "Type I"
    elif profile_choice == "Type II (S-Shape)":
        results, error_msg = solve_type_II_rabia(kop, target_tvd, target_dep, bur, dor, end_drop_tvd, final_inc)
        type_key = "Type II"
    elif profile_choice == "Type III (Deep Kickoff)":
        results, error_msg = solve_type_III_rabia(kop, target_tvd, target_dep, bur)
        type_key = "Type III"

    # --- Results Display ---
    if error_msg:
        st.error(f"❌ Design Failed: {error_msg}")
    elif results:
        # Generate the detailed survey
        df_survey = generate_full_survey(inputs, results, type_key)

        # Dashboard - Key Metrics
        total_md = df_survey['md'].iloc[-1]
        final_tvd = df_survey['tvd'].iloc[-1]
        final_dep = math.sqrt(df_survey['n'].iloc[-1] ** 2 + df_survey['e'].iloc[-1] ** 2)

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total Measured Depth", f"{total_md:,.0f} ft")
        col2.metric("Target TVD Reached", f"{final_tvd:,.0f} ft")
        col3.metric("Max Inclination", f"{df_survey['inc'].max():.2f}°")
        col4.metric("Departure Reached", f"{final_dep:,.0f} ft")

        # Tabs
        tab1, tab2, tab3 = st.tabs(["📊 Interactive Plots", "🔬 Engineering Data", "🚀 Innovation Hub"])

        with tab1:
            c1, c2 = st.columns([2, 1])
            with c1:
                fig3d = px.line_3d(
                    df_survey, x='e', y='n', z='tvd', color='section',
                    labels={'e': 'East (+)', 'n': 'North (+)', 'tvd': 'TVD (ft)'},
                    title="3D Well Trajectory Visualization"
                )
                fig3d.update_layout(scene=dict(zaxis=dict(autorange="reversed")))
                fig3d.add_scatter3d(
                    x=[df_survey['e'].iloc[-1]], y=[df_survey['n'].iloc[-1]], z=[df_survey['tvd'].iloc[-1]],
                    mode='markers', marker=dict(size=8, color='red'), name='Target'
                )
                st.plotly_chart(fig3d, use_container_width=True)

            with c2:
                df_survey['vs'] = np.sqrt(df_survey['n'] ** 2 + df_survey['e'] ** 2)
                fig_vs = px.line(
                    df_survey, x='vs', y='tvd', color='section',
                    title="Vertical Section View",
                    labels={'vs': 'Vertical Section (ft)', 'tvd': 'TVD (ft)'}
                )
                fig_vs.update_yaxes(autorange="reversed")
                st.plotly_chart(fig_vs, use_container_width=True)

        with tab2:
            st.subheader("📋 Directional Profile Parameters")

            # Formatted display of Engineering Data (Replacing raw JSON)
            ec1, ec2, ec3 = st.columns(3)

            with ec1:
                st.markdown("**Build Section**")
                if 'R' in results:
                    st.metric("Build Radius (R1)", f"{results['R']:.2f} ft")
                elif 'R1' in results:
                    st.metric("Build Radius (R1)", f"{results['R1']:.2f} ft")

                if 'md_build' in results:
                    st.metric("Build Section Length", f"{results['md_build']:.2f} ft")

            with ec2:
                st.markdown("**Tangent/Hold Section**")
                st.metric("Max Inclination (α)", f"{results.get('alpha', results.get('alpha_target', 0)):.2f}°")
                if 'md_tangent' in results:
                    st.metric("Tangent Length", f"{results['md_tangent']:.2f} ft")
                else:
                    st.metric("Tangent Length", "N/A")

            with ec3:
                st.markdown("**Drop/Final Section**")
                if 'R2' in results:
                    st.metric("Drop Radius (R2)", f"{results['R2']:.2f} ft")
                    st.metric("Drop Section Length", f"{results['md_drop']:.2f} ft")
                else:
                    st.caption("No drop section in this profile.")

            st.divider()

            st.subheader("📐 Governing Drilling Equations")
            st.markdown(
                "The following formulas are utilized based on the Minimum Curvature Method and Geometric Planning:")

            # General Radius Formula
            st.markdown("##### 1. Radius of Curvature")
            st.latex(r"R = \frac{180 \times 100}{\pi \times DLS_{deg}} \approx \frac{5729.58}{DLS}")

            # Specific Profile Formulas
            if type_key == "Type I":
                st.markdown("##### 2. Tangent Section Calculation (Type I)")
                st.latex(r"L_{tan} = \frac{\Delta TVD - R \cdot \sin(\alpha)}{\cos(\alpha)}")
                st.latex(r"\Delta Dep = R \cdot (1 - \cos(\alpha)) + L_{tan} \cdot \sin(\alpha)")

            elif type_key == "Type II":
                st.markdown("##### 2. S-Curve Calculation (Type II)")
                st.latex(r"\Delta TVD_{total} = \Delta TVD_{build} + \Delta TVD_{tangent} + \Delta TVD_{drop}")
                st.latex(
                    r"L_{tangent} = \frac{\Delta TVD - (R_1 \sin \alpha) - R_2 (\sin \alpha - \sin \theta_{final})}{\cos \alpha}")

            st.markdown("##### 3. Dogleg Severity (Minimum Curvature)")
            st.latex(r"\cos(DL) = \sin(I_1)\sin(I_2)\cos(A_2-A_1) + \cos(I_1)\cos(I_2)")
            st.latex(r"DLS = \frac{DL \times 100}{L_{course}} \times \frac{180}{\pi}")

            st.subheader("Survey Data Table")
            st.dataframe(df_survey.head(100), height=300)

            csv = df_survey.to_csv(index=False).encode('utf-8')
            st.download_button("Download Full Survey CSV", csv, "rabia_well_plan.csv", "text/csv")

        with tab3:
            st.subheader("Optimization & Risk")
            ic1, ic2 = st.columns(2)

            with ic1:
                st.markdown("### 💰 Cost Estimator")
                cost = calculate_cost_estimate(df_survey)
                st.metric("Estimated Cost", f"${cost:,.2f}")
                st.caption(f"Cost Factor: Includes complexity multiplier for {type_key} profile.")

            with ic2:
                st.markdown("### ⚠️ Collision Risk")
                min_dist, off_n, off_e = collision_risk_check(df_survey)
                if min_dist < 150:
                    st.error(f"CRITICAL: Path comes within {min_dist:.0f} ft of offset well!")
                else:
                    st.success(f"SAFE: Minimum separation is {min_dist:.0f} ft")
                st.caption(f"Checking against offset well at N:{off_n}, E:{off_e}")


if __name__ == "__main__":
    main()
