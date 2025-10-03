# -------------------------------
# Epitope Prediction App with 3D Structure and Disease Epitope Highlighting
# -------------------------------

import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
from Bio.PDB import PDBList

# -------------------------------
# Epitope Prediction Function
# -------------------------------
def predict_epitopes(seq, window=9, threshold=0.3):
    hydrophilic = set("KRDEQN")
    epitopes = []
    for i in range(len(seq) - window + 1):
        peptide = seq[i:i + window]
        score = sum(aa in hydrophilic for aa in peptide) / window
        if score >= threshold:
            epitopes.append({
                "peptide": peptide,
                "start": i + 1,
                "end": i + window,
                "score": round(score, 2)
            })
    return epitopes

# -------------------------------
# Streamlit Web App
# -------------------------------
st.title("🧬 Epitope Prediction App with 3D Structure")
st.write("Enter a protein sequence and optionally a PDB code to predict epitopes and highlight them in 3D.")

# User inputs
seq = st.text_area("Protein sequence (use single-letter amino acid codes):")
pdb_code = st.text_input("Enter PDB code for structure visualization (e.g., 1TUP):")

threshold = st.slider("Epitope score threshold", min_value=0.1, max_value=1.0, value=0.3, step=0.1)

# Prediction button
if st.button("Predict Epitopes"):
    if not seq:
        st.warning("Please enter a protein sequence")
    else:
        epitopes = predict_epitopes(seq.upper(), threshold=threshold)
        
        if not epitopes:
            st.info("No epitopes found with the given threshold.")
        else:
            st.success(f"Found {len(epitopes)} epitope(s):")
            for epi in epitopes:
                st.write(f"Peptide: {epi['peptide']}, Start: {epi['start']}, End: {epi['end']}, Score: {epi['score']}")

        # -------------------------------
        # 3D Structure Visualization
        # -------------------------------
        if pdb_code:
            st.subheader("3D Structure with Highlighted Epitopes")

            try:
                pdbl = PDBList()
                pdb_file = pdbl.retrieve_pdb_file(pdb_code, pdir=".", file_format="pdb")

                with open(pdb_file, "r") as f:
                    pdb_data = f.read()

                view = py3Dmol.view(width=800, height=500)
                view.addModel(pdb_data, "pdb")
                view.setStyle({"cartoon": {"color": "white"}})

                # Highlight epitopes
                for epi in epitopes:
                    resi_range = list(range(epi["start"], epi["end"] + 1))
                    view.addStyle({"resi": resi_range}, {"stick": {"colorscheme": "redCarbon"}})

                view.zoomTo()

                # Render structure in Streamlit
                components.html(view._make_html(), height=500)

            except Exception as e:
                st.error(f"Error loading structure: {e}")

        # -------------------------------
        # Disease-causing epitope example
        # -------------------------------
        st.subheader("Known Disease-Causing Epitopes Example")
        st.write("""
        Example: HIV-1 gp41 Epitope recognized by 2F5 antibody  
        Sequence: ELDKWAS  
        Location: Fusion peptide region of gp41  
        Disease significance: Critical for HIV entry into host cells, targeted by neutralizing antibodies.
        """)
