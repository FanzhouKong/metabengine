# Author: Hauxu Yu

# A module to support the data processing workflow by LC-Binbase

# Workflow
# 1. Load internal standard library with no retention time and MS/MS spectra
# 2. Get all method blank data
# 3. Run peak picking of method blank data
# 4. Targeted search for internal standards in method blank data
# 5. For observed internal standards, get retention time and MS/MS spectra, and add to the library (project-specific)
# 6. Generate aligned feature table for method blank data
# 7. Get all pooled quality control sample data
# 8. Run peak picking of pooled QC and correct retention time using internal standards
# 8. Align pooled QCs
# 9. Get all sample data
# 10. Run peak picking of sample data and correct retention time using internal standards
# 11. Align sample data

# Import modules
import json
import os

def lcb_workflow(data_dir):
    pass






def load_internal_standard_library(sample_type, ion_mode):


    if sample_type == "lipidomics":
        if ion_mode == "positive":
            fn = "lipid_istd_pos.json"
        elif ion_mode == "negative":
            fn = "lipid_istd_neg.json"
        
    elif sample_type == "metabolomics":
        if ion_mode == "positive":
            fn = "metabolite_istd_pos.json"
        elif ion_mode == "negative":
            fn = "metabolite_istd_neg.json"
    data_path = os.path.join(os.path.dirname(__file__), 'data', fn)

    with open(data_path, 'r') as f:
        return json.load(f)

