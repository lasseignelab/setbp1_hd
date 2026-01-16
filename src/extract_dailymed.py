
import requests
import csv
import time
from bs4 import BeautifulSoup

def clean_html(html_text):
    if not html_text:
        return "N/A"
    soup = BeautifulSoup(html_text, "html.parser")
    return soup.get_text(separator=' ').strip()

def get_drug_data(drug_name):
    print(f"Fetching data for: {drug_name}...")
    
    # 1. Search for Set ID
    search_url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={drug_name}"
    try:
        search_data = requests.get(search_url).json()
        if not search_data.get('data'):
            return None
        set_id = search_data['data'][0]['set_id']
    except Exception as e:
        print(f"Error finding {drug_name}: {e}")
        return None

    # 2. Define LOINC codes for sections
    # 34066-1: Boxed Warning
    # 42228-7: Pregnancy
    # 77290-5: Lactation
    # 43682-4: Pharmacokinetics (used for BBB info)
    sections_to_fetch = {
        "Boxed Warning": "34066-1",
        "Pregnancy": "42228-7",
        "Lactation": "77290-5",
        "Pharmacokinetics": "43682-4"
    }

    drug_info = {"Drug Name": drug_name, "Set ID": set_id}

    for section_name, loinc in sections_to_fetch.items():
        url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}/sections.json?section_loinc={loinc}"
        response = requests.get(url).json()
        
        raw_html = response['data'][0]['text'] if response.get('data') else ""
        clean_text = clean_html(raw_html)
        
        if section_name == "Pharmacokinetics":
            # Heuristic search for Blood-Brain Barrier info within the PK section
            bbb_keywords = ["blood-brain barrier", "crosses the BBB", "cerebrospinal fluid", "CNS penetration"]
            found_bbb = "Not explicitly mentioned in PK section."
            for phrase in bbb_keywords:
                if phrase in clean_text.lower():
                    # Take a snippet around the keyword for context
                    idx = clean_text.lower().find(phrase)
                    found_bbb = clean_text[max(0, idx-100) : idx+200] + "..."
                    break
            drug_info["BBB Information"] = found_bbb
        else:
            drug_info[section_name] = clean_text

    return drug_info

# --- Execution ---
drug_list = ["Sertraline", "Amoxicillin", "Lisinopril", "Ibuprofen", "Valproic Acid"]
output_file = "drug_safety_database.csv"

with open(output_file, mode='w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=["Drug Name", "Set ID", "Boxed Warning", "Pregnancy", "Lactation", "BBB Information"])
    writer.writeheader()

    for drug in drug_list:
        data = get_drug_data(drug)
        if data:
            writer.writerow(data)
        time.sleep(0.5) # Politeness delay for API

print(f"Data saved to {output_file}")
