import requests
import csv
import time
from bs4 import BeautifulSoup

def clean_html(html_text):
    if not html_text:
        return "N/A"
    # DailyMed HTML often contains complex tables; this keeps the text readable
    soup = BeautifulSoup(html_text, "html.parser")
    return soup.get_text(separator=' ').strip()

def get_drug_data(drug_name):
    print(f"Fetching data for: {drug_name}...")
    
    # 1. Search for Set ID (Corrected key: 'setid')
    search_url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={drug_name}"
    try:
        response = requests.get(search_url)
        response.raise_for_status()
        search_data = response.json()
        
        if not search_data.get('data'):
            print(f"No results found for {drug_name}")
            return None
            
        # We take the first result. Note: 'setid' has no underscore in the API response.
        set_id = search_data['data'][0]['setid']
    except Exception as e:
        print(f"Error finding {drug_name}: {e}")
        return None

    # 2. Define LOINC codes
    # 34066-1: Boxed Warning
    # 42228-7: Pregnancy (Modern PLLR)
    # 34092-7: Nursing Mothers (Legacy fallback)
    # 77290-5: Lactation (Modern PLLR)
    # 43682-4: Pharmacokinetics
    sections_to_fetch = {
        "Boxed Warning": "34066-1",
        "Pregnancy": "42228-7",
        "Lactation": "77290-5",
        "Pharmacokinetics": "43682-4"
    }

    drug_info = {"Drug Name": drug_name, "Set ID": set_id}

    for section_name, loinc in sections_to_fetch.items():
        url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}/sections.json?section_loinc={loinc}"
        try:
            sec_response = requests.get(url).json()
            
            # If modern section is missing, try legacy fallback for Pregnancy/Lactation
            if not sec_response.get('data'):
                if section_name == "Pregnancy":
                    # Try legacy 'Pregnancy' or 'Teratogenic Effects'
                    url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}/sections.json?section_loinc=34077-8"
                    sec_response = requests.get(url).json()
                elif section_name == "Lactation":
                    # Try legacy 'Nursing Mothers'
                    url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}/sections.json?section_loinc=34092-7"
                    sec_response = requests.get(url).json()

            raw_html = sec_response['data'][0]['text'] if sec_response.get('data') else ""
            clean_text = clean_html(raw_html)
            
            if section_name == "Pharmacokinetics":
                bbb_keywords = ["blood-brain barrier", "crosses the bbb", "cerebrospinal fluid", "cns penetration", "distribution"]
                found_bbb = "Not explicitly detailed in Pharmacokinetics section."
                lower_text = clean_text.lower()
                for phrase in bbb_keywords:
                    if phrase in lower_text:
                        idx = lower_text.find(phrase)
                        # Extract context around the match
                        found_bbb = clean_text[max(0, idx-50) : idx+200].replace('\n', ' ') + "..."
                        break
                drug_info["BBB Information"] = found_bbb
            else:
                drug_info[section_name] = clean_text if clean_text != "N/A" else "Section not found."
        except:
            drug_info[section_name] = "Error retrieving section."

    return drug_info

# --- Execution ---
# Add your list of medications here
drug_list = ["Sertraline", "Amoxicillin", "Lisinopril", "Ibuprofen", "Valproic Acid"]
output_file = "drug_safety_database.csv"

with open(output_file, mode='w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=["Drug Name", "Set ID", "Boxed Warning", "Pregnancy", "Lactation", "BBB Information"])
    writer.writeheader()

    for drug in drug_list:
        data = get_drug_data(drug)
        if data:
            writer.writerow(data)
        time.sleep(0.5) # Politeness delay

print(f"\nProcessing complete. You've got drugs! Saved at {output_file}")

