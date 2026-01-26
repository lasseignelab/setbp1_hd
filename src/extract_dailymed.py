import requests
import csv
import time
import os
import argparse
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

    # 2. Define LOINC codes with Fallbacks
    # Key = CSV Column Name : Value = List of LOINC codes to try in order
    sections_to_fetch = {
        "Boxed Warning": ["34066-1"],
        "Pregnancy": ["42228-7", "34077-8", "42229-5"], 
        "Lactation": ["77290-5", "34092-7"],
        "Pharmacokinetics": ["43682-4", "34090-1"] # 34090-1 is Clinical Pharmacology
    }

    drug_info = {"Drug Name": drug_name, "Set ID": set_id}

    for column_name, loinc_list in sections_to_fetch.items():
        section_content = "Section not found."
        
        for loinc in loinc_list:
            url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}/sections.json?section_loinc={loinc}"
            try:
                sec_response = requests.get(url).json()
                if sec_response.get('data'):
                    section_content = clean_html(sec_response['data'][0]['text'])
                    break # Stop once we find data for this section
            except:
                continue
        
        # Special handling for BBB info within the Pharmacokinetics section
        if column_name == "Pharmacokinetics":
            bbb_keywords = ["blood-brain barrier", "crosses the bbb", "cerebrospinal fluid", "cns penetration"]
            found_bbb = "Not explicitly detailed."
            lower_text = section_content.lower()
            for phrase in bbb_keywords:
                if phrase in lower_text:
                    idx = lower_text.find(phrase)
                    found_bbb = section_content[max(0, idx-50) : idx+200].replace('\n', ' ') + "..."
                    break
            drug_info["BBB Information"] = found_bbb
            # We keep the full Pharmacokinetics text in the dict because it's in fieldnames now
            drug_info["Pharmacokinetics"] = section_content
        else:
            drug_info[column_name] = section_content

    return drug_info

def main():
    # --- 1. SET UP THE USER ARGUMENT ---
    parser = argparse.ArgumentParser(description="Extract medication safety data from DailyMed.")
    parser.add_argument("-o", "--output_dir", help="Directory to save the CSV file", default=".")
    args = parser.parse_args()

    # --- 2. CREATE DIRECTORY IF IT DOESN'T EXIST ---
    if not os.path.exists(args.output_dir):
        print(f"Creating directory: {args.output_dir}")
        os.makedirs(args.output_dir)

    # --- 3. JOIN DIRECTORY WITH FILENAME ---
    output_path = os.path.join(args.output_dir, "drug_safety_database.csv")

    drug_list = ["Sertraline", "Amoxicillin", "Lisinopril", "Ibuprofen"]
    fieldnames = ["Drug Name", "Set ID", "Boxed Warning", "Pregnancy", "Lactation", "Pharmacokinetics", "BBB Information"]

    with open(output_path, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for drug in drug_list:
            data = get_drug_data(drug)
            if data:
                writer.writerow(data)
            time.sleep(0.5)

    absolute_path = os.path.abspath(output_path)
    print(f"\nFULL SYSTEM PATH: {absolute_path}")
    print(f"\nYou've got drugs! Data saved to: {output_path}")

if __name__ == "__main__":
    main()
