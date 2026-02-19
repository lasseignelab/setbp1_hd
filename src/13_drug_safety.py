# pylint: disable=import-error
# Script for accessing DailyMed API for safety data--constantly updated so outputs are NOT reproducible, see Zenodo for version used in original analysis # nolint
import requests
import csv
import time
import os
import argparse
from bs4 import BeautifulSoup

# DailyMed base API URL
DAILYMED_BASE = "https://dailymed.nlm.nih.gov/dailymed/services/v2/"


# Helper function for XML text parsing
def clean_xml_text(element):
    if element is None:
        return "Not found"
    return " ".join(element.get_text(separator=" ").split()).strip()


def get_drug_data(drug_name):
    print(f"Fetching data for: {drug_name}...")

    # 1. Search for Set ID (Corrected key: 'setid')
    search_url = f"{DAILYMED_BASE}spls.json?drug_name={drug_name}"
    try:
        response = requests.get(search_url, timeout=20)
        search_data = response.json()

        if not search_data.get("data"):
            print(f"No results found for {drug_name}")
            return None
        set_id = search_data["data"][0]["setid"]
    except Exception as e:
        print(f"Error finding {drug_name}: {e}")
        return None

    # 2. Fetch the Full SPL XML
    xml_url = f"{DAILYMED_BASE}spls/{set_id}.xml"
    try:
        xml_response = requests.get(xml_url, timeout=20)
        if xml_response.status_code != 200:
            print(f"    [!] XML API returned status {xml_response.status_code}")
            return None

        # Parse as XML
        soup = BeautifulSoup(xml_response.content, "xml")
    except Exception as e:
        print(f"    [!] XML Fetch/Parse Error: {e}")
        return None

    drug_info = {
        "Name": drug_name,
        "ID": set_id,
        "BBW": "Not found",
        "Pregnancy": "Not found",
        "Lact": "Not found",
        "Ph kin": "Not found",
        "BBB": "Not explicitly detailed",
    }

    # 3. Find all <section> tags in the SPL
    sections = soup.find_all("section")
    for section in sections:
        # Check the <title> tag inside the section
        title_tag = section.find("title")
        if not title_tag:
            continue

        title_text = title_tag.get_text().upper()

        # Determine which column this section belongs to
        target_key = None
        if "WARNING" in title_text and "BOX" in title_text:
            target_key = "BBW"
        elif "PREGNANCY" in title_text:
            target_key = "Pregnancy"
        elif "LACTATION" in title_text or "NURSING MOTHERS" in title_text:
            target_key = "Lact"
        elif "PHARMACOKINETICS" in title_text or "CLINICAL PHARMACOLOGY" in title_text:
            target_key = "Ph kin"

        if target_key:
            # Extract the <text> block which contains the clinical content
            text_block = section.find("text")
            content = clean_xml_text(text_block)
            drug_info[target_key] = content

            # BBB Logic specifically for the Pharmacology sections
            if target_key == "Pharmacokinetics":
                bbb_keywords = [
                    "blood-brain barrier",
                    "crosses the bbb",
                    "cerebrospinal fluid",
                    "cns penetration",
                ]
                lower_text = content.lower()
                for phrase in bbb_keywords:
                    if phrase in lower_text:
                        idx = lower_text.find(phrase)
                        start_idx = max(0, idx - 50)
                        end_idx = idx + 250
                        drug_info["BBB"] = content[start_idx:end_idx] + "..."

                        break

    return drug_info


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Extract medication safety data from DailyMed."
    )
    parser.add_argument(
        "-i",
        "--input_csv",
        help="Path to CSV file with drug names for querying in second column",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_dir", help="Directory to save the CSV file", default="."
    )
    args = parser.parse_args()

    # Parse input CSV
    drug_list = []
    try:
        with open(args.input_csv, mode="r", encoding="utf-8") as infile:
            reader = csv.reader(infile)
            for row in reader:
                if len(row) >= 2:
                    drug_name = row[1].strip().strip('"')
                    if drug_name:
                        drug_list.append(drug_name)
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        return

    print(f"Loaded {len(drug_list)} drugs for processing.")

    # Create output directory
    if not os.path.exists(args.output_dir):
        print(f"Creating directory: {args.output_dir}")
        os.makedirs(args.output_dir)

    output_path = os.path.join(args.output_dir, "drug_safety_data.csv")

    fieldnames = ["Name", "ID", "BBW", "Pregnancy", "Lact", "Ph kin", "BBB"]

    with open(output_path, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for drug in drug_list:
            data = get_drug_data(drug)
            if data:
                writer.writerow(data)
                print(f"    [+] Successfully processed {drug}")
            time.sleep(1.5)

    absolute_path = os.path.abspath(output_path)
    print(f"\nFULL SYSTEM PATH: {absolute_path}")
    print(f"\nYou've got drugs! Data saved to: {output_path}")


if __name__ == "__main__":
    main()
