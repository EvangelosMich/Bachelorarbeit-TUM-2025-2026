import xml.etree.ElementTree as ET
import csv
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
xml_file = os.path.join(script_dir, 'hmdb_metabolites.xml') 
output_csv = os.path.join(script_dir, 'hmdb_mini_reference.csv')

def extract_hmdb():
    print("Starting extraction with corrected tag spelling...")
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # We will use the correct spelling for our CSV header for your sanity
        writer.writerow(['accession', 'name', 'monoisotopic_mass', 'formula'])
        
        count = 0
        written_count = 0
        
        for event, elem in ET.iterparse(xml_file, events=('end',)):
            if elem.tag.endswith('metabolite'):
                data = {'accession': '', 'name': '', 'monisotopic_molecular_weight': '', 'chemical_formula': ''}
                
                for child in elem:
                    tag_name = child.tag.split('}')[-1]
                    if tag_name in data:
                        data[tag_name] = child.text if child.text else ''
                
                # Check for the misspelled tag specifically
                if data['accession'] and data['monisotopic_molecular_weight']:
                    writer.writerow([
                        data['accession'], 
                        data['name'], 
                        data['monisotopic_molecular_weight'], 
                        data['chemical_formula']
                    ])
                    written_count += 1
                
                count += 1
                if count % 5000 == 0:
                    print(f"Scanned {count}... Written {written_count} to CSV.")
                    f.flush()

                elem.clear()
                
    print(f"Done! Created {output_csv} with {written_count} entries.")


def discover_tags():
    print("Checking tag names...")
    context = ET.iterparse(xml_file, events=('end',))
    
    for event, elem in context:
        if elem.tag.endswith('metabolite'):
            print("\nFound a metabolite! Here are the tags inside it:")
            tags = []
            for child in elem:
                # Strip namespace
                tag_name = child.tag.split('}')[-1]
                tags.append(tag_name)
            
            print(", ".join(tags))
            # We only need to see one to diagnose the problem
            break    

if __name__ == "__main__":
    extract_hmdb()