#!/usr/bin/env python3
## example usage: python3 check_logs.py ../logs/

import os

def check_logs(directory, search_text="Gridpack created successfully"):
    """
    Recursively check all .log files in the given directory for a specific text.
    """
    print(f"Scanning log files in: {directory}\n{'-'*50}")
    
    for root, _, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".log"):
                filepath = os.path.join(root, filename)
                try:
                    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                        content = f.read()
                        if search_text in content:
                            print(f"Gridpacks success found in: {filepath}")
                        else:
                            print(f"Gridpack failed in: {filepath}")
                except Exception as e:
                    print(f"Error reading {filepath}: {e}")
    
    print(f"{'-'*50}\nScan complete.")

if __name__ == "__main__":
    import sys
    # Use current directory if none specified
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    check_logs(directory)

