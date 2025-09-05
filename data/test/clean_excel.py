import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python clean_excel.py <input_excel> <output_excel>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Step 1: Get original sheet name
xls = pd.ExcelFile(input_file)
original_sheet_name = xls.sheet_names[0]  # take the first sheet

# Step 2: Read Excel without header
df_original = pd.read_excel(input_file, sheet_name=original_sheet_name, header=None)

# Step 3: Remove rows where first column == "tag"
df_clean = df_original[df_original.iloc[:, 0] != "tag"]

# Step 4: Transpose
df_clean = df_clean.T

# Step 5: Replace empty strings with NaN
df_clean = df_clean.replace("", pd.NA)

# Step 6: Copy first row to 4th row where empty
row1 = df_clean.iloc[0]
row4 = df_clean.iloc[3]
df_clean.iloc[3] = row4.combine_first(row1)

# Step 7: Delete rows 1–3
df_clean = df_clean.drop(index=[0, 1, 2]).reset_index(drop=True)

# Step 8: Save both sheets (cleaned as Sheet1, original with same name)
with pd.ExcelWriter(output_file) as writer:
    df_clean.to_excel(writer, index=False, header=False, sheet_name="Sheet1")
    df_original.to_excel(writer, index=False, header=False, sheet_name=original_sheet_name)

print(f"Processed file saved to {output_file} with cleaned data in Sheet1 and original in '{original_sheet_name}'")

