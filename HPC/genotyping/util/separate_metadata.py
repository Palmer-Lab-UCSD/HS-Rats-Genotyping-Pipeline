import pandas as pd
import sys

origial_metadata = pd.read_csv(sys.argv[1], dtype=str)
metadata_cols = origial_metadata.columns.tolist()
origial_metadata = origial_metadata[origial_metadata['strain'] == 'Heterogenous stock'].reset_index(drop=True)
import_cols = {}
for col in metadata_cols:
	if 'rfid' == col.lower():
		import_cols['Sample_Name'] = col
		continue
	elif 'library_name' == col.lower():
		import_cols['Library_ID'] = col
		continue
	elif 'project_name' == col.lower():
		import_cols['Sample_Project'] = col
		continue
	elif 'barcode' == col.lower():
		import_cols['Sample_Barcode'] = col
		continue
	elif 'runid' == col.lower():
		import_cols['full_run_id'] = col
		continue
	elif 'fastq_files' == col.lower():
		import_cols['Fastq_Files'] = col
		continue
	elif 'pcr_barcode' == col.lower():
		import_cols['PCR_Barcode'] = col
		continue
	elif 'father' == col.lower():
		import_cols['Father'] = col
		continue
	elif 'mother' == col.lower():
		import_cols['Mother'] = col
		continue
	elif 'dames' == col.lower():
		import_cols['Mother'] = col
		continue
	elif 'sires' == col.lower():
		import_cols['Father'] = col
		continue
for new_col in ['Father', 'Mother', 'Fastq_Files', 'Sample_Barcode', 'Sample_Project', 'Library_ID', 'Sample_Name']:
	if new_col not in import_cols.keys():
		origial_metadata.insert(0, new_col, '')
		print("ERROR: " + new_col + " doesn't exist")
	else:
		origial_metadata.insert(0, new_col, origial_metadata[import_cols[new_col]])
temp_sample_ID = origial_metadata[[import_cols["Library_ID"], import_cols["Sample_Name"]]].agg('_'.join, axis=1)
origial_metadata.insert(0, "Sample_ID", temp_sample_ID)
origial_metadata.to_csv(sys.argv[2]+"/sample_sheet.csv", index=False, sep=',')

for unique in origial_metadata['Library_ID'].unique():
	temp_metadata = origial_metadata[origial_metadata['Library_ID'] == unique].reset_index(drop=True)
	if 'PCR_Barcode' not in import_cols.keys():
		pcr_barcode = "NONE"
		print("ERROR: pcr_barcode doesn't exist")
	else:
		pcr_barcode = str(temp_metadata.iloc[0][import_cols['PCR_Barcode']])
	if 'full_run_id' not in import_cols.keys():
		full_run_id = "NONE"
		print("ERROR: runid doesn't exist")
	else:
		full_run_id = str(temp_metadata.iloc[0][import_cols['full_run_id']])
	temp_metadata.to_csv(sys.argv[2]+"/SampleSheet_"+pcr_barcode+"_"+unique+"_"+full_run_id+".csv", index=False, sep=',')
