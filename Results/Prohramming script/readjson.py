import csv
import json
from pprint import pprint

class Patient_data:
	pat_count = 0;

	def __init__(self,file_name,patient_id):
		self.file_name = file_name
		self.patient_id = patient_id
		self.tp53_expression = 0
		Patient_data.pat_count += 1

	def display_count(self):
		print(Patient_data.pat_count)

	def display_data(self):
		print('Patient_id:',self.patient_id)
		print('Filename:',self.file_name)
		print('TP53 Expression:',self.tp53_expression)

	def retrieve_expression(self):
		with open('data/'+self.file_name,'r') as f:
			pat_dict = csv.reader(f, delimiter='\t')
			for row in pat_dict:
				if row[0].startswith('ENSG00000141510'):
					self.tp53_expression = row[1]
				else:
					pass
		f.close()

		

# data definition
list_patient = []
tp53_dict = {}


# open metadata and write Patient_data object to list_patient
with open('metadata.json','r') as f_data:
	data = json.load(f_data)



for index in range(len(data)):
	cdata = data[index]
	fn = cdata['file_name']
	if fn.endswith('FPKM-UQ.txt.gz'):
		fn = fn[0:48]
		pat_id = cdata['cases'][0]['submitter_id'][8:12]
		list_patient.append(Patient_data(fn,pat_id))
	else:
		pass

f_data.close()

# retrieve tp53 expression according to matched file in data folder,
# detail see the definition in class
for index in range(len(list_patient)):
	list_patient[index].retrieve_expression()
	tp53_dict[list_patient[index].patient_id] = list_patient[index].tp53_expression
	print('Index:',index,'retrieve successful')

# read from csv file and write to a new one
with open('p1.csv','r') as f_in:
	data = csv.reader(f_in)
	with open('p1_modified.csv','w') as f_out:
		f_out_csv = csv.writer(f_out)
		for row in data:
			if row[0] in tp53_dict:
				row[7] = tp53_dict[row[0]]
			else:
				pass
			f_out_csv.writerow(row)
	
