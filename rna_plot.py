import json
import os
import csv
from pprint import pprint
from xml.etree import ElementTree
from operator import itemgetter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib import pyplot as plt

class Patient_data(object):
	pat_count = 0;

	def __init__(self,patient_id,last_vital_status,overall_survival,disease_free_survival,censored):
		self.patient_id = patient_id
		self.last_vital_status = last_vital_status
		self.overall_survival = overall_survival
		self.disease_free_survival = disease_free_survival
		self.censored = censored
		self.rna_expressions = {}
		self.filename = None
		Patient_data.pat_count += 1

	def display_count(self):
		print(Patient_data.pat_count)

	def display_data(self):
		pprint(self.__dict__)

	def set_filename(self,fn):
		self.filename = fn

	def retrieve_rna_expression(self):
		if self.filename != None:
			with open('RNA-Seq/Exp_Data/'+self.filename,'r') as f:
				expressions = csv.reader(f, delimiter='\t')
				for row in expressions:
					self.rna_expressions[rna_names[row[0][:15]]] = float(row[1])

				return True
		else:
			return False



# data definition

dict_patient = {}

# initialize

with open('data.csv', 'r') as f:
	f_csv = csv.DictReader(f)
	for row in f_csv:
		dict_patient[row['patient_id']] = Patient_data(row['patient_id'],row['last_vital_status'],row['overall_survival'],row['disease_free_survival'],row['censored'])

with open('RNA-Seq/metadata_original.json','r') as f_data:
	data = json.load(f_data)

	for index in range(len(data)):
		cdata = data[index]
		fn = cdata['file_name']
		if fn.endswith('FPKM-UQ.txt.gz'):
			fn = fn[0:48]
			pat_id = cdata['cases'][0]['submitter_id'][8:12]
			if pat_id in dict_patient:
				dict_patient[pat_id].set_filename(fn)

rna_names = {}
i = 0

with open('gene_table.csv','r') as f:
	rows = csv.reader(f)
	for row in rows:
		rna_names[row[0]] = row[1]

i = 0
for key in dict_patient:
	result = dict_patient[key].retrieve_rna_expression()
	if result == True:
		i += 1
		print('Patient ID:',key,'retrieve successful')
	else:
		print('Patient ID:',key,'No data')

print(i,'found match')

expressions = []
for key in dict_patient:
	a = (key,dict_patient[key].retrieve_rna_expression())
	if a[1] != None:
		expressions.append(a)

expressions.sort(key=itemgetter(1))

T_high = []
E_high = []
T_low = []
E_low = []


print('Available partition types:\n0\t1/2\n1\t1/4')
selection = int(input('Choose:'))
if selection == 0:
	range_low = range(int(len(expressions)/2))
	range_high = range(int(len(expressions)/2),len(expressions)) 
	partition_title = '1/2_Partition'
elif selection == 1:
	range_low = range(int(len(expressions)/4))
	range_high = range(3*int(len(expressions)/4),len(expressions))
	partition_title = '1/4_Partition'
else:
	raise NameError('Unexpected partition type')

pprint(range_low)
pprint(range_high)

print('Available plot types:\n0\tOverall_Survival_Curve\n1\tDisease_Free_Survival_Curve')
selection = int(input('Choose:'))
if selection == 0:
	for index in range_low:
		if dict_patient[expressions[index][0]].overall_survival != 'unknown':
			T_low.append(int(dict_patient[expressions[index][0]].overall_survival))
			if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
				E_low.append(1)
			else:
				E_low.append(0)
	for index in range_high:
		if dict_patient[expressions[index][0]].overall_survival != 'unknown':
			T_high.append(int(dict_patient[expressions[index][0]].overall_survival))
			if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
				E_high.append(1)
			else:
				E_high.append(0)
	plot_title = 'Overall_Survival_Curve '
elif selection == 1:
	for index in range_low:
		if dict_patient[expressions[index][0]].disease_free_survival != 'unknown':
			T_low.append(int(dict_patient[expressions[index][0]].disease_free_survival))
			if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
				E_low.append(0)
			else:
				E_low.append(1)
	for index in range_high:
		if dict_patient[expressions[index][0]].disease_free_survival != 'unknown':
			T_high.append(int(dict_patient[expressions[index][0]].disease_free_survival))
			if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
				E_high.append(0)
			else:
				E_high.append(1)
	plot_title = 'Disease_Free_Survival_Curve '
else:
	raise NameError('Unexpected plot type')


results = logrank_test(T_high, T_low, E_high, E_low, alpha=0.95)
results.print_summary()

kmf = KaplanMeierFitter()

kmf.fit(T_high, event_observed=E_high, label='High')

kmf.survival_function_
kmf.median_
print(kmf)
ax = kmf.plot(show_censors=True)




kmf.fit(T_low, event_observed=E_low, label='Low')
kmf.plot(ax=ax,show_censors=True)

plt.title('Kaplan-Meier_'+plot_title+partition_title+' ensg_id:')
plt.show()

